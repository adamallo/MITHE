use warnings;
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Sort::Key::Natural;
use Sort::Key::Maker nat_i_sorter => qw(nat integer);

#Config
my $usage="Usage: $0 -i input_file -o output_file [options]\nOptions:\n-------------\n--filterINDELS: discards INDEL variants.\nThis script reads a vcf file generated with platypus, splits multi-SNV lines in independent single-SNV lines and small haplotype variants in independent SNVs\nATTENTION: GL and PL are kept as in the original VCF file for the different multi-SNV lines. All parameters are kept as original for the separation of haplotype variants in SNVs";
my $description_VCF="##Multi-SNV variants separated in different entries using separateMultipleSNVPlatypus.pl\n##Small haplotype variants separated in different SNV entries using separateMultipleSNVPlatypus.pl";
my $description_VCF_INDEL="##Multi-SNV variants separated in different entries using separateMultipleSNVPlatypus.pl\n##Small haplotype variants separated in different SNV entries using separateMultipleSNVPlatypus.pl\n##INDELS eliminated using separateMultipleSNVPlatypus.pl";
my $info_header_line='##INFO=<ID=MultiSNV,Number=4,Type=String,Description="pos, reference, |-separated alternative alleles of the original multi-SNV VCF line, |-separated original genotypes">'."\n".'##INFO=<ID=Haplotype,Number=3,Type=String,Description="pos, ref, alt">';
my $info_message_multisnv="MultiSNV=";
my $info_message_haplotype="Haplotype=";
my %forbidden_infos=map {$_=>1} qw(GL PL);
my %withref_infos=map {$_=>1} qw(AD);
my %genotype_infos=map {$_=>1} qw(GT);

#GETOPT
my ($input_file,$output_file)=("","");
my $delete_indels=0;
my $help;

(! GetOptions(
    'input_file|i=s' => \$input_file,
    'output_file|o=s' => \$output_file,
    'filterINDELS|f=i' => \$delete_indels,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($input_file eq "") || $help) and die $usage;

##Reading input file
open(my $INPUT_FH, $input_file) or die "ERROR opening the input file $input_file\n";
my @content=<$INPUT_FH>;
close($INPUT_FH);

open(my $OUTPUT, ">$output_file") or die "ERROR opening the output file $output_file\n";

print("Working on the file $input_file, to generate $output_file\n");

##Main loop
my $is_content=0;
my $isinfo=0;
my @columns;
my @ialts;
my $nout;
my $i;
my $j;
my @out;
my @in_info;
my $code;
my $skip=0;
my $ref_forbidden_formats;
my $ref_withref_formats;
my $ref_genotype_formats;
my $ref_out_genotypes;
my ($iref,$ialt,$ipos); #Input values = originals
my ($mref,$malt,$mpos); #MultiSNV value = simplified considering just one alt
my ($ref,$alt,$pos); #Final value = SNVs, after splitting haplotypes if needed
my ($refrefs,$refalts,$refposs); #Reference to arrays with the final values
my @outcontent;

foreach my $line (@content)
{
	if($is_content==1)
	{
		chomp($line);
		@columns=split("\t",$line);
		$ialt=$columns[4];
        $iref=$columns[3];
        $ipos=$columns[1];
		@in_info=split(";",$columns[7]);

        if($line =~ m/^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+,[^\t]+\t/) ##Multi-SNV line
		{
			##ADDING THE MAGIC HERE
			@ialts=split(",",$ialt);
			$nout=scalar @ialts;
			
			##Make masks for fields that must be treated in different ways
			$ref_forbidden_formats=get_isformat($columns[8],\%forbidden_infos);
			$ref_withref_formats=get_isformat($columns[8],\%withref_infos);
            $ref_genotype_formats=get_isformat($columns[8],\%genotype_infos);

			for ($i=0; $i<$nout; ++$i)
			{
                $skip=0;

				#Push constants and alt
				@out=();
                $ref_out_genotypes=[];
                ($mref,$malt,$mpos)=simplify_RApair($iref,$ialts[$i],$columns[1]);
                if($delete_indels==1 && length $mref != length $malt)
                {
                    next;
                }
				push(@out,$columns[0],$mpos,$columns[2],$mref,$malt,@columns[5..6]);
                ##Not adding INFO yet, I cannot use push anymore
				##Add format
				$out[8]=$columns[8];
				##Add samples
				for ($j=9; $j< scalar @columns; ++$j)
				{
					$out[$j]=get_processed_sample_fields([split(":",$columns[$j],$j)],$i,":",$ref_forbidden_formats,$ref_withref_formats,$ref_genotype_formats,$ref_out_genotypes);
                    ##If this variant is not in the genotype of this sample
                    if ($out[$j] =~ m/0\/0/)
                    {
                        if ($j==9) #If this is the first sample
                        {
                            $skip=1;
                        }
                        else
                        {
                            $skip=(1 & $skip);
                        }
                    }
                    else
                    {
                        $skip=0;
                    }
                    #If this variant was not in the genotype of any samples, it will not be in the output
				}

                if ($skip == 1)
                {
                    next;
                }

				##Add INFO, now with the original genotypes
                $out[7]=get_field_alt_addMultiSNVinfo(\@in_info,$i,";",$ipos, $iref, \@ialts,$ref_out_genotypes);
                
                unless(length $mref == length $malt && length $mref >1)
                {
                    push(@outcontent,join("\t",@out));
                }
                else
                {
                    ##Change info
                    $out[7]=get_field_alt_addHaplotypeinfo([split(";",$out[7])],$mpos,$mref,$malt);
                    ($refrefs,$refalts,$refposs)=split_haplotypes($mpos,$mref,$malt);
                    #loop, generate iirefs alts and pos and print them
                    for ($j=0; $j< scalar @{$refrefs}; ++$j)
                    {
                        $out[1]=$refposs->[$j];
                        @out[3..4]=($refrefs->[$j],$refalts->[$j]);
                        push(@outcontent,join("\t",@out));
                    }
                    
                }
			}
		}
		else ##Single variant
		{
            
            if($delete_indels==1 && length $iref != length $ialt)
            {
                    next;
            }
            elsif(length $iref == length $ialt && length $iref >1)
            {
				@out=();
                ($refrefs,$refalts,$refposs)=split_haplotypes($ipos,$iref,$ialt);
                $out[0]=$columns[0];
                $out[2]=$columns[2];
                @out[5..6]=@columns[5..6];
                $out[7]=get_field_alt_addHaplotypeinfo(\@in_info,$ipos,$iref,$ialt);
                splice(@columns,0,8);
                push(@out,@columns); #Pushes from 8 onwards

                for($i=0; $i<scalar @{$refrefs}; ++$i)
                {
                    $out[1]=$refposs->[$i];
                    @out[3..4]=($refrefs->[$i],$refalts->[$i]);
                    push(@outcontent,join("\t",@out));
                }
            }
            else
            {
                push(@outcontent,$line);
            }
		}
	}
	elsif($is_content==0 && $line =~m/^#[^#]/) ##Last VCF header line, add a comment and print
	{
		print($OUTPUT join("\n",$delete_indels==1?$description_VCF_INDEL:$description_VCF,$line));
		$is_content=1;
	}
	elsif($isinfo==1 && !($line=~m/^##INFO=/))
	{
		print($OUTPUT join("\n",$info_header_line,$line));
		$isinfo=0;
	}
	else ##VCF Header, just print it
	{
		if($line=~m/^##INFO=/ )	
		{
			$isinfo=1;
		}
		print($OUTPUT $line);
	}
}

my @int_v;

foreach my $outline (nat_i_sorter{@int_v=split("\t",$_);$int_v[0],$int_v[1]} @outcontent)
{
    print($OUTPUT $outline."\n");
}

close($OUTPUT);
exit;

sub get_field_alt
{
	my ($ref_input,$i,$char)=@_;
	my @columnout;
	my @altsthiscol;
	my $code;

	foreach my $helpi (@{$ref_input})
	{
	    if($helpi =~m/,/)
	    {
	        @altsthiscol=split(",",$helpi);
	        $code=$altsthiscol[0];
	        if($code=~s/^([^=]+=).*$/$1/ == 0)
	        {
	            $code="";
	        }
	        $altsthiscol[0]=~s/^[^=]+=//g;
	        push(@columnout,$code.$altsthiscol[$i]);
	    }
	    else
	    {
	        push(@columnout,$helpi);
	    }
	}
	return join($char,@columnout);
}

#Returns the sample fields for the $i genotype in this multi-SNV line
#It uses two references to mask arrays (binary arrays), $ref_skip_i, $ref_withref_i
#ref_skip_i[$i]==1 are not processed (i.e., are copied as they are for all $i)
#ref_withref_i[$i]==1 are processed taking into consideration that the first position
#refers to the REF allele, and the rest to ALT alleles.
#It also writes into an array ref the genotype of this i
sub get_processed_sample_fields
{
	my ($ref_input,$i,$char,$ref_skip_i,$ref_withref_i,$ref_genotype_i,$ref_out_genotypes)=@_;
	my @columnout;
	my @altsthiscol;
    my @genotype;
	my $code;

	#print("DEBUG: @{$ref_input}\n");
	#print("DEBUG: @{$ref_skip_i}\n");
	for (my $j=0; $j<scalar @{$ref_input}; ++$j)
	{
	    if($ref_skip_i->[$j]==0 && $ref_input->[$j] =~m/,/) ##If this should not be skipped and has an array of data
	    {
	        @altsthiscol=split(",",$ref_input->[$j]);
	        $code=$altsthiscol[0]; #Get the code (if present) and first element
            
            if($ref_withref_i->[$j]==1) ##The first element is the reference. We keep it and then add the needed alternative
            {
                shift(@altsthiscol);
                push(@columnout,join(",",$code,$altsthiscol[$i]));
            }
            else
            {
    	        if($code=~s/^([^=]+=).*$/$1/ == 0) ##We get the code (if present, otherwise clean it from previous iteractions)
    	        {
    	            $code="";
    	        }
    	        $altsthiscol[0]=~s/^[^=]+=//g; ##Eliminate the code
	            push(@columnout,$code.$altsthiscol[$i]);
            }
	    }
        elsif($ref_genotype_i->[$j] == 1)
        {
            push(@{$ref_out_genotypes},$ref_input->[$j]);
            @genotype=split("/",$ref_input->[$j]);
            if($genotype[0]==$i+1)
            {
                if($genotype[1]==$i+1)
                {
                    push(@columnout,"1/1");
                }
                else
                {
                    push(@columnout,"1/0");
                }
            }
            elsif($genotype[1]==$i+1)
            {
                push(@columnout,"0/1");
            }
            else
            {
                push(@columnout,"0/0");
            }
        }
	    else
	    {
	        push(@columnout,$ref_input->[$j]);
	    }
	}
	return join($char,@columnout);
}

sub get_isformat
{
	my @output;
	my ($data,$ref_forbidden_infos)=@_;
	my @infos=split(":",$data);

	for (my $i=0; $i< scalar @infos; ++$i)
	{
		$output[$i]=exists $ref_forbidden_infos->{$infos[$i]}?1:0;
	}
	return \@output;
}

##Function that modifies the info column in two ways:
# 1) Adds an additional field, properly sorted (alphabetically) with the structure $info_message.$pos,$ref,alt1|alt2|..|altn,ref1|ref2|..|refn
# 2) If some info fields have several comma-separated values, it splits them and gets only the one in the position indicated by $i
sub get_field_alt_addMultiSNVinfo
{
	my ($ref_input,$i,$char,$pos,$ref,$ref_alts,$ref_out_genotypes)=@_;
	my $info_to_add=$info_message_multisnv.join(",",$pos,$ref,join("|",@{$ref_alts}).",".join("|",@{$ref_out_genotypes}));
	my @columnout;
	my @altsthiscol;
	my $code;
	my $info_added=0;

	foreach my $helpi (@{$ref_input})
	{
		if($info_added==0 && ($info_to_add lt $helpi))
		{
			push(@columnout,$info_to_add);
			$info_added=1;
		}
	    if($helpi =~m/,/)
	    {
	        @altsthiscol=split(",",$helpi);
	        $code=$altsthiscol[0];
	        if($code=~s/^([^=]+=).*$/$1/ == 0)
	        {
	            $code="";
	        }
	        $altsthiscol[0]=~s/^[^=]+=//g;
	        push(@columnout,$code.$altsthiscol[$i]);
	    }
	    else
	    {
	        push(@columnout,$helpi);
	    }
	}
	return join($char,@columnout);
}

##Function that modifies the info column:
# 1) Adds an additional field, properly sorted (alphabetically) with the structure $info_message_haplotype.$pos,$ref,$alt1
sub get_field_alt_addHaplotypeinfo
{
	my ($ref_input,$pos,$ref,$alt)=@_;
	my $info_to_add=$info_message_haplotype.join(",",$pos,$ref,$alt);
	my @columnout;
	my $info_added=0;

	foreach my $helpi (@{$ref_input})
	{
		if($info_added==0 && ($info_to_add lt $helpi))
		{
			push(@columnout,$info_to_add);
			$info_added=1;
		}
	    else
	    {
	        push(@columnout,$helpi);
	    }
	}
	return join(";",@columnout);
}

#Function to simplify ref and alt when both of them are longer than one position long (from having been called as a multisnv)
#We eliminate common ground from the tail. If they are not ok yet, then from the front. It seems the latter never happens in practice, I don't really understan why, but I will leave it just in case
sub simplify_RApair
{
    my ($ref,$alt,$pos)=@_;
    
    if (length $ref ==1 || length $alt ==1)
    {
        return ($ref,$alt,$pos)
    }

    my @Cref=split("",$ref);
    my @Calt=split("",$alt);
    my $lenref=scalar @Cref;
    my $lenalt=scalar @Calt;

    my $minlength=$lenref<$lenalt?$lenref:$lenalt;
    my $ncommon=0;

    for (my $i=1; $i<=$minlength; ++$i)
    {
        if($Cref[$lenref-$i] eq $Calt[$lenalt-$i])
        {
            $ncommon+=1;
        }
        else
        {
            last;
        }
    }
    
    if ($ncommon==$lenref || $ncommon==$lenalt)
    {
        $ncommon-=1;
    }
    
    $lenref-=$ncommon;
    $lenalt-=$ncommon;
    splice(@Cref,$lenref);
    splice(@Calt,$lenalt);

    if($lenref==1 || $lenalt==1)
    {
        return(join("",@Cref),join("",@Calt),$pos);
    }
   
    $minlength=$lenref<$lenalt?$lenref:$lenalt;
    $ncommon=0;

    for (my $i=0; $i<$minlength; ++$i)
    {
        if($Cref[$i] eq $Calt[$i])
        {
            $ncommon+=1;
        }
        else
        {
            last;
        }
    }
    
    if ($ncommon==$lenref || $ncommon==$lenalt)
    {
        $ncommon-=1;
    }
    splice(@Cref,0,$ncommon);
    splice(@Calt,0,$ncommon);
    
    return(join("",@Cref),join("",@Calt),$pos+$ncommon);
}

#Splits the ref, alt, and pos of an haplotype into array references of the individual SNVs that compose it
sub split_haplotypes
{
    my ($pos,$ref,$alt)=@_;
    
    my @Cref=split("",$ref);
    my @Calt=split("",$alt);
    my $length=scalar @Cref;

    my @outref;
    my @outalt;
    my @outpos;
    
    for(my $i=0; $i<$length; ++$i)
    {
        if($Cref[$i] ne $Calt[$i])
        {
            push(@outref,$Cref[$i]);
            push(@outalt,$Calt[$i]);
            push(@outpos,$pos+$i);
        }
    }
    return (\@outref, \@outalt, \@outpos);
}
