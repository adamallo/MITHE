#!/usr/bin/env perl
use strict;
use warnings;
use local::lib;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;
use Bio::DB::HTS::Tabix;

my $usage="Usage: $0 [options] -i input_file -o output_file -f main_folder [-s suffix]\nThe variants file is the 1col result of tabulate_annovar_angelo.pl, filtered for the desired variants\nThis script adds the total number of reads, the number of alternative reads, the number of reads in the normal and the number of alternative reads in the normal, and the population allele frequency for each variant\nThe vcf parsed are in the form of A|B\${suffix}.*[different|common]. It requires the environment variable GNOMAD pointing to the GNOMAD database\n";
######################################################

######################################################
##MAIN
######################################################
##Configuration variables
my $IFS="\t";
my $OFS="\t";
my $FFS=","; #Separator of two cells of the same content

##Input variables
my $inputDir="";
my $output_file="";
my $input_file="";
my $suffix="filtcovBNABPAF";
my $normal='covN.*.tsv$';
my $help;
my $remove;

##Getopt
######################################################
(! GetOptions(
    'input_file|i=s' => \$input_file,
    'output_file|o=s' => \$output_file,
    'main_folder|f=s' => \$inputDir,
    'suffix|s=s' =>\$suffix ,
    'help|h' => \$help,
    'remove|r' => \$remove,
                )) or (($output_file eq "") || ($inputDir eq "")  || (! -s $input_file) || $help) and die $usage;


#Parse input_file
open(my $IFILE, $input_file);
my @variantsContent=<$IFILE>;
close($IFILE);

#Get environment variables
my $GNOMAD=$ENV{'GNOMAD'};

my @row;
my $case;
my $shared;
my $chr;
my $start;
my $end;
my $ref;
my $alt;
my $varid;
my $gene;
my $kind;
my $type;
my %variantdata;
my $header=$variantsContent[0];

print("Parsing input tabulated data...");
for (my $i=1; $i<scalar @variantsContent; ++$i){ ##Skips the header
    @row=split("$IFS",$variantsContent[$i]);
    scalar @row != 10 and die "The line $i of the variant file, $input_file does not conform with the expected format since it has ".scalar @row." columns instead of 10\n";
    ($case,$shared,$chr,$start,$end,$ref,$alt,$gene,$kind,$type)=@row;
    chomp($row[9]);
    $varid="${chr}$OFS${start}$OFS${ref}$OFS${alt}";

    if (exists $variantdata{$case}){
        if(exists $variantdata{$case}{$shared}){
            if(exists $variantdata{$case}{$shared}{$varid}){
                die "Error, several entries for the same variant";
            }else{
                $variantdata{$case}{$shared}{$varid}=join($OFS,@row);
            }
        }
        else{
            $variantdata{$case}{$shared}={$varid => join($OFS,@row)};
        }
    }else{
        $variantdata{$case}={$shared => {$varid => join($OFS,@row)}};
    }

}
print(" Done\n");

print("Obtaining population allele frequency data for parsed variants...");
my $ref_PAF_data=getPAFdata(\%variantdata);
print(" Done\n");

##Main loop, directories
opendir(my $IDIR, $inputDir) or die "can't opendir $inputDir: $!";
my @dirs = grep { /^[^.]/ && -d "$inputDir/$_" } readdir($IDIR);
closedir $IDIR;

print("Writting the results...\n");
#my @shared_states=("Common","A","B");
##Opening the output file and prining the header
open(my $OUTPUT, ">$output_file") or die "Impossible to write the output file\n";
chomp($header);
$header=~s/$IFS/$OFS/g;
print($OUTPUT join("$OFS",$header,"MeanVariantReads", "MeanTotalReads", "AVariantReads","ATotalReads","BVariantReads","BTotalReads", "NormalVariantReads", "NormalThisVariantReads", "NormalTotalReads", "PAF", "INDEL"),"\n");
my ($ref_a,$ref_b,$ref_n);
my $indel;

for (my $idir=0; $idir<scalar @dirs; ++$idir)
{
    $case=$dirs[$idir];
    print("\tCase $case, ".($idir+1)."/".scalar @dirs."\n");
    my @files=get_files($case);

    $ref_a=parseVCF($files[0]);
    $ref_b=parseVCF($files[1]);
    $ref_n=parseTSV($files[2]);

    #DEBUG
#    foreach my $key (keys %{$ref_n})
#    {
#        print("$key @{$ref_n->{$key}}\n");
#    }

    my $varnoalt;
    
    #print("DEBUG: ","case:$case",keys%{$variantdata{$case}->{"A"}});
    foreach my $varid (keys %{$variantdata{$case}->{"A"}})
    {
        $indel=isindel($varid);
        my $nvarid=updatevar($varid);
        
        unless(exists $ref_a->{$varid} && (exists $ref_n->{$nvarid} || $indel==1))
        {
            die "Private A variant $varid/$nvarid not found in the vcf files of the A or N samples for case $case\n";
        }
        
		print($OUTPUT join($OFS,$variantdata{$case}->{"A"}->{$varid},@{$ref_a->{$varid}},@{$ref_a->{$varid}},"NA","NA",exists $ref_n->{$nvarid}?@{$ref_n->{$nvarid}}:["NA","NA","NA"],$ref_PAF_data->{$varid}[1],$indel),"\n");
    }
    foreach my $varid (keys %{$variantdata{$case}->{"B"}})
    {
        $indel=isindel($varid);
        my $nvarid=updatevar($varid);
		unless(exists $ref_b->{$varid} && (exists $ref_n->{$nvarid} || $indel==1))
        {
            die "Private B variant $varid/$nvarid not found in the vcf files of the B or N samples for case $case\n";
        }
		print($OUTPUT join($OFS,$variantdata{$case}->{"B"}->{$varid},@{$ref_b->{$varid}},"NA","NA",@{$ref_b->{$varid}},exists $ref_n->{$nvarid}?@{$ref_n->{$nvarid}}:["NA","NA","NA"],$ref_PAF_data->{$varid}[1],$indel),"\n");
    }
    foreach my $varid (keys %{$variantdata{$case}->{"Common"}})
    {
        $indel=isindel($varid);
        my $nvarid=updatevar($varid);
		unless(exists $ref_a->{$varid} && exists $ref_b->{$varid} && (exists $ref_n->{$nvarid} || $indel==1))
        {
            die "Common variant $varid/$nvarid not found in A, B, or N vcf files for case $case\n";
        }
		print($OUTPUT join($OFS,$variantdata{$case}->{"Common"}->{$varid},($ref_a->{$varid}->[0]+$ref_b->{$varid}->[0])/2,($ref_a->{$varid}->[1]+$ref_b->{$varid}->[1])/2,$ref_a->{$varid}->[0],$ref_a->{$varid}->[1],$ref_b->{$varid}->[0],$ref_b->{$varid}->[1],exists $ref_n->{$nvarid}?@{$ref_n->{$nvarid}}:["NA","NA","NA"],$ref_PAF_data->{$varid}[1],$indel),"\n");
    }
}

close($OUTPUT);
print("Done\n");

if($remove == 1)
{
    print("Removing input file...");
    unlink $input_file;
    print("Done\n");
}

exit;

##Small functions to avoid repeating code
sub updatevar
{
    my ($ivar)=@_;
    if (! exists $ref_n->{$ivar})
    {
        $ivar=~s/$OFS[^$OFS]*$/$OFS\./g;
    }
    if (! exists $ref_n->{$ivar})
    {
        $ivar=~s/$OFS[^$OFS]*$OFS.$//g;
    }
    
    return $ivar;
}

sub isindel
{
    my ($key)=@_;
    my ($chr,$pos,$ref,$alt)=split($OFS,$key);

    $ref =~ s/-//g;
    $alt =~ s/-//g;

    if(length $ref != length $alt)
    {
        return 1;
    }

    return 0;
}

##Larger functions
sub get_files {
    my $directory=$_[0];
    ##We get the info from the ultimate source, the original vcf file.
    opendir(my $DIR, $inputDir."/".$directory) or die "can't opendir $directory: $!";
    my @dircontent = readdir($DIR);
    closedir($DIR);
    my @initial_variant_files;
    my @covN_files;
    my @outfiles;

    foreach my $file (@dircontent)
    {
        if ($file =~ m/^[A|B](#--[^#]+)+.vcf$/)
        {
            push(@initial_variant_files, $file);
        }
        elsif($file =~ m/$normal/)
        {
            push(@covN_files, $file);
        }
    }

    scalar @initial_variant_files != 2 or scalar @covN_files != 1 and die "ERROR: detected ".scalar @initial_variant_files." vcf files containing common variants and ".scalar @covN_files." files with information of coverage in the normal, while expecting 2 and 1, respectively\nFiles:".join(",",@initial_variant_files)."and ".join(",",@covN_files)."\n";

    if($initial_variant_files[0]=~m/^A/)
    {
        $outfiles[0]=$inputDir."/".$directory."/".$initial_variant_files[0];
        $outfiles[1]=$inputDir."/".$directory."/".$initial_variant_files[1];
    }
    else
    {
        $outfiles[0]=$inputDir."/".$directory."/".$initial_variant_files[1];
        $outfiles[1]=$inputDir."/".$directory."/".$initial_variant_files[0];
    }

    $outfiles[2]=$inputDir."/".$directory."/".$covN_files[0];

    foreach my $file (@outfiles)
    {
        if (!-f $file)
        {
            die "The file $file is not valid\n";
        }
    }

    #print("DEBUG: @outfiles\n");
    return @outfiles;

}

sub parseVCF
{
    my $file=$_[0];
    open (my $FILE, $file) or die "Error opening the file $file\n";
    #print("DEBUG: file $file\n");
    my @rawdata=<$FILE>;
    close($file);
    my %data;
    my @temp;
    my @format;
    my @sample;
    my ($chr, $pos, $ref, $alt);
    my ($newalt,$newref,$newpos); #IMPORTANT NOTE: Annovar uses a different ref/alt format for INDELS. They never contain repeated information. For example, if ref is A and alt AC, in annovar this will be noted as - C. This generates downstream problems. I am adding a second entry with this format to solve it.
    my ($nreads,$nvarreads);

    foreach my $line (@rawdata)
    {
        chomp($line);
        if (!($line =~ m/^#/))
        {
            @temp=split("\t", $line);
            ($chr, $pos, $ref, $alt)=@temp[0,1,3,4];
            @format=split(":", $temp[8]);
            @sample=split(":", $temp[9]);
            ($format[(scalar @format) -1] ne "NV") or ($format[(scalar @format) -2] ne "NR") and die "Format problem\n";
            ($nvarreads,$nreads)=@sample[scalar @sample -1, scalar @sample -2];
            $data{"$chr$OFS$pos$OFS$ref$OFS$alt"}=[$nvarreads,$nreads];
            if($alt=~m/^$ref/) ##Adding a second entry, as explained in the "important note" right above
            {
                ($newref,$newalt)=($ref,$alt);
                $newref="-";
                $newalt=~s/^$ref//;
                $data{"$chr$OFS$pos$OFS$newref$OFS$newalt"}=[$nvarreads,$nreads];
            }
            if($ref=~m/^$alt/) ##Adding a second entry, as explained in the "important note" right above
            {
                ($newref,$newalt)=($ref,$alt);
                $newalt="-";
                $newref=~s/^$alt//;
                $newpos=$pos+length($alt);
                #print("DEBUG: before $ref, $alt, $pos. after: $newref, $newalt, $newpos\n");
                $data{"$chr$OFS$newpos$OFS$newref$OFS$newalt"}=[$nvarreads,$nreads];
            }
            #print("DEBUG parseVCF: $chr$OFS$pos$OFS$ref$OFS$alt: $nvarreads,$nreads\n");
        }
    }
    #print("DEBUG:",join(",",keys %data),"\n");   
    return \%data;
}

#Parses a covN tsv file and returns a dictionary with variant ids and values \@[normalvariantreads,normalthisvariantreads,normaltotalreads]. Foreach position, at least thre entries are added, one with information on the specific variant, another with specific reference but flexible variant ("."), and the last with only position information. The first is the preferred one, the second is necessary for cases in which 0 variants are found in the normal (or a different variant than in the problem samples), the third is needed in INDELS due to the different annovar format. As noted below, an insertion from A to AC in annovar would be coded as - C. However, UnifyGenotyper would call it as A . if the INDEL is not called. There is no way to go from one format to the other without accessing the reference genome.
sub parseTSV
{
    my $file=$_[0];
    open (my $FILE, $file) or die "Error opening the file $file\n";
    #print("DEBUG: file $file\n");
    my @rawdata=<$FILE>;
    close($file);
    my %data;
    my @temp;
    my @format;
    my @sample;
    my ($chr, $pos, $ref, $alt, $nref, $nvarreads, $nreads);
    my ($newalt,$newref,$newpos); #IMPORTANT NOTE: Annovar uses a different ref/alt format for INDELS. They never contain repeated information. For example, if ref is A and alt AC, in annovar this will be noted as - C. This generates downstream problems. I am adding a second entry with this format to solve it.
    my @alts;
    my @array_nvarreads;

    foreach my $line (@rawdata)
    {
        chomp($line);
        if (!($line =~ m/^#/))
        {
            @temp=split("\t", $line);
            ($chr, $pos, $ref, $alt, $nref, $nvarreads)=@temp;

            @alts=split(",",$alt);
            @array_nvarreads=split(",",$nvarreads);
            ##We add a indetermined variant with 0 alternatives. If this was already present in the call, it will be eliminated since we are later using a hash, otherwise, it will provide info if the alternative in the problems is different than in the normal
            push(@alts,".");
            push(@array_nvarreads,0);

            $nvarreads=0;

            foreach my $thisvarreads (@array_nvarreads)
            {
                $nvarreads+=$thisvarreads;
            }

            $nreads = $nref+$nvarreads;

            for (my $ialt=0; $ialt<scalar @alts; ++$ialt)
            {
                $data{join($OFS,$chr,$pos,$ref,$alts[$ialt])}=[$nvarreads,$array_nvarreads[$ialt],$nreads];
                $data{join($OFS,$chr,$pos)}=[$nvarreads,$array_nvarreads[$ialt],$nreads];

                if($alts[$ialt]=~m/^$ref/) ##Adding a second entry, as explained in the "important note" right above
                {
                    ($newref,$newalt)=($ref,$alts[$ialt]);
                    $newref="-";
                    $newalt=~s/^\Q$ref//;
                    $data{join($OFS,$chr,$pos,$newref,$newalt)}=[$nvarreads,$array_nvarreads[$ialt],$nreads];
                }
                if($ref=~m/^$alts[$ialt]/) ##Adding a second entry, as explained in the "important note" right above
                {
                    ($newref,$newalt)=($ref,$alts[$ialt]);
                    $newalt="-";
                    $newref=~s/^\Q$alts[$ialt]//;
                    $newpos=$pos+length($alts[$ialt]);
                    #print("DEBUG: before $ref, $alt, $pos. after: $newref, $newalt, $newpos\n");
                    $data{join($OFS,$chr,$pos,$newref,$newalt)}=[$nvarreads,$array_nvarreads[$ialt],$nreads];
                }
            }
        }
    }
    #print("DEBUG:",join(",",keys %data),"\n");   
    return \%data;
    
}
 
##Returns a hash ref to the population allele frequency information for all variants present in the input reference hashes 
#The tabix query gets the data we are looking for and also some neighbors. To keep the output data clean from those, we handle two different hashes.
sub getPAFdata
{
    my $tabix = Bio::DB::HTS::Tabix->new( filename =>$GNOMAD );
    my %outdata; #key: CHROM${OFS}POS${OFS}REF${OFS}ALT value: [$tfilt,$taf];
    my $parsedPAFS;
    my %obtainedPAFs; #key = CHROM${OFS}POS${OFS}REF${OFS}ALT #value = [$tfilt, $taf]
    my $tabix_iter;
    my ($chr, $nstart, $ref, $alt);
    my $line;
    my $i;
    my ($ref_hash)=@_;
    my ($tstart,$talt,$tref,$tfilt,$taf);

    foreach my $case (keys %{$ref_hash})
    {
        foreach my $shared (keys %{$ref_hash->{$case}})
        {
            foreach my $key (keys %{$ref_hash->{$case}->{$shared}})
            {
                #key = CHROM${OFS}POS${OFS}REF${OFS}ALT
                unless (exists $obtainedPAFs{$key})
                {
                    ##Obtain gnomAD data for this genomic position if it has not been obtained before. MULI-SNVs are represented in different lines
                    ($chr,$nstart,$ref,$alt)=split($OFS,$key);
                    $tabix_iter=$tabix->query("$chr:$nstart-".($nstart+1));
                    if(defined $tabix_iter)
                    {
                        while($line=$tabix_iter->next)
                        {
                            #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  S1 ... SN
                            $line =~ s/^[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t[^\t]*AF=([^;]+).*$/$1\t$3\t$2\t$4\t$5/;
                            ($tstart,$talt,$tref,$tfilt,$taf)=split("\t",$line);
                            $obtainedPAFs{"$chr$OFS$tstart$OFS$tref$OFS$talt"}=[$tfilt,$taf];
                        }
        
                    }
    #DEBUG
    #                else
    #                {
    #                    warn "There is no data for $chr:$nstart-".($nstart+1);
    #                }
                }
                
                unless (exists $outdata{$key}) ##If it exists we don't need to do anything, otherwise, we copy the data if we have it, or add NA data since we know we have tried to obtain it and it is not available
                {
                    if(exists $obtainedPAFs{$key})
                    {
                        $outdata{$key}=$obtainedPAFs{$key};
                    }
                    else
                    {
                        $outdata{$key}=["NA","NA"];
                    }
                }

            }#foreach varid
        }#foreach shared
    }#foreach case

    return \%outdata;
}

sub simplify_varid
{
    my ($id)=@_;
    my ($chr,$pos,$ref,$alt)=split($OFS,$varid);
    my @aref=split("",$ref);
    my @aalt=split("",$alt);
    my $nshared=0;
    
    for (my $i=0; $i< scalar @aref; ++$i)
    {
        if($aref[$i] ne $aalt[$i])
        {
            last;
        }
        else
        {
            $nshared+=1;
        }
    }    
    
    my $newref;
    my $newalt;
    my $newpos;

    if($nshared>0)
    {
        $newref=join("",splice(@aref,$nshared));
        $newalt=join("",splice(@aalt,$nshared));
        $newpos=$pos+$nshared;
    }
    else
    {
        $newref=$ref;
        $newalt=$alt;
        $newpos=$pos;
    }

    if($newref eq "")
    {
        $newref="-";
    }
    
    if($newalt eq "")
    {
        $newalt="-";
    }
    
    return(join($OFS,$chr,$newpos,$newref,$newalt));
}
