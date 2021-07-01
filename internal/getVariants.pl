use warnings;
use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use File::Path;
use File::Basename;
use Sort::Key::Natural;
use Sort::Key::Maker nat_i_sorter => qw(nat integer);
use Bio::DB::HTS::Tabix;

#IO
my $indir="";
my $covB="covB";
my $PAF="PAF";
my $filt_file="";
STDOUT->autoflush(1);

#Configuration variables
my $OFS="\t";
my $IFS="\t";
my $MFS="," ;
my $min_coverage=20;
my $max_alternate_reads=5;
my $nsamples=0;
my $output_sufix="";

#Flags
my $help;

my $usage="Usage: $0 -d directory -n number_of_samples -o output_prefix [options]\nThis script parses the output of MITHE_control.pl from a single patient and generates a list of variants for them taking into consideration the information of all samples. It generates a VCF file per sample, a sample file with information of number of variants in each sample, variants file with the presence of each in the different samples and classification of clonal, subclonal, and private, and a CSV file with the genotype of each sample for each variant. FASTA and CloneFinder files are quick and dirty solutions at this point and need to be improved\n-d output_dir for this patient\n-n number of tumor samples for this patient\n-o sufix of the directory in which this script will generate a vcf file per sample and an additional statistics file (output directory will be samplename_sufix)\n\nOptions:\n--------\n\t--covB: suffix to select the covB results. Default:\"$covB\"\n\t--PAF: suffix to select the PAF results. Default:\"$PAF\"\n\t--min_coverage: minimum number of reads to recall a variant. Default:\"$min_coverage\"Below that threshold, a variant in a sample will be considered ?/?\n\t--min_alternative_reads: minimum number of reads to call a non-called variant present in another sample. Default:\"$max_alternate_reads\"\n\t--filt_file: file with a list of variants to keep, the rest are not taken into consideration. If this file is not used, all of them are taken into consideration\n\nRecall process: if a variant has not been called in the sample is considered ?/? if it does not have enough coverage (min_coverage). If it does have enough coverage, it will be considered 0/0 if the number of alternative reads is n<min_alternative_reads. Otherwise, the variant will be recalled. In the non-recalling output, the latter would be ?/? (never call something that has not been called, unless it is a reference 0/0)\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
    'input_dir|d=s' => \$indir,
    'nsamples|n=i' => \$nsamples,
    'output_sufix|o=s' => \$output_sufix,
    'covB|b=s' =>\$covB,
    'PAF|p=s' =>\$PAF,
    'min_coverage|c=i' =>\$min_coverage,
    'min_alternative_reads|a=i' =>\$max_alternate_reads,
    'filt_file=s' => \$filt_file,
    'help|h' => \$help,
                )) or (($output_sufix eq "") || ($indir eq "") || (! -d $indir) || ($nsamples == 0) || $help) and die $usage;

#Get environment variables
my $GNOMAD=$ENV{'GNOMAD'};

##Main data structures
my @dicts; #Array of dict arrays, one per sample, with information on all their variants
my @fdicts; #Array of dict arrays, one per sample, with information on all the variants filtered (kept) by MITHE.
my @sdicts; #Array of dict arrays, one per sample, with information on all the variants filtered (kept) by MITHE, after syncing without recalling (non-present variants can only be 0/0 or ?/?)
my @rdicts; #Array of dict arrays, one per samle, with information on all the variants filtered (kept) by MITHE, after re-calling considering all other samples
my %vars; #Dict of variants, with information on the number of samples that carry them
my %rvars; #Dict of variants, after re-calling considering all other samples
my %samples; #Dict of samples (SXX format), with refs to arrays with the number of clonal, subclonal, and private mutations
my %rsamples; #Dict of samples (SXX format), with refs to arrays with the number of clonal, subclonal, and private mutations after re-calling considering all other samples
my %nvars; #Dict of variant information in the normal

my @ovcfs; #Array of the original vcf filenames for each sample
my @folders;
my $ndigits;

$ndigits=length $nsamples;
my %varsToKeep;
my $ref_parse_function;

##Parsing filterfile
if ($filt_file ne "")
{
    mopen(my $FH_FILTER, $filt_file);
    %varsToKeep=map{s/\n$//r=>1} <$FH_FILTER>;
    close($FH_FILTER);

    $ref_parse_function=\&parse_selected_vcf;
}
else
{
    $ref_parse_function=\&parse_all_vcf;
}

#Parsing original vcf files
my @samplename;
print("\nParsing original vcf input files from all $nsamples samples in the folder $indir:\n");
#Information to find the files that contain the variants we want
open(my $vcfdictFH, "$indir/vcfdict.csv") or die "ERROR: $indir/vcfdict.csv cannot be opened\n";
my @vcfdictcontent=<$vcfdictFH>;
close($vcfdictFH);

for (my $i=0; $i<$nsamples; ++$i)
{
    my $sample=sprintf("S%0*d",$ndigits,$i);
    
    #Detecting sample, and finding the file to parse and store to use as vcf input for output
    ($ovcfs[$i],$samplename[$i])=get_vcf_samplename(",.*$sample#",\@vcfdictcontent);

    print("\tParsing $sample name $samplename[$i], original vcf file $ovcfs[$i]...");
    
    $dicts[$i]=$ref_parse_function->($ovcfs[$i]);
    
    print("Done\n");
}
print("Done\n");

#Parsing final vcf input files to select among the original samples
print("\nParsing filtered vcf input files from all comparisons of $nsamples samples in the folder $indir...");

#Parsing vcf files
my $name=basename(Cwd::abs_path("$indir"));
my $outdir=join("_",$indir,$output_sufix);

#Initializing fdicts
for (my $i=0; $i<$nsamples; ++$i)
{
   $fdicts[$i]={};
}

##Main parsing loop
for (my $i=0; $i<$nsamples-1; ++$i)
{
    my $Asample=sprintf("S%0*d",$ndigits,$i);
    my $ovcf;
    
    for (my $j=$i+1; $j<$nsamples; ++$j)
    {
        #Identify A and B samples
        my $Bsample=sprintf("S%0*d",$ndigits,$j);
        my $folder=join("_",$name,$Asample,$Bsample);
        push(@folders,$folder);
        
        #Look for the files that containt the A and B variants we want
        open(my $vcfdictFH, "$indir/$folder/vcfdict.csv") or die "ERROR: $indir/$folder/vcfdict.csv cannot be opened\n";
        my @vcfdictcontent=<$vcfdictFH>;
        close($vcfdictFH);

        my @afiles=get_files(\@vcfdictcontent,"Afilt${covB}NAB${PAF}#.*different.vcf");
        my @bfiles=get_files(\@vcfdictcontent,"Bfilt${covB}NAB${PAF}#.*different.vcf");
        my @commonfiles=get_files(\@vcfdictcontent,"^filt${covB}NABU${PAF}#.*common.vcf");

        scalar  @afiles != 1 || scalar @bfiles != 1 || scalar @commonfiles != 1 and die "Error detecting input vcf files in folder $folder\n";
    
        $fdicts[$i]=add_vcfvariants($fdicts[$i],\@afiles); ##POSSIBLE PROBLEM
        $fdicts[$j]=add_vcfvariants($fdicts[$j],\@bfiles); ##POSSIBLE PROBLEM
       
        #I don't have a subroutine for two but I could 
        my $refhash=$ref_parse_function->($commonfiles[0]);
        $fdicts[$i]={%{$fdicts[$i]},%$refhash}; ##POSSIBLE PROBLEM
        $fdicts[$j]={%{$fdicts[$j]},%$refhash}; ##POSSIBLE PROBLEM    
    }
}

##Loop to update information graving it from @dicts
##When parsing fdicts, the value information with read counts is not correct for common variants. Here, we fix this copying that information from @dicts, which is a superset of variants with correct information

for (my $isample=0; $isample<$nsamples; ++$isample)
{
    foreach my $var (keys %{$fdicts[$isample]})
    {
        ! defined $dicts[$isample]->{$var} and die "Variant present in a filtered file is not present in the original file. This should not be possible and shows that there is something wrong with the files or this script";
        $fdicts[$isample]->{$var}=[@{$dicts[$isample]->{$var}}];#Deep copy of the variant information from $dict to $fdicts
    }
}

print("Done\n");

##Generate statistics with original variants
print("\nAnalyzing and summarizing original variants... ");
##Dictionary of variants, with number of samples that bear them
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    foreach my $var (keys %{$fdicts[$isample]})
    {
        $vars{$var}+=1;
    }
    
} #foreach my $sample

##DEBUG
#for (my $isample=0; $isample<$nsamples; ++$isample)
#{
#    foreach my $var (keys %{$fdicts[$isample]})
#    {
#        print("DEBUG: Sample $isample, variant $var, info ".join(" ",@{$dicts[$isample]->{$var}})."\n");
#    }
#    
#} #foreach my $sample

##Dictionary of samples, with the number of private, shared, clonal mutation
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    my $Asample=sprintf("S%0*d",$ndigits,$isample);
    $samples{$Asample}=[0,0,0];
    foreach my $var (keys %{$fdicts[$isample]})
    {
        if($vars{$var}>1)
        {
            if($vars{$var}==$nsamples)
            {
                $samples{$Asample}->[0]+=1;
            }
            else
            {
                $samples{$Asample}->[1]+=1;
            }
        }
        else
        {
            $samples{$Asample}->[2]+=1;
        }
    }
    
} #foreach my $sample
print("Done\n");

print("\nObtaining population allele frequency data for all parsed filtered variants... ");
my %PAF_data=%{getPAFdata(\%vars)};
print("Done\n");

print("\nParsing covB files to get cross-sample basic information... ");

opendir(my $dirhandler,$indir) or die "ERROR opening the directory $indir\n";
my @covB_files=grep {/covB.*\.tsv/} readdir($dirhandler);
closedir($dirhandler);
my @covBdata;

for (my $i=0; $i<$nsamples; ++$i)
{
    my $sample=sprintf("S%0*d",$ndigits,$i);
    my @temp=grep {/^$sample/} @covB_files;
    scalar @temp != 1 and die "Error detecting covBfiles for sample $sample, detected ".scalar @temp." instead of one\n";
    my $covBfile=$temp[0];
    $covBdata[$i]=parse_tsv("$indir/$covBfile"); 
}

print("Done\n");

#Parsing N
print("\nParsing covN file to get control sample information... ");
opendir($dirhandler,$indir) or die "ERROR opening the directory $indir\n";
my @covN_files=grep {/covN.*\.tsv$/} readdir($dirhandler);
closedir($dirhandler);
scalar @covN_files != 1 and die "Error detecting the covN file, detected ".scalar @covN_files." instead of one\n";
my $nfile=join("/",$indir,$covN_files[0]);
%nvars=%{parse_tsv($nfile)};
print("Done\n");

#Re-calling
print("\nRe-calling variants considering multiple samples...\n\tINFO: Minimum number of reads for call: $min_coverage, Minimun number of alternative reads for call: $max_alternate_reads\n");

for (my $isample=0; $isample<$nsamples; ++$isample)
{
    $rdicts[$isample]={};
    $sdicts[$isample]={};
    my $refdictr=$rdicts[$isample];
    my $refdicts=$sdicts[$isample];
    foreach my $var (keys %vars)
    {
        my $realvar=$var;
        
        #Existing filtered variant, just copy it and we are done with it (deep copy)
        if(exists $fdicts[$isample]->{$var})
        {
            $refdictr->{$var}=[@{$fdicts[$isample]->{$var}}];
            $refdicts->{$var}=[@{$fdicts[$isample]->{$var}}];
            next;
        }
        else 
        {
            if (exists $dicts[$isample]->{$var}) ##Variant in original vcf file for this sample, but removed by MITHE
            {
                #Copy the information
                $refdictr->{$var}=[@{$dicts[$isample]->{$var}}];
                $refdicts->{$var}=[@{$dicts[$isample]->{$var}}];
                
                #But change the genotype to undetermined, since it was removed by MITHE. It will be recalled below if needed
                $refdictr->{$var}->[0]="?/?";
                $refdicts->{$var}->[0]="?/?";

                print("\n\tDEBUG: variant $var in sample $isample was discarded by MITHE\n");
            }
            else ##We get the information from covB. It must be there, otherwise there is a problem with the data
            {
                if(! exists $covBdata[$isample]->{$realvar}) #The original variable is not present in covB, so we need to simplify it to get depth information
                {
                    $var=altvar($realvar);
                    exists $covBdata[$isample]->{$var} or die "ERROR: no covB information for $realvar or $var in sample $isample\n"; #Neither variants are present, but they should
                }
                
                #The variant was never called for this sample, so we considered it 0/0 for the non-recalling option. For the re-calling option, we initiate it as ?/?, because we will re-call it if we can.
                #covBdataformat:\@[number of total reads, number of reference reads, number of alternative reads supporting this variant, number of total alternative reads]
                #dictformat: genotype, total reads, reference reads, this alternative reads
                $refdicts->{$realvar}=["0/0",@{$covBdata[$isample]->{$var}}[0..1],0]; 
                $refdictr->{$realvar}=["?/?",@{$covBdata[$isample]->{$var}}[0..1],0];
                
                #If the alternative is the same, we copy the number of alternatives, otherwise it should be 0
                if($var eq $realvar)
                {
                    $refdicts->{$realvar}->[3]=$covBdata[$isample]->{$var}->[2];
                    $refdictr->{$realvar}->[3]=$covBdata[$isample]->{$var}->[2];
                }

                print("\n\tDEBUG: variant $realvar, here as $var, in sample $isample was never present, so we are using covB information and setting it as 0/0 for the regular output and ?/? for the recall output. Variant information: ".join(" ",@{$refdicts->{$realvar}})."\n");
            }
        }
        #There are some problems here
        #Re-calling. Only modifies refdictr
        if($refdicts->{$realvar}->[1]>=$min_coverage) ##Enough information to call the variant
        {
            if($refdicts->{$realvar}->[3]>=$max_alternate_reads) #We recall the variant here for the recalling output, ?/? for the other (we don't do anything). ATTENTION: recalling never generates an homozygous variant call. We may want to modify this if all the reads are variants?
            {
                print("\n\tDEBUG: recalling $var\n");
                $refdictr->{$realvar}->[0]=$covBdata[$isample]->{$var}->[2]>$covBdata[$isample]->{$var}->[1]?"1/0":"0/1";
            }
            elsif($refdicts->{$realvar}->[1]-$refdicts->{$realvar}->[2]<$max_alternate_reads)#Difference between the total number of reads and the number of reference reads, to calculate the number of alternative reads of any kind #There was enough depth and not too many alternate reads, so we call this as reference. WARNING: this may be too stringent if the depth is extreme. It is probably better to use proportion of reads instead of an absolute number. I am using any variant, instead of the variant of interest, because this will work much better with indels that are difficult to call, and would otherwise be called 0/0 here when it is unclear if it is the same indel but with the breakpoints not detected properly or a different one.
            {
                $refdictr->{$realvar}->[0]="0/0";
            }
        }

    }#foreach variant
}

for (my $isample=0; $isample<$nsamples; ++$isample)
{
    foreach my $var (keys %{$rdicts[$isample]})
    {
        if($rdicts[$isample]->{$var}->[0] ne "0/0" && $rdicts[$isample]->{$var}->[0] ne "?/?")
        {
            $rvars{$var}+=1;
        }
    }
    
} #foreach my $sample

##Dictionary of samples, with the number of private, shared, clonal mutation
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    my $Asample=sprintf("S%0*d",$ndigits,$isample);
    $rsamples{$Asample}=[0,0,0];
    foreach my $var (keys %{$rdicts[$isample]})
    {
        if(exists $rvars{$var} && $rdicts[$isample]->{$var}->[0] ne "0/0" && $rdicts[$isample]->{$var}->[0] ne "?/?")
        {
            if($rvars{$var}>1)
            {
                if($rvars{$var}==$nsamples)
                {
                    $rsamples{$Asample}->[0]+=1;
                }
                else
                {
                    $rsamples{$Asample}->[1]+=1;
                }
            }
            else
            {
                $rsamples{$Asample}->[2]+=1;
            }
        }
    }
    
} #foreach my $sample

print("Done\n");

print("\nWritting output files...");
##Outputs
mkpath("$outdir/vcfs");
mkpath("$outdir/stats");
mkpath("$outdir/alignments");

#VCF per file, only with original variants
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    my $sample=sprintf("S%0*d",$ndigits,$isample);
    my $outfile="$outdir/vcfs/${sample}_final.vcf";

    # 1 vcf file per sample
    write_variant_vcf($dicts[$isample],$outfile,$ovcfs[$isample],"#Variants filtered using MITHE");
}

#Samples files, with the number of clonal, subclonal, and private variants
write_sample_file("$outdir/stats/sample.csv",\%samples);
write_sample_file("$outdir/stats/sample_recalled.csv",\%rsamples);

#Variant files, with the number of samples that contain them and a classification of clonal, subclonal, and private
write_variant_file("$outdir/stats/var.csv",\%vars);
write_variant_file("$outdir/stats/var_recalled.csv",\%rvars);

#Comprehensive variant file, similar to ITHE results
write_full_variantlist_file("$outdir/stats/final_variants.csv",\%vars,\@sdicts);
write_full_variantlist_file("$outdir/stats/final_variants_recalled.csv",\%rvars,\@rdicts);

#CSV files with the genotype in each sample for each variant
write_genotype_file("$outdir/stats/genotypes.csv",\@sdicts,\%vars);
write_genotype_file("$outdir/stats/genotypes_recalled.csv",\@rdicts,\%rvars);

#FASTA files TODO:WARNING:WORKING HERE this is a quick and dirty solution. I am only considering the mutated allele, and eliminating INDELS. The normal is a fake all 0/0, when it should be called from covN
write_FASTA_file("$outdir/alignments/alignment.fas",\@sdicts,\%vars);
write_FASTA_file("$outdir/alignments/alignment_recalled.fas",\@rdicts,\%rvars);

#CloneFinder files TODO:WARNING:WORKING HERE this is a quick and dirty solution. I am considerint all variants independently, and liminating INDELS.
write_CloneFinder_file("$outdir/stats/cloneFinder.tsv",\@rdicts,\%rvars,\@samplename); ##In this case, we only want the recalled one, since it has gathered more information. We are outputing all the info on reference, alternate reads, even if the variant has not been called for our purpose. We can use the same filters (minimum number of reads and number of alternate reads) in cloneFinder.

print(" Done\n");

exit;

##FUNCTIONS

sub write_CloneFinder_file
{
    my ($filename,$refdata,$refvars,$refnames)=@_;

    mopen(my $OUT_CloneFinder, ">$filename");
    my @int_v;
    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refvars;
    my (@finalvars,@refs,@alts);
    my @splitvar;
    foreach my $var (@sortedvars)
    {
        @splitvar=split($OFS,$var);
        if(length $splitvar[2] == length $splitvar[3]) #INDEL filter
        {
            push(@finalvars,$var);
            push(@refs,$splitvar[2]);
            push(@alts,$splitvar[3]);
        }
    }
    
    #Building the header
    my @header=("SNVID","Wild","Mut");
    push(@header,map{("$_:ref","$_:alt")} @$refnames);
    print($OUT_CloneFinder join($OFS,@header),"\n");
    my @outcontent;
    my $var;
    for (my $ivar=0; $ivar<scalar @finalvars; $ivar++)
    {
        $var=$finalvars[$ivar];
        @outcontent=("S$ivar",$refs[$ivar],$alts[$ivar]);
        for (my $iname=0; $iname<scalar @$refnames; ++$iname)
        {
            if(defined $refdata->[$iname]->{$var}->[2])
            {
                push(@outcontent,@{$refdata->[$iname]->{$var}}[2,3]);
            }
            else
            {
                push(@outcontent,0,0);
            }
        }
        print($OUT_CloneFinder join($OFS,@outcontent),"\n");
    }
}

sub write_genotype_file
{
    my ($filename,$refdata,$refvars)=@_; 
    mopen(my $OUT_GENOTYPES, ">$filename");

    #Building the multi-line header, with information on the CHROM, POS, REF, and ALT
    my @outcontent=("CHROM","POS","REF","ALT");
    my @int_v;
    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refvars;
    foreach my $var (@sortedvars)
    {
        my @thisdata=split($OFS,$var);
        for (my $idata=0; $idata<scalar @thisdata; ++$idata)
        {
            $outcontent[$idata].=$OFS.$thisdata[$idata];
        }
    }
    my $nheaders=scalar @outcontent;

    #Building the actual lines with genotypes
    for (my $isample=0; $isample<$nsamples; ++$isample)
    {
        $outcontent[$isample+$nheaders]=$samplename[$isample];
        for (my $ivar=0; $ivar<scalar @sortedvars; ++$ivar)
        {
            $outcontent[$isample+$nheaders].=$OFS.$refdata->[$isample]->{$sortedvars[$ivar]}->[0];
        }
    }
    
    foreach my $line (@outcontent)
    {
        print($OUT_GENOTYPES $line."\n");
    }
    
    close($OUT_GENOTYPES);
}

sub write_FASTA_file
{
    my ($filename,$refdata,$refvars)=@_; 
    mopen(my $OUT_MSA, ">$filename");

    my @int_v;
    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refvars;
    my @refs;
    my @alts;
    my @finalvars;
    my @thisdata;

    #Making a list of references and alternates to be used to print later
    #Also filtering out INDELS  
    foreach my $var (@sortedvars)
    {
        @thisdata=split($OFS,$var);
        if(length $thisdata[2] == length $thisdata[3]) #INDEL filter
        {
            push(@finalvars,$var);
            push(@refs,$thisdata[2]);
            push(@alts,$thisdata[3]);
        }
    }

    #Fake normal
    print($OUT_MSA ">Normal\n");
    for (my $ivar=0; $ivar<scalar @finalvars; ++$ivar)
    {
        print($OUT_MSA $refs[$ivar]);
    }
    
    for (my $isample=0; $isample<$nsamples; ++$isample)
    {
        print($OUT_MSA "\n>$samplename[$isample]\n");
        for (my $ivar=0; $ivar<scalar @finalvars; ++$ivar)
        {
            if($refdata->[$isample]->{$finalvars[$ivar]}->[0]=~/1/)
            {
                print($OUT_MSA $alts[$ivar]);
            }
            elsif($refdata->[$isample]->{$finalvars[$ivar]}->[0]=~/\?/)
            {
                print($OUT_MSA "?");
            }
            else
            {
                print($OUT_MSA $refs[$ivar]);
            }
        }
    }
    
    close($OUT_MSA);
}

sub write_variant_file
{
    my ($filename,$refdata)=@_;
    mopen(my $OUT_VAR_STATS, ">$filename");
    print($OUT_VAR_STATS "CHROM,POS,ID,REF,ALT,NSAMPLES,KIND\n");
    my @int_v;
    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refdata;
    foreach my $var (@sortedvars)
    {
        my @temp=split($OFS,$var);
        my $kind="Private";
    
        if($refdata->{$var}>1)
        {
            $kind=$refdata->{$var}==$nsamples?"Clonal":"Subclonal";
        }
        print($OUT_VAR_STATS join(",",@temp[0..1],$var,@temp[2..3],$refdata->{$var},$kind)."\n");
    }
    close($OUT_VAR_STATS);
}

sub write_full_variantlist_file
{
    my ($filename,$refvar,$refdict)=@_;
    mopen(my $OUT_FULLVAR, ">$filename");
    #print($OUT_VAR_STATS "CHROM,POS,ID,REF,ALT,NSAMPLES,KIND\n");
    print($OUT_FULLVAR "CHROM,POS,REF,ALT,SAMPLE,GT,Depth,VariantN,NormalDepth,NormalVariantN,NormalThisVariantN,PAF,INDEL\n");
    my @int_v;
    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refvar;
    my ($genotype,$sample_depth,$sample_variant_reads,$normal_depth,$normal_variant_reads,$normal_this_variant_reads,$paf,$indel);
    my $sample;
    my @keyinfo;
    my @varinfo;
    my $nvar;

    foreach my $var (@sortedvars)
    {
        #Identify if it is an indel
        $indel=isindel($var);

        #Get normal information. In the normal, it may be that the variant does not have the same alternative information, here substituted by a "." and the reference simplified to just one position if it was an indel
        if(! defined $nvars{$var})
        {
            $nvar=altvar($var);
        }
        else
        {
            $nvar=$var;
        }
        
        ! defined $nvars{$nvar} and warn "$var information not found in the normal\n";
        ($normal_depth,$normal_variant_reads,$normal_this_variant_reads)=@{$nvars{$nvar}}[0,3,2];

        #print("DEBUG: original variant: \"$var\", modified variant: \"$nvar\". Information:".join(",",@{$nvars{$nvar}}[2,0,1])."\n");

        #Get PAF information
        ! defined $PAF_data{$var} and warn "$var information not found in PAF\n";
        $paf=@{$PAF_data{$var}}[1];

        @keyinfo=split($OFS,$var);

        #Get Sample-Specific information    
        for (my $isample=0; $isample<$nsamples; ++$isample)
        {
            my $sample=sprintf("S%0*d",$ndigits,$isample);
            
            if(defined $refdict->[$isample]->{$var})
            {
                @varinfo=@{$refdict->[$isample]->{$var}};
                ($genotype,$sample_depth,$sample_variant_reads)=@varinfo[0,1,3];
                print($OUT_FULLVAR join(",",@keyinfo[0..3],$samplename[$isample],$genotype,$sample_depth,$sample_variant_reads,$normal_depth,$normal_variant_reads,$normal_this_variant_reads,$paf,$indel."\n"));
            }
#            else
#            {
#                print("DEBUG: variant $var not present in sample $sample\n");
#            }
        }
    }
    close($OUT_FULLVAR);
}

sub write_sample_file
{
    my ($filename,$refdata)=@_;
    mopen(my $OUT_SAMPLE_STATS, ">$filename");
    print($OUT_SAMPLE_STATS "Sample,Clonal,Subclonal,Private\n");
    for (my $isample=0; $isample<$nsamples; ++$isample)
    {
        my $sample=sprintf("S%0*d",$ndigits,$isample);
        # 1 line per sample summary
        print($OUT_SAMPLE_STATS join(",", $samplename[$isample], @{$refdata->{$sample}})."\n");
    }
    close($OUT_SAMPLE_STATS);
}

sub altvar
{
    my ($var)=@_;
    my $indel=isindel($var);
    if($var=~m/$OFS\.$/)
    {
        return undef
    }
    if ($indel eq 1)
    {
        $var=~s/$OFS([^$OFS])([^$OFS]*)$OFS([^$OFS]*)$/$OFS$1$OFS$3/;
    }
    $var=~s/$OFS[^$OFS]*$/$OFS./;
    return $var;
}

sub get_files
{
    my ($refcontent,$regex,$dir)=@_;
    my @output;
    my @temp;

    if(defined $dir)
    {
        $dir=$dir."/";
    }
    else
    {
        $dir="";
    }

    for my $line (@{$refcontent})
    {
        chomp($line);
        #print("DEBUG: $line\n");
        if ($line =~ m/$regex/)
        {
            @temp=split(",",$line);
            push(@output,$dir.$temp[1]);
           #print("\tDEBUG: pushing $temp[1] in common files\n");
        }
    }
    
    return @output;
}

# Parses the variants in vcf files pointed by the array @{$filesref} and adds them to the hash %{$hashref}
# It calls a parsing function pointed by ref_parse_function
# The value indicates the genotype
# ########################################################################################################

sub add_vcfvariants
{
    my ($hashref,$filesref)=@_;

    foreach my $file (@{$filesref})
    {
        my $irefhash=$ref_parse_function->($file);
        $hashref={%$hashref,%$irefhash};
    }
    return $hashref;

}

# This function uses the parsed info from a vcfdict file to obtain the filename of
# the original vcf file for a sample and the name of that sample, returning them
sub get_vcf_samplename
{
    my ($regex,$refvcfdictcontent,$dir)=@_;
    my @retdata;
    my @temp=get_files($refvcfdictcontent,$regex,$dir);
    scalar @temp != 1 and die "ERROR: detected more than one candidate vcf file for sample $regex\n";
    push(@retdata,$temp[0]);
    @temp=get_sample_names_vcf($retdata[0]);
    scalar @temp != 1 and die "ERROR: detected more than one sample in the vcf file $retdata[0]\n";
    $temp[0]=~s/H_SL-DCIS-//g; ##Simplification to make them easier to read #TODO: HARDCODED
    push(@retdata,$temp[0]);
    return @retdata;
}


##THESE SUBRUTINES SHOULD PROBABLY BE PART OF MY OWN MODULE/LIBRARY

#ATTENTION TODO WARNING: THIS SUBRUTINE IS ERROR PRONE. I AM MODIFYING IT IN A RUSH TO GET SOME DATA QUICKLY
# Parse vcf returning a hash reference with keys in the form CHROM$OFSPOS$OFSREF$OFSALT and value = array_ref [genotype (coded as 0/0, 0/1, 1/0, or 1/1),$nref+$nalt, $nref, $nalt]
##########################################################################
sub parse_selected_vcf
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file cannot be opened";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $value;
    my $filterkey;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            chomp($vcf1[$i]);
            $filterkey=$key=$value=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            $filterkey=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
            if (exists $varsToKeep{$filterkey})
            {
                $value=[(split(":",(split($IFS,$value))[9]))[0,4,5]];
                $value->[3]=$value->[2];
                $value->[2]=$value->[1]-$value->[3];
                $hash{$key}=$value;
            }
        }
    }
    return \%hash;
}

# Parse vcf returning a hash reference with keys in the form CHROM$OFSPOS$OFSREF$OFSALT and value = array_ref [genotype (coded as 0/0, 0/1, 1/0, or 1/1),$nref+$nalt, $nref, $nalt]
##########################################################################
sub parse_all_vcf
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file cannot be opened";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $value;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            chomp($vcf1[$i]);
            $key=$value=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            $value=[(split(":",(split($IFS,$value))[9]))[0,4,5]];
            $value->[3]=$value->[2];
            $value->[2]=$value->[1]-$value->[3];
            $hash{$key}=$value;
        }
    }
    return \%hash;
}

#WARNING the order of information in the value of this dict is different than other version of this function used in other components of ITHE/MITHE

#Parses a TSV file (covN or covB) and returns a dictionary with variant ids and values \@[number of total reads, number of reference reads, number of alternative reads supporting this variant, number of total alternative reads]. Foreach position, at least three entries are added: one with information on the specific variant, another with specific reference but flexible variant ("."), and the last one with only position information. The first is the preferred one, the second is necessary for cases in which 0 variants are found in the normal (or a different variant than in the problem samples), the third is needed in INDELS due to the different annovar format. As noted below, an insertion from A to AC in annovar would be coded as - C. However, UnifyGenotyper would call it as A . if the INDEL is not called. There is no way to go from one format to the other without accessing the reference genome.
#Here I am not using the output of annovar yet, but I will.
sub parse_tsv
{
    my ($file)=@_;
    mopen(my $FT, $file);
    my @content=<$FT>;
    close($FT);
    my %data;
    my @temp;
    my @format;
    my @sample;
    my ($chr, $pos, $ref, $alt, $nref, $nvarreads, $nreads);
    my ($newalt,$newref,$newpos); #IMPORTANT NOTE: Annovar uses a different ref/alt format for INDELS. They never contain repeated information. For example, if ref is A and alt AC, in annovar this will be noted as - C. This generates downstream problems. I am adding a second entry with this format to solve it.
    my @alts;
    my @array_nvarreads;

    foreach my $line (@content)
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
                $data{join($OFS,$chr,$pos,$ref,$alts[$ialt])}=[$nreads,$nreads-$nvarreads,$array_nvarreads[$ialt],$nvarreads];
                $data{join($OFS,$chr,$pos)}=[$nreads,$nreads-$nvarreads,$array_nvarreads[$ialt],$nvarreads];

                if($alts[$ialt]=~m/^$ref/) ##Adding a second entry, as explained in the "important note" right above
                {
                    ($newref,$newalt)=($ref,$alts[$ialt]);
                    $newref="-";
                    $newalt=~s/^\Q$ref//;
                    $data{join($OFS,$chr,$pos,$newref,$newalt)}=[$nreads,$nreads-$nvarreads,$array_nvarreads[$ialt],$nvarreads];
                }
                if($ref=~m/^$alts[$ialt]/) ##Adding a second entry, as explained in the "important note" right above
                {
                    ($newref,$newalt)=($ref,$alts[$ialt]);
                    $newalt="-";
                    $newref=~s/^\Q$alts[$ialt]//;
                    $newpos=$pos+length($alts[$ialt]);
                    #print("DEBUG: before $ref, $alt, $pos. after: $newref, $newalt, $newpos\n");
                    $data{join($OFS,$chr,$pos,$newref,$newalt)}=[$nreads,$nreads-$nvarreads,$array_nvarreads[$ialt],$nvarreads];
                }
            }
        }
    }
    #print("DEBUG:",join(",",keys %data),"\n");
    return \%data;

}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################

sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment)=@_;
    open(my $OFILE, ">$filename") or die "ERROR: The file $filename could not be opened for writting";
    open(my $IFILE, "$vcf");
    my @icontent= <$IFILE>;
    close($IFILE);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;

    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for(my $i=0;$i< scalar @icontent; ++$i)
    {
        if ($flag==0 and $icontent[$i]=~/^##/)
        {
            print($OFILE $icontent[$i]);
        }
        elsif($flag==0)
        {
            print($OFILE "$comment\n$icontent[$i]");
            $flag=1;
        }
        else
        {
            $key=$icontent[$i];
            chomp($key);
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                print($OFILE $icontent[$i]);
                delete($hash{$key});
            }
            if(scalar keys %hash == 0)
            {
                last;
            }

        }

    }
    close $OFILE;
}

# Parses the sample names from a VCF and returns an array with them

sub get_sample_names_vcf
{
    my ($vcfname)=@_;
    my @names;
    
    mopen(my $FT, $vcfname);
    my @content=<$FT>;
    close($FT);
    
    foreach my $line (@content)
    {
        if($line =~m/^#[^#]/)
        {
            chomp($line);
            @names=split($IFS,$line);
            splice(@names,0,9);
            last;
        }
    }
    return @names;
}

# Open with error message
###################################################################################
sub mopen
{
    open($_[0],$_[1]) or die "ERROR: impossible to open the file $_[1] in ".($_[1]=~m/^>/?"write":"read")."mode\n";
    return $_[0];
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

    foreach my $key (keys %{$ref_hash})
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
#             else
#             {
#                 warn "There is no data for $chr:$nstart-".($nstart+1);
#             }
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

    return \%outdata;
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
