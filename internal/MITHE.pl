#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;
use Sort::Key::Natural;
use Sort::Key::Maker nat_i_sorter => qw(nat integer);
use Bio::DB::HTS::Tabix;

##ATTENTION: Implementation notes:
#---------------------------------
#Variants are represented in two ways across this program, genotype-specific: CHROM$OFSPOS$OFSREF$OFSALT and location-specific: CHROM$OFSPOS
#Variants read from a VCF file are parsed in the genotype-specific format, while those parsed from a TSV file in the location-specific one (but all information is retained in the value). 
#
#The function get_todelete_variants_covfiltering generates (and returns) a hash of variants to eliminate based on filtering options. These are a mixture of genotype-specific and location-specific variants (see below why).
#Two functions (vcf_prune_tsv_vars and vcf_prune_covB) use this mixture of variants to eliminate VCF variants. This is done as follows:
#    - For each vcf variant, a location-specific variant is calculated and checked against the list of variants to eliminate. Eliminated if found
#    - If the location-specific variant does not exist, check for the genotype-specific
# The variants to delete are a mixture of genotype-specific and location-specific variants because:
#   - Depth filter is location-specific. In these cases the alternative is often reference, and would not be matched by the genotype-specific variant.
#   - Filter of alternative variants is variant-specific, we don't want to eliminate variants because there is presence of another alternative in that position (i.e., germline heterozygous)

##Configuration variables
######################################################

my $sep_param="#";
my $sep_value=":";
my $OFS=",";
my $FS=",";
my $output_vcf=1;
my $output_list=1;
my $output_comprehensive=1;
my $n_cores=1;
my $MAX_OPEN_IT=10;
my $MAX_SYSTEM_IT=10;

######################################################

##Global Variables
######################################################
my %dictrealnamelist;
my @listnames;
my %dictrealnamevcf;
my @vcfnames;
my @namedvcffiles;
my @outfiles;

##IO Variables
######################################################
my $nameA="";
my $nameB="";
my $outputdir="";
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile="";
my $covBfiltercond_inputfile="";
my $popAFfiltercond_inputfile="";
my $output_file="";
my $MITHE_INT=$ENV{'MITHE_INT'};
my $GNOMAD=$ENV{'MITHE_GNOMAD'};

#Flags
my $help;
my $usage="Usage: $0 [options] -a nameA -b nameB -o output_file \n\nWARNING: This is a secondary script that is not inteded to be executed directly.\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the control sample\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--popAF_cond_inputfile: input file for the filter of population allele frequencies using gnomAD\n\t--n_cores : number of cores to execute some steps in parallel (requires the perl package Parallel::Loops)\n\t--output_vcf: (bool) generate resulting vcf files or not\n\t--output_list: (bool) generate resulting list of variants or not\n\t--comp: (int) indicating the comprehensiveness of the output, 0=no files, 1=only needed files to call variants, 2= all intermediate variants\n\n";
######################################################

######################################################
##MAIN
######################################################
#TODO: check if there is a dictionary with files, read it and rename all files to number.vcf_bkp. Then, when I am going to filter I check if there is the need if the file is there with the _bkp, otherwise overwrite. If there is no need to recalculate, mv to the new name. I may want to rewrite the dict file every time I write a vcf file, to make sure that if I re-start a lost job everything is working fine
##Getopt
######################################################
(! GetOptions(
    'fileA|a=s' => \$nameA,
    'fileB|b=s' => \$nameB,
    'exec_cond_inputfile|e=s' => \$execond_inputfile,
	'filt_cond_inputfile|f=s' => \$filtercond_inputfile,
	'NABfilt_cond_inputfile=s' => \$NABfiltercond_inputfile,
    'covaltB_cond_inputfile=s' => \$covBfiltercond_inputfile,
    'popAF_cond_inputfile=s' => \$popAFfiltercond_inputfile,
	'output_file|o=s' => \$output_file,
    'n_cores=i' => \$n_cores,
    'output_vcf=i'=>\$output_vcf,
    'output_list=i' => \$output_list,
    'comp=i' => \$output_comprehensive,
    'help|h' => \$help,
                )) or (($output_file eq "") || $help) and die $usage;

##Load Parallel::Loops if it is available and it's needed
#########################################################

if($n_cores>1)
{
    eval "use Parallel::Loops";
    if($@)
    {
        print "\n\nWARNING: You are asking to execute this script using $n_cores cores, but the required module \"Parallel::Loops\" has not been found in \@INC\n\n";
         $n_cores=1;
   }
    else
    {
        print "\nUsing Parallel::Loops with $n_cores cores\n\n";
    }
}

##Input file parsing and directory creation
######################################################

my @exe_parameters=("input");
my @filtering_parameters=("");
my @NABfiltering_parameters=("");
my @covBfiltering_parameters=("");
my @popAFfiltering_parameters=("");

my @exe_param_values=([("")]);
my @filtering_param_values=([("")]);
my @NABfiltering_param_values=([("")]);
my @covBfiltering_param_values=([("")]);
my @popAFfiltering_param_values=([("")]);

##Input files

if ($execond_inputfile ne "")
{
    @exe_parameters=();
    parse_parameters_values($execond_inputfile,\@exe_parameters,\@exe_param_values);
}

if ($filtercond_inputfile ne "")
{
    @filtering_parameters=();
    parse_parameters_values($filtercond_inputfile,\@filtering_parameters,\@filtering_param_values);
}

if ($NABfiltercond_inputfile ne "")
{
    @NABfiltering_parameters=();
    parse_parameters_values($NABfiltercond_inputfile,\@NABfiltering_parameters,\@NABfiltering_param_values);
}

if ($covBfiltercond_inputfile ne "")
{
    @covBfiltering_parameters=();
    parse_parameters_values($covBfiltercond_inputfile,\@covBfiltering_parameters,\@covBfiltering_param_values);
}

if ($popAFfiltercond_inputfile ne "")
{
    @popAFfiltering_parameters=();
    parse_parameters_values($popAFfiltercond_inputfile,\@popAFfiltering_parameters,\@popAFfiltering_param_values);
}

my $vcf_filt_exe="$MITHE_INT/vcf_filtering.sh";

if (! -f "$vcf_filt_exe")
{
    die "The executable vcf_filtering.sh can't be located. Please, make sure that the environment variable MITHE_INT indicates the directory with this package of scripts\n";
}

my $gnomAD_prec_file;
my $GNOMAD_MIN;
if (! defined "$GNOMAD" || ! -s "$GNOMAD" || ! -s "$GNOMAD.tbi")
{
    die "The gnomAD data can't be found. Make sure you hae configured the MITHE_GNOMAD configuration variale properly. Files: $GNOMAD, and $GNOMAD.tbi\n"
}
else
{
    $gnomAD_prec_file=dirname($GNOMAD)."/n_genomes.info";
    $GNOMAD_MIN=0;
    if (! -s $gnomAD_prec_file)
    {
        die "The gnomAD data can't be found at $gnomAD_prec_file\n";
    }
    else
    {
        open(my $gnomAD_prec_filehandle, $gnomAD_prec_file) or die "Problems opening $gnomAD_prec_file\n";
        my @tempcontent=<$gnomAD_prec_filehandle>;
        $GNOMAD_MIN=1.0/($tempcontent[0]+1);
        close($gnomAD_prec_filehandle);
    }
}

#Parsing inputs generated by MITHE_control and shared among several comparisons
%dictrealnamevcf=parse_vcfdict("vcfdict.csv");
refine_vcfdict(\%dictrealnamevcf);

#Changing pwd and output information accordingly
$outputdir=basename($output_file,(".csv"));
open(my $tempfh,">$output_file");
close($tempfh);
$output_file=Cwd::abs_path("$output_file");
$outputdir=~s/\.[^\.]+$//;
mkdir $outputdir;
$outputdir=Cwd::abs_path("$outputdir");
my $wd=Cwd::abs_path(".");
chdir $outputdir or die "The output directory $outputdir is not accesible";

##The working directory is changed here to the dir where the script will write files
### $wd contains the previous directory, to which relative paths relate

## Parsing static variants
############################################################################################

my %N;
my ($ref_n,$nameN);

if ( -f "$wd/N.vcf")
{
    ($ref_n,$nameN)=parse_vcf_name("$wd/N.vcf");
    #print("DEBUG: N genotype name $nameN\n");
    %N=%{$ref_n};
}
else
{
    die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
}

## Main conditions loop
#######################

my $name_condition;
my @filtering_conditions;
my @NABfiltering_conditions;
my @covBfiltering_conditions;
my @popAFfiltering_conditions;
my @exe_conditions;
#my %results;
my $bamfiles;
my $parallel;
my %constExeCondContentRefs;
my %nofiltResultsRef;

##Generate output information
my $sample=$output_file;
$sample=basename($sample);
$sample=~s/([^.]*)\..*/$1/;

##Generate all combinations of proposed values for execution and filtering parameters
print("Generating all combinations of filtering values to be explored by this run...");
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);
combs(0,"",\@NABfiltering_parameters,\@NABfiltering_param_values,\@NABfiltering_conditions);
combs(0,"",\@covBfiltering_parameters,\@covBfiltering_param_values,\@covBfiltering_conditions);
combs(0,"",\@popAFfiltering_parameters,\@popAFfiltering_param_values,\@popAFfiltering_conditions);
print(" Done\n");

##Generate headers of statistics that change in number and name depending on the input filtering conditions
my @HeaderStatsfiltNPAF,
my @HeaderStatsfiltNcovBPAF;
my @HeaderStatsfiltNABPAF;
my @HeaderStatsfiltNABcovBPAF;
my @HeaderStatsfiltNABcovBPAFPAF;
my @filtNPAFnames=qw(AfiltN_Priv_PAF_ BfiltN_Priv_PAF_ filtN_U_PAF_ );
my @filtNcovBPAFnames=qw(AfiltNcovB_Priv_PAF_ BfiltNcovB_Priv_PAF_ filtNcovB_U_PAF_);
my @filtNABPAFnames=qw(AfiltNAB_Priv_PAF_ BfiltNAB_Priv_PAF_ filtNAB_U_PAF_);
my @filtNABcovBPAFnames=qw(AfiltNABcovB_Priv_PAF_ BfiltNABcovB_Priv_PAF_ filtNABcovB_U_PAF_);
my @filtNABcovBPAFPAFnames=qw(AfiltNABcovBPAF_Priv_PAF_ BfiltNABcovBPAF_Priv_PAF_ filtNABcovBPAF_U_PAF_);

@HeaderStatsfiltNPAF=makeHeaderPAF(\@filtNPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNcovBPAF=makeHeaderPAF(\@filtNcovBPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNABPAF=makeHeaderPAF(\@filtNABPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNABcovBPAF=makeHeaderPAF(\@filtNABcovBPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNABcovBPAFPAF=makeHeaderPAF(\@filtNABcovBPAFPAFnames,\@popAFfiltering_conditions);

#print("DEBUG: @exe_parameters, @exe_param_values");
#print("DEBUG: @exe_conditions,@filtering_conditions,@NABfiltering_conditions\n");

if ($n_cores>1)
{
    $parallel = Parallel::Loops->new($n_cores);
    $parallel->share(\%dictrealnamelist, @listnames, \%dictrealnamevcf, @vcfnames, \@namedvcffiles, \%constExeCondContentRefs, \%nofiltResultsRef, \@outfiles);
}

#Parsing covX files and population allele frequency data for variants of each exe_condition
#TODO:I am not 100% sure this is worth it since the forking process takes quite a lot of time to duble-check if I have time
print("Parsing input data...");
if($n_cores>1 && scalar @exe_conditions > 1)
{
    $parallel->foreach(\@exe_conditions,\&parse_const_execond);
}
else
{
    foreach my $exe_condition (@exe_conditions)
    {
        parse_const_execond($exe_condition);
    }
}
print(" Done\n");

my ($Aexecondname,$Bexecondname);

foreach my $exe_condition (@exe_conditions) ##Options that require to call variants again 
{
    $Aexecondname="$nameA$sep_param$exe_condition.vcf"; ##So far I am not using numbers for execondnames. If I change this, I will need to get the number names here and in the filter function
    $Bexecondname="$nameB$sep_param$exe_condition.vcf";

#print("DEBUG: Exe condition loop $exe_condition\n");
    if ((!-f "$wd/$Aexecondname") || (!-f "$wd/$Bexecondname")) ##VCF file we will use for filtering, if it does not exist we have to perform the variant calling
    {
        die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the two scripts are using the same filtering options\n";
    }

    my @current_conditions;
    for (my $i=0; $i<scalar @filtering_conditions;++$i)
    {
        $current_conditions[$i]=[$exe_condition,$filtering_conditions[$i]];
    }
#    ##DEBUG
#    foreach my $cond (@current_conditions)
#    {
#        print("DEBUG: ", join(" ",@{$cond}, "\n"));
#    }

    if($n_cores>1 && scalar @current_conditions > 1)
    {
        $parallel->foreach(\@current_conditions,\&filter);
    }
    else
    {
        foreach my $filtering_condition (@current_conditions)
        {
            filter($filtering_condition);
        }
    }
}

## Output
#############################################################
open(my $OFILE,">$output_file") or die "ERROR opening the output file $output_file";

##Making a block here to use the no warnings qw. It complains because of the # but they are valid
{
    no warnings 'qw';
    print($OFILE join(",",qw(Sample Condition A_# B_# N_# AN_# BN_# Afilt_prop Afilt_N Afilt_# Bfilt_prop Bfilt_N Bfilt_# filt_propU filt_NU filt_#U filt_propI filt_NI filt_#I filt_prop_mean filt_N_mean filt_#_mean AfiltN_prop AfiltN_N AfiltN_# BfiltN_prop BfiltN_N BfiltN_# filtN_propU filtN_NU filtN_#U filtN_propI filtN_NI filtN_#I filtN_prop_mean filtN_N_mean filtN_#_mean),@HeaderStatsfiltNPAF,qw(AfiltNcovB_prop AfiltNcovB_N AfiltNcovB_# BfiltNcovB_prop BfiltNcovB_N BfiltNcovB_# filtNcovB_propU filtNcovB_NU filtNcovB_#U filtNcovB_propI filtNcovB_NI filtNcovB_#I filtNcovB_prop_mean filtNcovB_N_mean filtNcovB_#_mean),@HeaderStatsfiltNcovBPAF,qw(AfiltNAB_prop AfiltNAB_N AfiltNAB_# BfiltNAB_prop BfiltNAB_N BfiltNAB_# filtNAB_propU filtNAB_NU filtNAB_#U filtNAB_propI filtNAB_NI filtNAB_#I filtNAB_prop_mean filtNAB_N_mean filtNAB_#_mean),@HeaderStatsfiltNABPAF,qw(AfiltNABcovB_prop AfiltNABcovB_N AfiltNABcovB_# BfiltNABcovB_prop BfiltNABcovB_N BfiltNABcovB_# filtNABcovB_propU filtNABcovB_NU filtNABcovB_#U filtNABcovB_propI filtNABcovB_NI filtNABcovB_#I filtNABcovB_prop_mean filtNABcovB_N_mean filtNABcovB_#_mean),@HeaderStatsfiltNABcovBPAF,qw(AfiltNABcovBPAF_prop AfiltNABcovBPAF_N AfiltNABcovBPAF_# BfiltNABcovBPAF_prop BfiltNABcovBPAF_N BfiltNABcovBPAF_# filtNABcovBPAF_propU filtNABcovBPAF_NU filtNABcovBPAF_#U filtNABcovBPAF_propI filtNABcovBPAF_NI filtNABcovBPAF_#I filtNABcovBPAF_prop_mean filtNABcovBPAF_N_mean filtNABcovBPAF_#_mean),@HeaderStatsfiltNABcovBPAFPAF),"\n");
}
    #Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,(nPAFfreqs) x (AfiltN_Priv_PAF_XX,BfiltN_Priv_PAF_XX,filtN_U_PAF_XX),AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,(nPAFfreqs) x (AfiltNcovB_Priv_PAF_XX,BfiltNcovB_Priv_PAF_XX,filtNcovB_U_PAF_XX),AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,(nPAFfreqs) x (AfiltNAB_Priv_PAF_XX,BfiltNAB_Priv_PAF_XX,filtNAB_U_PAF_XX),AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean,(nPAFfreqs) x (AfiltNABcovB_Priv_PAF_XX,BfiltNABcovB_Priv_PAF_XX,filtNABcovB_U_PAF_XX),AfiltNABcovBPAF_prop,AfiltNABcovBPAF_N,AfiltNABcovBPAF_#,BfiltNABcovBPAF_prop,BfiltNABcovBPAF_N,BfiltNABcovBPAF_#,filtNABcovBPAF_propU,filtNABcovBPAF_NU,filtNABcovBPAF_#U,filtNABcovBPAF_propI,filtNABcovBPAF_NI,filtNABcovBPAF_#I,filtNABcovBPAF_prop_mean,filtNABcovBPAF_N_mean,filtNABcovBPAF_#_mean,(nPAFfreqs) x (AfiltNABcovBPAF_Priv_PAF_XX,BfiltNABcovBPAF_Priv_PAF_XX,filtNABcovBPAF_U_PAF_XX)

close($OFILE);

## Concatenate and eliminate temporary output files
my $OFH;
myopen(\$OFH,">>$output_file") or die "Error concatenating the temporary output files to $output_file\n";

foreach my $tempOfilename (@outfiles)
{
    my $FH;
    myopen(\$FH,$tempOfilename) or die "Error opening the temporary output file $tempOfilename\n";
    local $/=\65536; # read 64kb chunks
    while ( my $chunk = <$FH>) #Avoids to slurp the whole content to save memory
    {
        print $OFH $chunk;
    }
    close($FH);
    unlink($tempOfilename) or die "Error deleting the temporary file $tempOfilename\n";    
}
close($OFH);

my $tag=$exe_conditions[0];
$tag=~s/[^0-9]+//g;
write_vcfrefname_dict("vcfdict.$tag.csv");
write_listrefname_dict("listdict.$tag.csv");
print("Done!\n");

## MAIN BODY AS FUNCTION FOR PARALLELISM
########################################
sub filter
{
    my $exe_condition;
    my $filtering_condition;
    
    if (defined $_)
    {
        ($exe_condition,$filtering_condition)=@{$_};
    }
    else
    {
        ($exe_condition,$filtering_condition)=@{$_[0]};
    }

    my $condition="$exe_condition$sep_param$filtering_condition";

    my $Aexecondname="$nameA$sep_param$exe_condition.vcf"; ##So far I am not using numbers for execondnames. If I change this I need to obtain the names here from the array (I may make a hash then)
    my $Bexecondname="$nameB$sep_param$exe_condition.vcf";
    my $AcovBname="$nameA$sep_param$exe_condition${sep_param}covBfiltering.tsv"; #Same here, although these are tsvs
    my $BcovBname="$nameB$sep_param$exe_condition${sep_param}covBfiltering.tsv"; 
    my $covNname="covN$sep_param$exe_condition.tsv";

    ##Reading pre-calculated shared data    
    my($refA,$refB)=($constExeCondContentRefs{$Aexecondname}, $constExeCondContentRefs{$Bexecondname});

    my $thisAname=vcf_refname("A$sep_param$condition.vcf",$condition);
    my $thisBname=vcf_refname("B$sep_param$condition.vcf",$condition);

    ##Output for this subset of paremeter conditions
    my $thisOname=$thisAname;
    $thisOname=~s/.vcf/_temp.out/;
    my %results;
    
    if (!-f $thisAname)
    {
        die("ERROR: the file $thisAname has not been generated by MITHE_control.pl\n");
    }
    if (!-f $thisBname)
    {
        die("ERROR: the file $thisBname has not been generated by MITHE_control.pl\n");
    }
    
    my %Afilt=%{parse_vcf($thisAname)};
    my %Bfilt=%{parse_vcf($thisBname)};
    my $refdataAfiltcovB=$constExeCondContentRefs{$AcovBname};
    my $refdataBfiltcovB=$constExeCondContentRefs{$BcovBname};
    my $refdatacovN=$constExeCondContentRefs{$covNname};

    #print("DEBUG: reference B $refB, reference A $refA, $filtering_condition\n");
    
    #Compare Afilt with B without filter --> Common variants + %
    my @statsAfilt;
    my ($ref_common_variantsAfilt,$ref_different_variantsAfilt)=vcf_compare_parsed($refB,\%Afilt,\@statsAfilt); ##I have to generate two hashes. One with common variants, the other with non common. Thus, the consecutive filter I can do it towards these new (smallest) hashes.
    
    #Compare Bfiltered with A without filter --> Common variants + %
    my @statsBfilt;
    my ($ref_common_variantsBfilt,$ref_different_variantsBfilt)=vcf_compare_parsed($refA,\%Bfilt,\@statsBfilt);
    #print("DEBUG: statsA of $filtering_condition @statsAfilt\n");

    #Stats filter
    my @statsfiltU;
    my @statsfiltI;
    my ($ref_common_variantsfiltU,$ref_different_variantsfiltU)=vcf_unite_parsed($ref_common_variantsAfilt,$ref_different_variantsAfilt,$ref_common_variantsBfilt,$ref_different_variantsBfilt,\@statsfiltU);
    my ($ref_common_variantsfiltI,$ref_different_variantsfiltI)=vcf_intersect_parsed($ref_common_variantsAfilt,$ref_different_variantsAfilt,$ref_common_variantsBfilt,$ref_different_variantsBfilt,\@statsfiltI);
    my @statsfiltmean=(($statsAfilt[0]+$statsBfilt[0])/2.0,($statsAfilt[1]+$statsBfilt[1])/2.0,($statsAfilt[2]+$statsBfilt[2])/2.0);

    #Substract N from Afilt. Compare the result to B without filter --> Common variants + %
    my @statsAfiltN;
    my ($ref_common_variantsAfiltN,$ref_different_variantsAfiltN)=vcf_prune($ref_common_variantsAfilt,$ref_different_variantsAfilt,\%N,\@statsAfiltN);

    #Substract N from Bfilt. Compare the result to A without filter --> Common variants + %
    my @statsBfiltN;
    my ($ref_common_variantsBfiltN,$ref_different_variantsBfiltN)=vcf_prune($ref_common_variantsBfilt,$ref_different_variantsBfilt,\%N,\@statsBfiltN);

    #Mean stats filterN
    my @statsfiltNU;
    my @statsfiltNI;
    my ($ref_common_variantsfiltNU,$ref_different_variantsfiltNU)=vcf_prune($ref_common_variantsfiltU,$ref_different_variantsfiltU,\%N,\@statsfiltNU);
    my ($ref_common_variantsfiltNI,$ref_different_variantsfiltNI)=vcf_prune($ref_common_variantsfiltI,$ref_different_variantsfiltI,\%N,\@statsfiltNI);

    my @statsfiltNmean=(($statsAfiltN[0]+$statsBfiltN[0])/2.0,($statsAfiltN[1]+$statsBfiltN[1])/2.0,($statsAfiltN[2]+$statsBfiltN[2])/2.0);

    #print(join(",","DEBUG nofilt:",@{$nofiltResultsRef{$exe_condition}},"\n"));
    #print(join(",","DEBUG filt and filtN:",@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,"\n")); #DEBUG

    #First group of PAF statistics, for filtN
    my @statsfiltNPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltN, $ref_different_variantsBfiltN, $ref_common_variantsfiltNU);
    
    #Variables that will be used just in the first iteraction of ncovB and reused for all NABs
    my @statsAfiltNAB;
    my @statsBfiltNAB;
    my @statsfiltNABU;
    my @statsfiltNABI;
    my @statsfiltNABPAF;
    my ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB);
    my ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB);
    my ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU);
    my ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI);
    my @statsfiltNABmean;
    my @statsfiltNABrefs; #NAB is inside the covB loop while they are independent, in order to have just one set of nested loops instead of two single loops and one set of nested loops (since we have, covB, NAB and their interaction, covBNAB). I save here the ref to an array with stats for each NAB condition
    my @statsfiltNABPAFrefs; #Like above
 
    for (my $ncovB=0; $ncovB < scalar @covBfiltering_conditions; ++$ncovB)
    {   
        #CovB loop variables
        my $covBfiltering_condition;
    	my $NcovBcondition;        	
        my @statsAfiltNcovB;
        my @statsBfiltNcovB;
        my @statsfiltNcovBU;
        my @statsfiltNcovBI;
        my @statsfiltNcovBPAF;
        my ($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB);
        my ($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB);
        my ($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU);
        my ($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI);
        my @statsfiltNcovBmean;

        $covBfiltering_condition=$covBfiltering_conditions[$ncovB];
    	$NcovBcondition="$exe_condition$sep_param$filtering_condition${sep_param}covB$sep_param$covBfiltering_condition";        	
     
        #AfiltNcovB
        ($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB)=vcf_prune_covB($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,$refdataBfiltcovB,\@statsAfiltNcovB,$covBfiltering_condition);
    
        #BfiltNcovB
        ($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB)=vcf_prune_covB($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,$refdataAfiltcovB,\@statsBfiltNcovB,$covBfiltering_condition);
    
        #Union
        ($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU)=vcf_unite_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBU);
        
    	#Intersection
        ($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI)=vcf_intersect_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBI);
    
        #Means
        @statsfiltNcovBmean=(($statsAfiltNcovB[0]+$statsBfiltNcovB[0])/2.0,($statsAfiltNcovB[1]+$statsBfiltNcovB[1])/2.0,($statsAfiltNcovB[2]+$statsBfiltNcovB[2])/2.0);
        
        #PAF
        @statsfiltNcovBPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltNcovB, $ref_different_variantsBfiltNcovB, $ref_common_variantsfiltNcovBU);

        #print(join(",","DEBUG filtNcovB:",@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsfiltNcovBmean,"\n")); #DEBUG
    
        #Output of list of variants and/or intermediate vcf files
        if($output_comprehensive>0)
        {
            if ($output_list)
            {
                write_variant_list($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.list",$condition);
                write_variant_list($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.list",$condition);
                write_variant_list($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.list",$condition);
            }
            if ($output_vcf)
            {
                write_variant_2vcf($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                write_variant_vcf($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                write_variant_vcf($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
            }
            if($output_comprehensive>1)
            {
                if ($output_list)
                {
                    write_variant_list($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.list",$condition);
                }
                if ($output_vcf)
                {
                    write_variant_vcf($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.vcf","$thisAname","$thisBname","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.vcf","$thisAname","$thisBname","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.vcf","$thisAname","$thisBname","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.vcf","$thisAname","$thisBname","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                }
            }
        }
    
        for (my $nNAB=0; $nNAB< scalar @NABfiltering_conditions; ++$nNAB)
        { 
            #covB*NAB loop variables
            my $NAB_condition;
            my $NAB;
            my @statsAfiltcovBNAB;
            my ($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB);
            my @statsBfiltcovBNAB;
            my ($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB);
            my @statsfiltcovBNABU;
            my @statsfiltcovBNABI;
            my @statsfiltcovBNABPAF;
            my ($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU);
            my ($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI);
            my @statsfiltcovBNABmean;

            $NAB_condition=$NABfiltering_conditions[$nNAB];
            $NAB=$constExeCondContentRefs{"${exe_condition}_${NAB_condition}"};

            if ($ncovB==0)
            {    
        #Substract NAB from AfiltN. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
                ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB)=vcf_prune_tsv_vars($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,$NAB,\@statsAfiltNAB);
        
        #Substract NAB from BfiltN. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
                ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB)=vcf_prune_tsv_vars($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,$NAB,\@statsBfiltNAB);
        
        #Mean stats filter NAB		
                ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU)=vcf_prune_tsv_vars($ref_common_variantsfiltNU,$ref_different_variantsfiltNU,$NAB,\@statsfiltNABU);
                ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI)=vcf_prune_tsv_vars($ref_common_variantsfiltNI,$ref_different_variantsfiltNI,$NAB,\@statsfiltNABI);
                @statsfiltNABmean=(($statsAfiltNAB[0]+$statsBfiltNAB[0])/2.0,($statsAfiltNAB[1]+$statsBfiltNAB[1])/2.0,($statsAfiltNAB[2]+$statsBfiltNAB[2])/2.0);
                $statsfiltNABrefs[$nNAB]=[@statsAfiltNAB,@statsBfiltNAB,@statsfiltNABU,@statsfiltNABI,@statsfiltNABmean];
        
                @statsfiltNABPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltNAB, $ref_different_variantsBfiltNAB, $ref_common_variantsfiltNABU);
                $statsfiltNABPAFrefs[$nNAB]=\@statsfiltNABPAF;

                if($output_comprehensive>1)
                {
                    if($output_list)
                    {
                        write_variant_list($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_vcf($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition); 
        			    write_variant_2vcf($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    }
                }
            }

            #print(join(",","DEBUG filtNAB:",@{$statsfiltNABrefs[$nNAB]},"\n")); #DEBUG
            
            #Substract NAB from AfiltNcovB. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
            ($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB)=vcf_prune_tsv_vars($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$NAB,\@statsAfiltcovBNAB);
    
            #Substract NAB from BfiltNcovB. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
            ($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB)=vcf_prune_tsv_vars($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,$NAB,\@statsBfiltcovBNAB);
    
            #Mean stats filter covBNAB    
            ($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU)=vcf_prune_tsv_vars($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU,$NAB,\@statsfiltcovBNABU);
            ($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI)=vcf_prune_tsv_vars($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI,$NAB,\@statsfiltcovBNABI);
    
            @statsfiltcovBNABmean=(($statsAfiltcovBNAB[0]+$statsBfiltcovBNAB[0])/2.0,($statsAfiltcovBNAB[1]+$statsBfiltcovBNAB[1])/2.0,($statsAfiltcovBNAB[2]+$statsBfiltcovBNAB[2])/2.0);
    
            @statsfiltcovBNABPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltcovBNAB, $ref_different_variantsBfiltcovBNAB, $ref_common_variantsfiltcovBNABU);

            #print(join(",","DEBUG filtcovBNAB:",@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI,@statsfiltcovBNABmean,"\n")); #DEBUG

            #Output of list of variants and/or intermediate vcf files
            
            if($output_comprehensive>0)
            {
                if($output_list)
                {
                    write_variant_list($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                }
                if($output_vcf)
                {
                    write_variant_2vcf($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                }
            }
            if($output_comprehensive>1)
            {
                if($output_list)
                {
    				write_variant_list($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                }
                if($output_vcf)
                {
                    write_variant_vcf($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                }
            }
            
            for (my $nPAF=0; $nPAF< scalar @popAFfiltering_conditions; ++$nPAF)
            {
                my $PAF_condition;
                
                my @statsAfiltcovBNABPAF;                
                my ($ref_common_variantsAfiltcovBNABPAF,$ref_different_variantsAfiltcovBNABPAF);
                
                my @statsBfiltcovBNABPAF;
                my ($ref_common_variantsBfiltcovBNABPAF,$ref_different_variantsBfiltcovBNABPAF);

                my @statsfiltcovBNABUPAF;
                my @statsfiltcovBNABIPAF;
                my ($ref_common_variantsfiltcovBNABUPAF,$ref_different_variantsfiltcovBNABUPAF);
                my ($ref_common_variantsfiltcovBNABIPAF,$ref_different_variantsfiltcovBNABIPAF);
                my @statsfiltcovBNABmeanPAF;
                my @statsfiltcovBNABPAFPAF;
    
                $PAF_condition=$popAFfiltering_conditions[$nPAF];

                #Substract PAF from AfiltcovBNAB. Compare the results to B without filter --> Common variants + %
                ($ref_common_variantsAfiltcovBNABPAF,$ref_different_variantsAfiltcovBNABPAF)=filter_with_PAF($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsAfiltcovBNABPAF);

                #Substract PAF from BfiltcovBNAB. Compare the results to A without filter --> Common variants + %
                ($ref_common_variantsBfiltcovBNABPAF,$ref_different_variantsBfiltcovBNABPAF)=filter_with_PAF($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsBfiltcovBNABPAF);

        
                #Mean stats filter covBNAB
                ($ref_common_variantsfiltcovBNABUPAF,$ref_different_variantsfiltcovBNABUPAF)=filter_with_PAF($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsfiltcovBNABUPAF);
                
                ($ref_common_variantsfiltcovBNABIPAF,$ref_different_variantsfiltcovBNABIPAF)=filter_with_PAF($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsfiltcovBNABIPAF);
        
                @statsfiltcovBNABmeanPAF=(($statsAfiltcovBNABPAF[0]+$statsBfiltcovBNABPAF[0])/2.0,($statsAfiltcovBNABPAF[1]+$statsBfiltcovBNABPAF[1])/2.0,($statsAfiltcovBNABPAF[2]+$statsBfiltcovBNABPAF[2])/2.0);

                @statsfiltcovBNABPAFPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},$ref_different_variantsAfiltcovBNABPAF, $ref_different_variantsBfiltcovBNABPAF, $ref_common_variantsfiltcovBNABUPAF);

#                #Output of list of variants and/or intermediate vcf files
                
                if($output_comprehensive>0)
                {
                    if($output_list)
                    {
                        write_variant_list($ref_common_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_2vcf($ref_common_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","$thisBname","## Somatic covBNABPAF U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_different_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_different_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                    }
                }
                if($output_comprehensive>1)
                {
                    if($output_list)
                    {
        				write_variant_list($ref_common_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_common_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_vcf($ref_common_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_common_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_different_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","$thisBname","## Somatic covBNABPAF U(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_different_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","$thisBname","## Somatic covBNABPAF I(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_common_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","$thisBname","## Somatic covBNABPAF I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                    }
                }
                #Final variables
                my @statistics;
                #Store and/or print
                @statistics=(@{$nofiltResultsRef{$exe_condition}},@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsfiltmean,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,@statsfiltNPAF,@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsfiltNcovBmean,@statsfiltNcovBPAF,@{$statsfiltNABrefs[$nNAB]},@{$statsfiltNABPAFrefs[$nNAB]},@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI,@statsfiltcovBNABmean,@statsfiltcovBNABPAF,@statsAfiltcovBNABPAF, @statsBfiltcovBNABPAF, @statsfiltcovBNABUPAF, @statsfiltcovBNABIPAF, @statsfiltcovBNABmeanPAF,@statsfiltcovBNABPAFPAF);
        #Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,(nPAFfreqs) x (AfiltN_Priv_PAF_XX,BfiltN_Priv_PAF_XX,filtN_U_PAF_XX),AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,(nPAFfreqs) x (AfiltNcovB_Priv_PAF_XX,BfiltNcovB_Priv_PAF_XX,filtNcovB_U_PAF_XX),AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,(nPAFfreqs) x (AfiltNAB_Priv_PAF_XX,BfiltNAB_Priv_PAF_XX,filtNAB_U_PAF_XX),AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean,(nPAFfreqs) x (AfiltNABcovB_Priv_PAF_XX,BfiltNABcovB_Priv_PAF_XX,filtNABcovB_U_PAF_XX),AfiltNABcovBPAF_prop,AfiltNABcovBPAF_N,AfiltNABcovBPAF_#,BfiltNABcovBPAF_prop,BfiltNABcovBPAF_N,BfiltNABcovBPAF_#,filtNABcovBPAF_propU,filtNABcovBPAF_NU,filtNABcovBPAF_#U,filtNABcovBPAF_propI,filtNABcovBPAF_NI,filtNABcovBPAF_#I,filtNABcovBPAF_prop_mean,filtNABcovBPAF_N_mean,filtNABcovBPAF_#_mean
        
        		$results{"$NcovBcondition${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param$PAF_condition"}=\@statistics;
                #print("DEBUG:$condition$OFS",array_to_string(@statistics),"\n");
            }#foreach PAF condition
        }#foreach NAB condition
    }#foreach covB condition
    
    my $thisOFILE;
    myopen(\$thisOFILE,">$thisOname") or die "ERROR opening the output file $thisOname";
    foreach my $condition (keys %results)
    {
        print($thisOFILE "$sample$OFS$condition$OFS",array_to_string(@{$results{$condition}}),"\n");
    }
    close($thisOFILE);
    push(@outfiles,$thisOname);
}#filter

###################################################################################
###FUNCTIONS
###################################################################################

#Recursive function to generate parameter combinations
#One recursion per parameter
#########################################################
sub combs 
{
    my ($id,$string,$ref_array1,$ref_array2,$ref_array_output)=@_;
    my $priv_sep_value=$sep_value;
#print("DEBUG: $id,$string,$ref_array1,$ref_array2,$ref_array_output\n");
    my @values=@{@{$ref_array2}[$id]};
    my $parameter=@{$ref_array1}[$id];
    if ($parameter=~/=$/ )
    {
        $priv_sep_value="";
    }   

    if($id<scalar(@{$ref_array1})-1) #If there are more params recurse
    {
        foreach my $value (@values)
        {
            combs($id+1,"$string$parameter$priv_sep_value${value}$sep_param",$ref_array1,$ref_array2,$ref_array_output);
        }
    }
    else #Store the string in the output array
    {
        foreach my $value (@values)
        {
            my $new_string="$string$parameter$priv_sep_value$value";
            push(@{$ref_array_output},$new_string);
        }
    }
}

#Parses parameters file, into an array of parameters
#and a matrix of values per parameter
#####################################################
sub parse_parameters_values
{
    my ($file,$ref_params,$ref_values)=@_;
    open (my $FILE, $file) or die "Error opening the file $file";
    my @content=<$FILE>;
    my $i=0;
    foreach my $line (@content)
    {
        chomp($line);
        my @temp=split($FS,$line);
        push(@{$ref_params},splice(@temp,0,1));
        ${$ref_values}[$i]=\@temp;
        ++$i;
    }
    close($FILE);
}

#Generates an string concatenating array members with OFS
#####################################################
sub array_to_string
{
    my $outstring="";
    for (my $i=0;$i<scalar(@_);++$i)
    {
        $outstring.="$_[$i]$OFS";
    }
    chop($outstring);
    return $outstring;
}

# Parse vcf, splitting multi-snv lines in different lines
#########################################################
sub parse_vcf
{
    my ($vcf1_file)=@_;
    my $VCF1;
    myopen(\$VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    
    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            ##KEY will be CHROM$OFSPOS$OFSREF$OFSALT
            $key=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            chomp($key);
            $hash{$key}=1;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return \%hash;
}

#Parse tsv
# WARNING: genotype information in TSV-related structures is different than in vcf-related ones. Here I used an array of array refs, with a different set of information. The format for each arry is REF ALT1,ALTN READS READSALT1,READSALTN.
################################################################################################################################
sub parse_tsv
{
    my ($tsv1_file)=@_;
    my $TSV1;
    myopen(\$TSV1,$tsv1_file) or die "The file $tsv1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @tsv1=<$TSV1>;
    close($TSV1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;

    for ($i=0;$i<scalar @tsv1;$i++)
    {
        chomp($tsv1[$i]);
        my @values=split("\t",$tsv1[$i]);
        $key="$values[0]$OFS$values[1]";
        splice(@values,0,2);
        if(!defined $hash{$key})
        {
            $hash{$key}=();
        }
        push(@{$hash{$key}},\@values);
        #print("DEBUG: New tsv variant data being hashed $key\n");
    }
    return \%hash;

}

#Parse vcf and get name of the sample (for only one sample)
###########################################################
sub parse_vcf_name
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $name;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $name=(split("\t",$vcf1[$i-1]))[9];
                chomp($name);
#print("DEBUG: Name $name\n");
                $flag=1;
            }
            chomp($vcf1[$i]);
            $key=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            $hash{$key}=1;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return (\%hash,$name);

}

#Compare two variant hashes, return a hash with commonvariants and another one with variants in the problem
# that aren't present in the reference and generate statistics
##############################################
sub vcf_compare_parsed
{
    my ($vcf_reference,$vcf_problem,$ref_statistics)=@_;
    my %variants2=%{$vcf_problem};
    my %commonvariants;
    my $n1=scalar(keys %{$vcf_reference});
    my $n2=scalar(keys %variants2);


    foreach my $variant (keys %{$vcf_reference})
    {
        if (exists $variants2{$variant})
        {
            $commonvariants{$variant}=1;
            delete $variants2{$variant}; ###This is a local copy of the original hash. I cut it down to reduce the number of comparisons.
#print("DEBUG: New common variant $variant");
        }

        if(scalar(keys %variants2)== 0)
        {
            #print("DEBUG: Variants2 is empty\n");
            last; ##No more pending comparisons
        }
    }	

    if ($n1<$n2)
    {
        print("WARNING: The number of variants in the unfiltered sample is smaller than the one of the filtered\n");		
    }

    my $n_common=scalar(keys %commonvariants);
    @{ $ref_statistics }=($n2>0?$n_common/$n2:"NaN",$n_common,$n2); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%commonvariants,\%variants2); ##Common, not_in_ref

}

#Filter out variants in one hash from two others and generate statistics
########################################################################
sub vcf_prune
{
    my ($vcf_1,$vcf_2,$vcf_todelete,$ref_statistics)=@_;
    my %variants1=%{$vcf_1};
    my %variants2=%{$vcf_2};
    my $tag1=0;
    my $tag2=0;

    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
        #print("DEBUG vcf_prune: Variant $variant_to_remove\n ");
		if ($tag1==0 and exists($variants1{$variant_to_remove}))
		{
			delete($variants1{$variant_to_remove});
            #print("\tDEBUG vcf_prune: deleted in vcf1");
		}
		if ($tag2==0 and exists($variants2{$variant_to_remove}))
		{
			delete($variants2{$variant_to_remove});
            #print("deleted in vcf2");
		}
		#print("\n");
		if($tag1==0 and scalar(keys %variants1)==0)
		{
			$tag1=1;
            #print("DEBUG \tvcf_prune: No more variants in vcf1\n");
		}	
		if($tag2==0 and scalar(keys %variants2)==0)
		{
			$tag2=1;
            #print("DEBUG \tvcf_prune: No more variants in vcf2\n");
		}
		
		if($tag1==1 and $tag2==1)
		{
            #print("DEBUG \tvcf_prune: No more variants\n");
			last;
		}
		
	}
	my $n_selvariants=scalar keys %variants1;
	my $n_filtvariants=scalar keys %variants2;

	@{ $ref_statistics }=($n_selvariants+$n_filtvariants>0?$n_selvariants/($n_selvariants+$n_filtvariants):"NaN",$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%variants1,\%variants2);
}

#Filter out variants in one hash depending on certain filters from two others and generate statistics
#####################################################################################################
sub vcf_prune_covB
{
    my ($vcf_1,$vcf_2,$tsv_todelete,$ref_statistics,$covBfiltering_options)=@_;
    my $refvariants_todelete=get_todelete_variants_covfiltering($covBfiltering_options,$tsv_todelete,0);
    my $ref_filtered_private_variants=remove_variants($vcf_2,$refvariants_todelete);

    
    my $n_selvariants=scalar keys %{$vcf_1};
    my $n_filtvariants=scalar keys %$ref_filtered_private_variants;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_statistics }=(0,$n_selvariants,$n_selvariants+$n_filtvariants);
    }
    else
    {
        @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    }
    return ($vcf_1,$ref_filtered_private_variants);
}

#Filter out variants in two hashes based on a hash of variables to delete.
#####################################################################################################
sub vcf_prune_tsv_vars
{
    my ($ref_common,$ref_private,$tsv_todelete,$ref_statistics)=@_;

    my $ref_filtered_common=remove_variants($ref_common,$tsv_todelete);
    my $ref_filtered_private=remove_variants($ref_private,$tsv_todelete);
    
    my $n_selvariants=scalar keys %$ref_filtered_common;
    my $n_filtvariants=scalar keys %$ref_filtered_private;
    
    if($n_selvariants+$n_filtvariants==0)
    {
        @{ $ref_statistics }=(0,$n_selvariants,$n_selvariants+$n_filtvariants);
    }
    else
    {
        @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    }
    return ($ref_filtered_common,$ref_filtered_private);
}

##Generates and returns a dictionary of variants to eliminate, given a list of parsed tsv variants and a string with filtering options. The third argument, and, indicates if the conditions of proportion and number of alternatives should be logical and or or. CovN uses && while covB || (for backwards compatibility). We should probably use the same in both.
sub get_todelete_variants_covfiltering
{
    my ($cov_filtering_options,$ref_variants,$and)=@_;
    my $min_coverage=0;
    my $max_alternative=9**9**9; ##pseudo +inf
    my $max_propalt=1;
    my %outvariants;

    #print("DEBUG: cov filtering options $cov_filtering_options\n");    

    my @covfiltering_params=split($sep_param,$cov_filtering_options);
    #Parsing filtering options
    foreach my $covoption (@covfiltering_params)
    {
        my ($param, $value)= split($sep_value,$covoption);
        $value=~s/0_//; ##covN values start with 0_ for backwards compatiblity

        #print("DEBUG: cov option $covoption, param: $param, value: $value\n");
        
        if ($param =~ /--min_coverage/i || $param =~ /--max_coverage/i) ##The second term is included for covN backwards compatibility and it is misleading
        {
           $min_coverage=$value;
        }
        elsif ($param =~ /--max_alternative/i || $param =~ /--min_reads_alternate/i)
        {
           $max_alternative=$value;
        }
        elsif ($param =~ /--max_propalt/i || $param =~ /--min_freq_alt/i)
        {
           $max_propalt=$value;
        }
        elsif ($param =~ /--default/i)
        {
            next;
        }
        else
        {
            die "The param $param has not been recognized properly in a cov(N|B) filtering step\n";
        }
    }

    #Selecting the comparison function
    my $comparefunction;
    if ($and==1)
    {
        $comparefunction=\&remove_alts_covN;
    }
    else
    {
        $comparefunction=\&remove_alts_covB;
    }

    my $ref;
    my $alt_string;
    my $ref_reads;
    my $alt_reads_string;
    my $totalreads;
    my @alts;
    my @reads;
    my $shortvariant;
    
    foreach my $variant (keys %{$ref_variants})
    {
        foreach my $genotype (@{$ref_variants->{$variant}})
        {
            #value structure: (REF,"ALT0,ALTN",#REFREADS,"#ALT1READS,#ALTNREADS")
            ($ref,$alt_string,$ref_reads,$alt_reads_string)=@{$genotype};
            @alts=();
            @reads=();
            $shortvariant=$variant;
            $shortvariant=~s/$OFS[^$OFS]+$OFS[^$OFS]+$//;
    
            #Prepare lists of alternatives I may not need it for the loop, but I still need to calculate the number of total reads
            if($alt_string=~/,/) ##Multiple alternatives
            {
                @alts=split(",",$alt_string);
                @reads=split(",",$alt_reads_string);
                $totalreads=$ref_reads;
                $totalreads+=$_ foreach @reads; #Oneline for loop
            }
            else
            {
                $totalreads=$alt_reads_string+$ref_reads;
                $alts[0]=$alt_string;
                $reads[0]=$alt_reads_string;
            }
    
            #If there are coverage problems, we eliminate the position and we are done with this position
            if($totalreads<$min_coverage)
            {
                $outvariants{$shortvariant}=1;
            }
            else ##Presence of variant alleles
            {
                foreach (my $ialt=0; $ialt<scalar @alts; ++$ialt)
                {
                    if($comparefunction->($max_alternative,$max_propalt,$reads[$ialt],$totalreads==0?0:$reads[$ialt]/($totalreads+0.0))) #genotype eliminated due to alternative variants
                    {
                       $outvariants{join($OFS,$shortvariant,$ref,$alts[$ialt])}=1; 
                    }
                }
            }
        }#genotype
    }#variant

    return \%outvariants;
}

##Short snippets

#Comparison to determine if a variant should be deleted or not, covN style
sub remove_alts_covN
{
    #my ($max_alternative, $max_propalt, $altreads, $altprop)=@_;
    #return $altreads>$max_alternative && $altprop>$max_propalt;
    return $_[2]>$_[0] && $_[3]>$_[1];
    
}

#Comparison to determine if a variant should be deleted or not, covB style
sub remove_alts_covB
{
    #my ($max_alternative, $max_propalt, $altreads, $altprop)=@_;
    #return $altreads>$max_alternative || $altprop>$max_propalt;
    return $_[2]>$_[0] || $_[3]>$_[1];
}

#Retries the opening of a file a number of times, instead of failing at the first attempt
#This seems to be needed when running on a network attached system in parallel with many processes working on it
sub myopen
{
    my ($refFH,$name)=@_;
    my $i=0;
    while (!open($$refFH,$name))
    {
        warn "Retrying to open the file $name, attempt $i\n";
        sleep 10;
        $i+=1;
        if ($i > $MAX_OPEN_IT)
        {
            return 0;
        }
    }
    return 1;
}

#Retries the execution of a system command MAX_SYSTEM_IT times
sub mysystem
{
    my ($command)=@_;
    my $signal=1;
    my $i=0;

    while ($signal != 0)
    {
        if ($i > $MAX_SYSTEM_IT)
        {
            return $signal;
        }
        elsif ($i>1)
        {
            warn "Retrying to execute the system command $command, attempt $i\n";
        }
        $signal=system($command);
        $i+=1;
        sleep 10;
    }
    return 0;
}

#Filters out vcf variants using variants originated in a tsv file after filtering using get_todelete_variants_covBfiltering
#It returns the resulting variants in a new array ref
#First checks for the shortvar version of the variant, then for the full version of it
sub remove_variants
{
    my ($ref_problem_vars,$ref_todelete_vars)=@_;
    my %outvariants=%{$ref_problem_vars}; #Shallow copy
    my $shortvar;

#    print("\nDEBUG: execution of remove_variants\n\tDEBUG: variants to consider\n");
    foreach my $variant (keys %outvariants)
    {
        $shortvar=$variant;
        $shortvar=~s/$OFS[^$OFS]+$OFS[^$OFS]+$//;
        if(exists $ref_todelete_vars->{$shortvar} or exists $ref_todelete_vars->{$variant})
        {
            delete($outvariants{$variant});
        }
#        print("\tDEBUG: $variant, $shortvar\n");
    }
#        print("DEBUG: variants to delete\n");
#    foreach my $variant (keys %{$ref_todelete_vars})
#    {
#        print("\tDEBUG: $variant\n");
#    }
#    print("DEBUG: Done\n\n");
    return \%outvariants;
}

#Filter out variants in one hash from another
#############################################
sub vcf_prune_single
{
    my ($vcf_1,$vcf_todelete)=@_;
    my %variants1=%{$vcf_1};

    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
#print("DEBUG: Variant $variant_to_remove ");
		if (exists($variants1{$variant_to_remove}))
		{
			delete($variants1{$variant_to_remove});
            #print("deleted in vcf1");
		}
		#print("\n");
		if(scalar(keys %variants1)==0)
		{
            last;
            #print("DEBUG: No more variants in vcf1\n");
		}
		
	}
    return (\%variants1);
}

#Combines the variants of two samples. Calculates the union of the two samples
#(union of the commons and union of the differences) and remove spurious
#differences (that can potentially appear in the union).Then calculate the stats.
########################################################################
sub vcf_unite_parsed
{
    my ($ref_common_variantsA,$ref_different_variantsA,$ref_common_variantsB,$ref_different_variantsB,$ref_statistics)=@_;
    my %common_variants=(%{$ref_common_variantsA},%{$ref_common_variantsB}); ##Union
    my %different_variants=(%{$ref_different_variantsA},%{$ref_different_variantsB});#Union

    foreach my $common_variant (keys %common_variants) ##We have to delete possible variants that are not different any more
    {
        delete($different_variants{$common_variant}); 
    }
    
    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %different_variants;
    @{ $ref_statistics }=($n_selvariants+$n_filtvariants>0?$n_selvariants/($n_selvariants+$n_filtvariants):"NaN",$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants

    return (\%common_variants,\%different_variants);
}

#Combines the variants of two samples. Calculates the intersection of the two samples
#(intersection of the commons and recalculates the differences from all the variants) and remove spurious
#Then calculate the stats.
########################################################################
sub vcf_intersect_parsed
{
    my ($ref_common_variantsA,$ref_different_variantsA,$ref_common_variantsB,$ref_different_variantsB,$ref_statistics)=@_;
    my %common_variants;
    
    foreach my $key (keys %$ref_common_variantsA)
    {
        if (exists $ref_common_variantsB->{$key})
        {
            $common_variants{$key}=$ref_common_variantsA->{$key};
        }
    }
    
    my %different_variants=(%{$ref_different_variantsA},%{$ref_different_variantsB},%{$ref_common_variantsA},%{$ref_common_variantsB});#Union of all variants

    foreach my $common_variant (keys %common_variants) ##We recalculate the different variants (total-common)
    {
        delete($different_variants{$common_variant});
    }
   
    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %different_variants;
    @{ $ref_statistics }=($n_selvariants+$n_filtvariants>0?$n_selvariants/($n_selvariants+$n_filtvariants):"NaN",$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%common_variants,\%different_variants);
}

#Removes the VCF header
##############################################
sub clean_vcfcontent
{
	my ($array)=@_;
	my @array;
	my $last_comment=0;
	for (my $i=0;$i<scalar @{$array};$i++)
	{
		if(${$array}[$i]=~/^##/)
		{
			++$last_comment;
		}
		else
		{
			last;
		}
	}
	splice(@{ $array },0,$last_comment);
}

# Generates a hash of variants from an array with the content of a VCF without the header
#########################################################################################
sub variants_to_hash
{
	my ($array)=@_;
	my %hash;
	my $key;
	for (my $i=1; $i<scalar @{$array};$i++) ###The first line is the header
	{
        chomp(${$array}[$i]);
		$key=${$array}[$i];
        $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
		$hash{$key}=1;
		#print("DEBUG: New variant being hashed $key\n");
	}
	return \%hash;
}


# Writes a list of variants in csv format contained in a hash
# ############################################################
sub write_variant_list
{
    my ($ref_hash,$filename,$filter)=@_;
    my $FILE=openwrite_list_refname($filename,$filter);
    (! defined $FILE) and return;
    print($FILE "#CHROM,POS\n");
    foreach my $variant (keys %{$ref_hash})
    {
        print($FILE join($OFS,split(",",$variant)),"\n");
    }
    close $FILE;
}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################

sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment,$filter)=@_;
    my $OFILE=openwrite_vcf_refname($filename,$filter);
    (! defined $OFILE) and return;
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

# Writes a vcf with the variables contained in a hash selected from two VCF files
# ATTENTION: Variants in the first VCF file have higher priority when they are 
# present in the two VCF files. The header comes also from the first VCF file. ONLY one sample information, comming from the first VCF file with higher priority.
# This subrutine is compatible with genotypes. If a given genotype is present in the two
# VCF files, the information will come from the first VCF file. THIS IS SOME KIND OF FRANKENSTEIN'S
# MONSTER VCF FILE. IT WILL DEFINETLY CONTAIN INCONSISTENCIES. USE IT ONLY AS A LIST OF VARIANTS
# ##################################################################################
sub write_variant_2vcf
{
    my ($ref_hash,$filename,$vcf,$vcf2,$comment,$filter)=@_;
    my $OFILE=openwrite_vcf_refname($filename,$filter);
    my %outcontent;
    (! defined $OFILE) and return;
    open(my $IFILE, $vcf);
    open(my $IFILE2, $vcf2);
    my @icontent= <$IFILE>;
    my @icontent2= <$IFILE2>;
    close($IFILE);
    close($IFILE2);
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
        else ##Most lines
        {
            $key=$icontent[$i];
            chomp($key);
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                $outcontent{$key}=$icontent[$i];
                delete($hash{$key});
            }
            if(scalar keys %hash == 0)
            {
                last;
            }
            
        }

    }
    if (scalar keys %hash !=0) #Not present at all in vcf1 (the needed alleles), needs to be printed from vcf2
    {
        for(my $i=0;$i< scalar @icontent2; ++$i)
    	{
		    unless ($icontent2[$i]=~/^#/)
		    {
			    $key=$icontent2[$i];
                chomp($key);
                $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            	#print("DEBUG: Key $key\n");
            	if(exists $hash{$key})
            	{
                    $outcontent{$key}=$icontent2[$i];
                    delete($hash{$key});
            	}
            	if(scalar keys %hash == 0)
            	{
                	last;
            	}
		    }
	    }
    }
    
    ##Sort outcontent and print to OFILE

    my @int_v;

    foreach $key (nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %outcontent)
    {
        print($OFILE $outcontent{$key});
    }

    
    close $OFILE;
}

## Gets the position in a vcf file of the sample with the given name
####################################################################

sub get_col_name_vcf
{
    my ($vcf,$name)=@_;
    open(my $VCF,$vcf) or die "The file $vcf is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf=<$VCF>;
    close($VCF);
    my $flag=0;
    my $i;
    my $j;
    my @cols;
    my $rescol=-1;;

    for ($i=0;$i<scalar @vcf;$i++)
    {
        unless($flag==0 and $vcf[$i]=~/^#/)
        {
            if ($flag==0)
            {
                        @cols=split("\t",$vcf[$i-1]);
                        for ($j=0;$j<scalar @vcf;$j++)
                        {
                            if($cols[$j]=~/$name/)
                            {
                                $rescol=$j-9; ### VCF files have 9 columns before the first genotype
                                last;
                            }
                        }
                        #print("DEBUG: Name $name\n");
                        $flag=1;
                        last;
            }
        }
    } 
    if ($rescol < 0)
    {
        die "Column with genotype $name not detected \n";
    }
    return $rescol;
}

# Change the genotype column from 0 to the one indicated by $pos. 
# Pos is obtained with parse_vcf_name + get_col_name_vcf
################################################################

sub filtervalues_changepost
{
    my ($ref_array,$pos)=@_;
    my @array=@{$ref_array};
    for (my $i=0; $i< scalar @array; $i++)
    {
        my @temp=@{$array[$i]};
    #    print("DEBUG: before modification ",join(",",@temp));
        for (my $j=0; $j<scalar @temp; $j++)
        {
            $temp[$j]=~s/0_/${pos}_/;
        }
        $array[$i]=\@temp;
    #    print(" after modification ",join(",",@{$array[$i]}),"\n");
    }
    return @array;
   
}

# Opens a vcf file to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names is printed by write_vcfrefname_dict
# If the file is already in the hash/array system, returns undef
# ######################################################################

sub openwrite_vcf_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamevcf{$filename} && basename($dictrealnamevcf{$filename})=~/^$filter/)
    {
        #warn "DEBUG: skipping file $filename since it had been previously written by this process\n";
        return undef;
    }
    my $n_name=push (@vcfnames,$filename);
    my $realfilename=$outputdir . "/" . $filter . $n_name . ".vcf";
    $dictrealnamevcf{$filename}=$realfilename;
    open (my $FILE, ">$realfilename") or die "Error opening $realfilename, in place of $filename\n";
    return $FILE;
}

# Returns the new vcf filename to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names has to be then printed using write_vcfrefname_dict
# If the file is already in the hash/array system, returns the name previously stored
# ###############################################################################################

sub vcf_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamevcf{$filename} && basename($dictrealnamevcf{$filename})=~/^$filter/)
    {
        #warn "DEBUG: this file, $filename had been previously written by this process\n";
        return $dictrealnamevcf{$filename}
    }
    my $n_name=push (@vcfnames,$filename);
    my $realfilename=$outputdir . "/" . $filter . $n_name . ".vcf";
    $dictrealnamevcf{$filename}=$realfilename; 
    return $realfilename;
}

# Opens a list file to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names is printed by write_listrefname_dict
# ######################################################################
sub openwrite_list_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamelist{$filename} && $dictrealnamelist{$filename}=~/^$filter/)
    {
        #warn "DEBUG: skipping file $filename since it had been previously written by this process\n";
        return undef; 
    }
    my $n_name=push(@listnames,$filename);
    my $realfilename=$filter . $n_name . ".list";
    $dictrealnamelist{$filename}=$realfilename;
    open (my $FILE, ">$realfilename") or die "Error opening $realfilename, in place of $filename\n";
    return $FILE;
}


# Generates the output file with the mapping of vcf file numbers and ideal names
# ##############################################################################
sub write_vcfrefname_dict
{
    my ($filename)=@_;
    open (my $FILE, ">$filename");
    foreach my $filename (keys %dictrealnamevcf)
    {
        print($FILE "$filename$OFS$dictrealnamevcf{$filename}\n");
    }
    foreach my $file (@namedvcffiles)
    {
        print($FILE "$file$OFS$outputdir/$file\n");
    }
    close($FILE);
}

# Generates the output file with the mapping of list file numbers and ideal names
# ##############################################################################

sub write_listrefname_dict
{
    my ($filename)=@_;
    open (my $FILE, ">$filename");
    foreach my $filename (keys %dictrealnamelist)
    {
        print($FILE "$filename$OFS$dictrealnamelist{$filename}\n");
    }
    close($FILE);
}

# Parses all const_execond
# ###################################################################

sub parse_const_execond
{
    my $exe_condition;
    if (defined $_)
    {
        $exe_condition=$_;
    }
    else
    {
        $exe_condition=$_[0];
    }
    
    my $Aexecondname="$nameA$sep_param$exe_condition.vcf";
    my $Bexecondname="$nameB$sep_param$exe_condition.vcf";
    my $AcovBname="$nameA$sep_param$exe_condition${sep_param}covBfiltering.tsv";
    my $BcovBname="$nameB$sep_param$exe_condition${sep_param}covBfiltering.tsv";
    my $covNname="covN$sep_param$exe_condition.tsv";
    
    (! -f "$wd/$covNname") and die "Error calculating coverage of N";
    (! -f "$wd/$AcovBname") or (! -f "$wd/$BcovBname") and die "Error calculating coverages of A and B";
    (! -f "$wd/$Aexecondname") or  (! -f "$wd/$Bexecondname") and die "Missing vcf files for constant files. Check the variant calling step";

    $constExeCondContentRefs{$AcovBname}=parse_tsv("$wd/$AcovBname");
    $constExeCondContentRefs{$BcovBname}=parse_tsv("$wd/$BcovBname");
    $constExeCondContentRefs{$covNname}=parse_tsv("$wd/$covNname");
    
    #TODO: I am not parallelizing this, but I am getting a lot of benefits from reusing it. I think that parallalelizing it would spend more extra time spliting the fork than the actual benefit due to the several processes
    ##NAB constants
    foreach my $cond (@NABfiltering_conditions)
    {
        $constExeCondContentRefs{"${exe_condition}_${cond}"}=get_todelete_variants_covfiltering($cond,$constExeCondContentRefs{$covNname},1);
    }

    $constExeCondContentRefs{$Aexecondname}=parse_vcf("$wd/$Aexecondname");
    $constExeCondContentRefs{$Bexecondname}=parse_vcf("$wd/$Bexecondname");
    
    ##PAF constants
#    print("DEBUG: parsing PAF data for exe_condition: $exe_condition\n");
    $constExeCondContentRefs{"${exe_condition}_PAF"}=getPAFdata($constExeCondContentRefs{$Aexecondname},$constExeCondContentRefs{$Bexecondname});
#    print("DEBUG: done parsing PAF data for exe_condition: $exe_condition\n");

    ##TODO: Do I want to parallelize this independently? Same question about the benefits and costs of forking
    foreach my $cond (@popAFfiltering_conditions)
    {
        $constExeCondContentRefs{"${exe_condition}_${cond}"}=filt_PAF($constExeCondContentRefs{"${exe_condition}_PAF"},$cond);
    } 
    
    $nofiltResultsRef{$exe_condition}=[scalar keys %{$constExeCondContentRefs{$Aexecondname}}, scalar keys %{$constExeCondContentRefs{$Bexecondname}}, scalar keys %N, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Aexecondname},\%N)}, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Bexecondname},\%N)}];
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
    my ($tstart,$talt,$tref,$tfilt,$taf);

    foreach my $ref_hash (@_)
    {
        foreach my $key (keys %{$ref_hash}) ##Foreach key
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
        }##foreach key
    }##foreach hash
    return \%outdata;
}

##This subrutine gets an array of PAF conditions and generates the frequency of private A, B and common C variants with PAF<=condition for each condition plus the mean PAF for A, B, and C respectively.
sub getPAFStats
{
    my ($ref_conditions, $ref_pAFfiltN, $ref_A, $ref_B, $ref_C)=@_;
    my @output;
    my @filt_values;
    my $ngood;
    my $ntotal;
    my $key;
    my $filtvalue;
    my $cumPAF;
    my $meanPAFA;
    my $meanPAFB;
    my $meanPAFC;

    #ref_pAFfiltN; #key: CHROM${OFS}POS${OFS}REF${OFS}ALT value: [$tfilt,$taf];

    foreach my $filtcond (@{$ref_conditions})
    {
        $filtvalue=$filtcond;
        ($filtvalue=~s/^--max_pAF://) == 1 or die "The PAF parameter file is not specified properly, $filtcond. The only implemented parameter is --max_pAF";
        push(@filt_values,$filtvalue);
    }

    for (my $i=0; $i<scalar@filt_values; ++$i)
    {
        my $value=$filt_values[$i];

        #PrivA
        $ngood=0;
        $ntotal=0;
        $cumPAF=0;

        foreach $key (keys %{$ref_A})
        {
            if(exists $ref_pAFfiltN->{$key})
            {
                if($ref_pAFfiltN->{$key}->[1] eq "NA")
                {
                    $i==0 and $cumPAF+=$GNOMAD_MIN;
                    ++$ngood;
                }
                else
                {
                    if($ref_pAFfiltN->{$key}->[1] <= $value)
                    {
                        ++$ngood;
                    }
#                    else
#                    {
#                        ##Bad variant
#                    }
#
                    $i==0 and $cumPAF+=$ref_pAFfiltN->{$key}->[1];
                }
            }
            else
            {
                $i==0 and $cumPAF+=$GNOMAD_MIN;
                ++$ngood;
                warn "No population data for $key\n";
            }
            ++$ntotal;
        }

        if($ntotal != 0)
        {
            push(@output,$ngood*1.0/$ntotal);
            $i==0 and $meanPAFA=$cumPAF*1.0/$ntotal;
        }
        else
        {
            push(@output,"NaN");
            $i==0 and $meanPAFA="NaN";
        }
        
        #PrivB
        $ngood=0;
        $ntotal=0;
        $cumPAF=0;
        
        foreach $key (keys %{$ref_B})
        {
            if(exists $ref_pAFfiltN->{$key})
            {
                if($ref_pAFfiltN->{$key}->[1] eq "NA")
                {
                    $i==0 and $cumPAF+=$GNOMAD_MIN; 
                    ++$ngood;
                }
                else
                {
                    if($ref_pAFfiltN->{$key}->[1] <= $value)
                    {
                        ++$ngood;
                    }
                    $i==0 and $cumPAF+=$ref_pAFfiltN->{$key}->[1]; 
                }
            }
            else
            {
                $i==0 and $cumPAF+=$GNOMAD_MIN;
                ++$ngood;
                warn "No population data for $key\n";
            }
            ++$ntotal;
        }
        if($ntotal != 0)
        {
            push(@output,$ngood*1.0/$ntotal);
            $i==0 and $meanPAFB=$cumPAF*1.0/$ntotal;
        }
        else
        {
            push(@output,"NaN");
            $i==0 and $meanPAFB="NaN";
        }
        
        #Common
        $ngood=0;
        $ntotal=0;
        $cumPAF=0;

        foreach $key (keys %{$ref_C})
        {
            if(exists $ref_pAFfiltN->{$key})
            {
                if($ref_pAFfiltN->{$key}->[1] eq "NA")
                {
                    $i==0 and $cumPAF+=$GNOMAD_MIN;
                    ++$ngood;
                }
                else
                {
                    if($ref_pAFfiltN->{$key}->[1] <= $value)
                    {
                        ++$ngood;
                    }
                    $i==0 and $cumPAF+=$ref_pAFfiltN->{$key}->[1];
                }
            }
            else
            {
                $i==0 and $cumPAF+=$GNOMAD_MIN;
                ++$ngood;
                warn "No population data for $key\n";
            }
            ++$ntotal;
        }
        
        if($ntotal != 0)
        {
            push(@output,$ngood*1.0/$ntotal);
            $i==0 and $meanPAFC=$cumPAF*1.0/$ntotal;
        }
        else
        {
            push(@output,"NaN");
            $i==0 and $meanPAFC="NaN";
        }

    } ##Foreach filter

    push(@output,$meanPAFA,$meanPAFB,$meanPAFC);
    return @output;
}

##Subrutine to generate the header of the PAF stats for the output file
##Inputs: two array references, to a) names of the different output columns
## b) array of PAF filtering conditions
## Return: array of column headers that combines the names with the filtering values
sub makeHeaderPAF
{
    my ($ref_names,$ref_conditions)=@_;
    my @outcolumns;
    my $outcondition;
    foreach my $condition (@{$ref_conditions})
    {
        $outcondition=$condition;
        $outcondition =~s/--[^:]+://;
        foreach my $name (@{$ref_names})
        {
            push(@outcolumns,$name.$outcondition);
        }
    }
    foreach my $name (@{$ref_names})
    {
        push(@outcolumns, $name."Mean");
    }
    return @outcolumns;
}

##Subrutine that generates a list of variants that are valid (<= condition) considering their population allele frequency
sub filt_PAF
{
    my ($ref_data,$condition)=@_;
    $condition =~s/--[^:]+://; ##Get the maximum population allele frequency
    my %outdata;
    foreach my $key (keys %{$ref_data})
    {
        if($ref_data->{$key}[1] eq "NA" || $ref_data->{$key}[1]<=$condition)
        {
            $outdata{$key}=$ref_data->{$key};
        }
    }

    return \%outdata;
}

##Filters two hashes of variants using a hash of variants that are valid (this will be a superset of the variants that will be valid)  and generates stats like vcf_prune
sub filter_with_PAF
{
    my ($ref_common, $ref_different, $ref_PAF, $ref_stats)=@_;
    my (%out_common, %out_different);
   
#    print("DEBUG: filter_with_PAF\n"); 
    ##The number of ref_PAF variants will be much larger than ref_common or ref_different. Thus, I will loop along those two, instead of the first.
    foreach my $key (keys %{$ref_common})
    {
#        print("DEBUG:\t variant $key\n");
        if(exists $ref_PAF->{$key}) ##I don't need to check the value of this since this list of variants has been pre-made
        {
            $out_common{$key}=1;
#            print("DEBUG:\t\t with PAF information ",join(",",@{$ref_PAF->{$key}}),"\n");
        }
#        else
#        {
#            print("DEBUG:\t\t This variant had been filtered out and therefore is not present in PAF information\n");
#        }
    }
    foreach my $key (keys %{$ref_different})
    {
#        print("DEBUG:\t variant $key\n");
        if(exists $ref_PAF->{$key}) ##I don't need to check the value of this since this list of variants has been pre-made
        {
            $out_different{$key}=1;
#            print("DEBUG:\t\t with PAF information ",join(",",@{$ref_PAF->{$key}}),"\n");
        }
#        else
#        {
#            print("DEBUG:\t\t This variant had been filtered out and therefore is not present in PAF information\n");
#        }
    }
    	
    my $n_selvariants=scalar keys %out_common;
	my $n_filtvariants=scalar keys %out_different;
	@{ $ref_stats }=($n_selvariants+$n_filtvariants>0?$n_selvariants/($n_selvariants+$n_filtvariants):"NaN",$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
#    print("DEBUG stats:",join(",",@{$ref_stats}),"\n");
    return (\%out_common, \%out_different);
}

##Parses a previously generated vcfdict file into a hash (returned)
sub parse_vcfdict
{
    my ($vcfdict)=@_;
    open(my $FT, $vcfdict) or die "ERROR: the file $vcfdict cannot be opened\n";
    my @content=<$FT>;
    my ($realname,$keyname);
    my %dict;
    foreach my $line (@content)
    {
        chomp($line);
        ($keyname,$realname)=split(",",$line);
        $dict{$keyname}=$realname;
    }
    close($FT);
    return %dict;
}

##Refines previously generated vcfdict for this specific comparison, eliminating samples that are not being compared here, and changing their names from S0 ... SN to A and B.
sub refine_vcfdict
{
    my ($refdict)=@_;
    my $newkey;
    foreach my $key (keys %$refdict)
    {
        if($key =~ m/^S/) #Needs to be either modified or deleted, because it is a sample
        {
            if($key =~ m/^$nameA/) #New A sample
            {
                $newkey=$key;
                $newkey =~ s/^$nameA/A/g;
                $refdict->{$newkey}=$refdict->{$key};
            }
            elsif($key =~ m/^$nameB/) #New B sample
            {
                $newkey=$key;
                $newkey =~ s/^$nameB/B/g;
                $refdict->{$newkey}=$refdict->{$key};
            }
            delete $refdict->{$key};
        }
    }
}
