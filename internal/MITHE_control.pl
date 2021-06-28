#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;

##Configuration variables
######################################################

#Workload manager
my $queue="";
my $queue_argument=$ENV{'MITHE_SUBMIT_PAR'};
my $mithe_time=$ENV{'MITHE_MAX_TIME'};
my $n_cores=1;
my $scheduler=$ENV{'MITHE_SUBMIT_CMD'};
my $qsub=$ENV{'MITHE_SUBMIT_MUL'};
my $qsub_noparallel=''; 
my $sed=$ENV{'MITHE_SUBMIT_SED'};
my $sep_dep=":";
my $dep_prefix=$ENV{'MITHE_SUBMIT_DEP'};

#Internal dependencies
my $MITHE_INT=$ENV{'MITHE_INT'};
my $variant_caller="platypus";
my $variant_calling_sh;
my $helper_sh="perl.sh";
my $helper_pl="MITHE.pl";
my $tstv_sh="tstv.sh";
my $annotation_sh="annovar_loop.sh";
my $covBED_sh="covBED.sh";
my $covB_sh="covB_half.sh";
my $covN_sh="covN.sh";
my $reduce_sh="reduce.sh";
my $vcf_filt_exe="vcf_filtering.sh";

#IO
my $sep_param="#";
my $sep_value=":";
my $OFS=",";
my $FS=",";
my $sleep=60;
######################################################

##Global variables
#######################################################
my %dictrealnamevcf;
my @vcfnames;
my @namedvcffiles;

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile="";
my $covBfiltercond_inputfile="";
my $popAFfiltercond_inputfile="";
my $output_dir="vcf_outputdir";
my $output_file="";
my $normal_bam="";
my @sample_bams;
my $output_vcf=0;
my $output_list=0;
my $comp=0;
my $filterINDELS=0;

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file --normal_bamfile bamfile_normal_sample --output_vcf --output_list --comp [--queue] file1 file2 ... fileN\n--output_vcf: (bool) generate resulting vcf files or not\n--output_list: (bool) generate resulting list of variants or not\n--comp: (int) indicating the comprehensiveness of the output, 0=no files, 1=only needed files to call variants, 2= all intermediate variants\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the control sample\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--popAF_cond_inputfile: input file for the filter of population allele frequencies using gnomAD\n\t--output_dir: output directory for vcf files\n\t--n_cores: number of cores to execute some steps in parallel\n\t--queue: optional, name of the queue jobs should be submitted\n\t--filterINDELS: if activated, INDELS are discarded from the whole process\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
    'exec_cond_inputfile|e=s' => \$execond_inputfile,
	'filt_cond_inputfile|f=s' => \$filtercond_inputfile,
	'NABfilt_cond_inputfile=s' => \$NABfiltercond_inputfile,
    'covaltB_cond_inputfile=s' => \$covBfiltercond_inputfile,
    'popAF_cond_inputfile=s' => \$popAFfiltercond_inputfile,
    'output_dir=s' => \$output_dir,
	'output_file|o=s' => \$output_file,
	'normal_bamfile=s' => \$normal_bam,
    'n_cores=i' => \$n_cores,
    'output_vcf=i'=>\$output_vcf,
    'output_list=i' => \$output_list,
    'comp=i' => \$comp,
    'queue=s' => \$queue,
    'filterINDELS=i' => \$filterINDELS,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($normal_bam eq "") || $help) and die $usage;

##Parsing input bam files
#########################

@sample_bams=@ARGV;

foreach my $samplefile (@sample_bams)
{
    if ( ! -f $samplefile)
    {
        die "The BAM file $samplefile is not accesible. Please, check input options.\n";
    }
    else
    {
        $samplefile=Cwd::abs_path($samplefile);
    }

}

##Input file parsing and directory creation
######################################################

my @exe_parameters=("input");
my @filtering_parameters=("");
my @exe_param_values=([("")]);
my @filtering_param_values=([("")]);

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

##BAMs
if ( ! -f $normal_bam)
{
	die "The BAM file $normal_bam is not accesible. Please, check input options.\n";
}
else
{
	$normal_bam=Cwd::abs_path($normal_bam);
}


mkdir $output_dir;
my $oefile=Cwd::abs_path($execond_inputfile);
my $offile=Cwd::abs_path($filtercond_inputfile);
my $onfile=Cwd::abs_path($NABfiltercond_inputfile);
my $ocfile=Cwd::abs_path($covBfiltercond_inputfile);
my $opfile=Cwd::abs_path($popAFfiltercond_inputfile);

chdir $output_dir or die "The output directory $output_dir is not accesible";
$output_dir=Cwd::abs_path($output_dir);

if (-f "$MITHE_INT/$variant_caller.sh")
{
    $variant_calling_sh="$MITHE_INT/$variant_caller.sh";
}
else
{
    die "Error, the sh file for the variant caller $variant_caller can't be located. Please, make sure that the environment variable MITHE_INT is indicating the folder with this package of scripts\n"
}

if (-f "$MITHE_INT/$helper_sh" && -f "$MITHE_INT/$helper_pl" && -f "$MITHE_INT/$reduce_sh" && -f "$MITHE_INT/$covB_sh" && -f "$MITHE_INT/$covN_sh" && -f "$MITHE_INT/$covBED_sh")
{
    $helper_sh="$MITHE_INT/$helper_sh";
    $helper_pl="$MITHE_INT/$helper_pl";
    $reduce_sh="$MITHE_INT/$reduce_sh";
    $covB_sh="$MITHE_INT/$covB_sh";
    $covN_sh="$MITHE_INT/$covN_sh";
    $covBED_sh="$MITHE_INT/$covBED_sh";
}
else
{
    die "Error, the files for the secondary analysis of the data, $helper_sh, $helper_pl, $reduce_sh, $covB_sh, $covN_sh, and $covBED_sh can't be located. Please, make sure that the environment variable MITHE_INT is indicating the folder with this package of scripts\n"
}

if (-f "$MITHE_INT/$tstv_sh")
{
    $tstv_sh="$MITHE_INT/$tstv_sh";
}
else
{
  die "Error, the sh file to get the TsTv data, $tstv_sh can't be located. Please, make sure that the environment variable MITHE_INT is indicating the folder with this package of scripts\n"
}

if (-f "$MITHE_INT/$annotation_sh")
{
    $annotation_sh="$MITHE_INT/$annotation_sh";
}
else
{
  die "Error, the sh file to annotate variants using annovar, $annotation_sh can't be located. Please, make sure that the environment variable MITHE_INT is indicating the folder with this package of scripts\n"
}


if (-f "$MITHE_INT/$vcf_filt_exe")
{
    $vcf_filt_exe="$MITHE_INT/$vcf_filt_exe";
}
else
{
    die "The executable vcf_filtering.sh can't be located. Please, make sure that the environment variable MITHE_INT indicates the directory with this package of scripts\n";
}


if($output_vcf==0 || $comp != 2)
{
    print "Warning: not generating output vcf files or restricting them (--comp<2) will eliminate the tstv posterior analyses\n"
}

##Job submitting command finalization
if ($queue ne "")
{
    $qsub=join(" ",$scheduler,$queue_argument.$queue,$qsub.$n_cores);
    $qsub_noparallel=join(" ",$scheduler,$queue_argument.$queue,$qsub_noparallel);
}
else
{
    $qsub=join(" ",$scheduler,$qsub.$n_cores);
    $qsub_noparallel=join(" ",$scheduler,$qsub_noparallel);
}

#
## Main conditions loop
#######################

my $name_condition;
my @exe_conditions;
my @filtering_conditions;
my $bamfiles;

###Generate all combinations of proposed values for execution
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);

my %job_ids;
my $exe_condition;

mkdir("e_logs");
mkdir("o_logs");

### Default calling parameters. This may be changed in the future
#########################################################

print("Unfiltered variant calling detection/execution:\n");	
if(! -f "N.vcf")
{
	$bamfiles=$normal_bam;
	$exe_condition="";
    my $job_id=submit_job("$qsub $variant_calling_sh $bamfiles N.vcf N_platypus.log $filterINDELS $exe_condition");
    print("\tNormal tissue variant calling submited with job_id $job_id\n");
}
else
{
	print("\tNormal tissue variant calling already present, skipping it\n");
}
push(@namedvcffiles,"N.vcf");

my $samplename;

for (my $i=0; $i< scalar @sample_bams; ++$i)
{

    $samplename=sprintf("S%0*d", length scalar @sample_bams, $i);
    if(! -f "${samplename}.vcf")
    {
        $bamfiles=$sample_bams[$i];
        $exe_condition="";
        my $job_id=submit_job("$qsub $variant_calling_sh $bamfiles ${samplename}.vcf ${samplename}_platypus.log $filterINDELS $exe_condition");
        print("\tSample $samplename variant calling submited with job_id $job_id\n");
    }
    else
    {
        print("\tSample $samplename variant calling already present, skipping it\n");
    }
    push(@namedvcffiles,"${samplename}.vcf");
}

my $cdeps=compile_dependencies(); #Dependencies for the N and SX files, for MITHE

### Calling with filters (Only cancer samples so far)
##############################################################

my $job_id;
my $actual_exe_conditions;
print("Filtered variant calling detection/execution:\n");
my $covNname;
my $BEDname;
my $fileexecondname;
my $filecovBname;
my @rep_deps; ##Dependencies of a specific group of samples with an specific exe_condition --> For covN
my %exe_deps; ##Dependencies of a specific group of samples with a specific exe_condition --> For MITHE (includes covN)
my $deps;
my $tempofile;
my $tag;
my $EXETOFILE;
my $exetcontent;
my $casename=basename($output_dir,(".csv"));
my $condition;
my $this_filtered_sample;
my $filtering_options;

foreach my $exe_condition (@exe_conditions)
{
    $BEDname="ALL$exe_condition.bed";
    $covNname="covN$sep_param$exe_condition.tsv";
    my @sample_vcfs;
    @rep_deps=();
    %exe_deps=();

    my @current_conditions;
    ##Current_conditions: for naming, includes info on exe and filtering conditions
    ##Filtering_conditions: needed to filter the files generated using the exe conditions
    for (my $i=0; $i<scalar @filtering_conditions;++$i)
    {
        $current_conditions[$i]="$exe_condition$sep_param$filtering_conditions[$i]";
    }

    for (my $nfile=0; $nfile<scalar @sample_bams; ++$nfile)
    {
        $samplename=sprintf("S%0*d", length scalar @sample_bams, $nfile);
        $fileexecondname="$samplename$sep_param$exe_condition.vcf";
        push(@sample_vcfs,$fileexecondname);
        push(@namedvcffiles,$fileexecondname);
        $filecovBname="$samplename$sep_param$exe_condition${sep_param}covBfiltering.tsv";
        if (! -f "$fileexecondname")
        {
            $bamfiles=$sample_bams[$nfile];
            if ($exe_condition eq "Default${sep_value}1")
            {
                $actual_exe_conditions="";
                #print("DEBUG: A: Default exe conditions\n");
            }
            else
            {
                $actual_exe_conditions=join(" ",split("$sep_value",join(" ",split("$sep_param",$exe_condition))));
                #print("DEBUG: A: Real exe_conditions\n");
            }

            $job_id=submit_job("$qsub $variant_calling_sh $bamfiles $fileexecondname $samplename$sep_param${exe_condition}_platypus.log $filterINDELS $actual_exe_conditions",\%exe_deps);
            push(@rep_deps,$job_id);
            $deps="$dep_prefix$sep_dep$job_id"; #For all the subsequent filters
            print("\tSample $samplename variant calling for conditions $exe_condition submited with job_id $job_id\n");

        }
        else
        {
            $deps="";
            print("\tSample $samplename variant calling for conditions $exe_condition already present, skipping it\n");
        }

        #WORKING HERE
        #TODO Improvement
        #These jobs are very quick, I could just calculate them with a single script? I cannot run them in this script since they depend on jobs submitted from it
        for (my $i=0; $i<scalar @current_conditions; ++$i)
        {
            $condition=$current_conditions[$i];
            $this_filtered_sample=vcf_refname("$samplename$sep_param$condition.vcf",$condition);
            $filtering_options=join(" ",split("$sep_value",join(" ",split("$sep_param",$filtering_conditions[$i]))));
            if (! -f $this_filtered_sample)
            {
                $job_id=submit_job("$qsub_noparallel $deps -J Filtering.$samplename.$i $helper_sh $vcf_filt_exe $filtering_options -i $fileexecondname -o $this_filtered_sample",\%exe_deps);

                print("\tFiltering $samplename to generate $this_filtered_sample submited with job_id $job_id\n");
            }
            else
            {
                print("$this_filtered_sample has been previously generated and it will be recycled. WARNING: this may generate inconsistencies if the number or order of parameters are changed between executions\n");
            }
        }
    }
    
    $deps=join($sep_dep,@rep_deps);
    
    if ($deps ne "")
    {
        $deps="$dep_prefix$sep_dep$deps";
    }
    
    #Generates a BED file with all the variants for this exe_cond. This will be used for covB and covN and it depends on @rep_deps
    if ( ! -f "$BEDname")
    {
        $job_id=submit_job("$qsub_noparallel $deps $covBED_sh $BEDname ".join(" ",@sample_vcfs),\%exe_deps);
        print("\tThe generation of the BED file for all sample variants for exe_conditions $exe_condition submitted with job_id $job_id\n");
        $deps="$dep_prefix$sep_dep$job_id";
    }
    else
    {
        $deps="";
        print("\tReusing previously generated $BEDname\n");
    }

    
    #covN depends on cdeps and cosvBED_sh
    my $covNdeps=$deps;
    if($cdeps ne "") ##There are cdeps
    {
        if ($covNdeps ne "") ##And there are previous deps, just concatenate
        {
            $covNdeps=$covNdeps.$sep_dep.$cdeps;
        }
        else ##No previous deps, we need to format them and concatenate cdeps
        {
            $covNdeps="$dep_prefix$sep_dep$cdeps";
        }
    }

    #Call covN variants just once
    if ( ! -f "$covNname")
    {
        $job_id=submit_job("$qsub_noparallel $covNdeps $covN_sh $BEDname $covNname $normal_bam",\%exe_deps);
        print("\tThe variant calling of N using the $BEDname BED file for exe_conditions $exe_condition submitted with job_id $job_id\n");
    }
    else
    {
        print("\tReusing previously generated $covNname\n");
    }


    #Call covB variants in each sample. These depend on covBED_sh
    for (my $nfile=0; $nfile<scalar @sample_bams; ++$nfile)
    {
        $samplename=sprintf("S%0*d", length scalar @sample_bams, $nfile);
        $filecovBname="$samplename$sep_param$exe_condition${sep_param}covBfiltering.tsv";
        if(! -f $filecovBname)
        {
            $job_id=submit_job("$qsub_noparallel $deps $covB_sh $BEDname $filecovBname $sample_bams[$nfile]",\%exe_deps);
            print("\tCalling variants on $samplename for covB calculations for exe_conditions $exe_condition submitted with job_id $job_id\n");
        }
        else
        {
            print("\tReusing previously generated $filecovBname\n");
        }

    }

    $deps=dependencies_string(\%exe_deps);

    if($cdeps ne "")
    {
        $deps=$deps.$sep_dep.$cdeps;
    }

    #Generating temprorary file with the exe-params of this exe-condition
    $tag=$exe_condition;
    $tag=~s/[^0-9]*//g; #Unique identifier generated with the numbers of the parameter values of this exe-condition
    open($EXETOFILE,">exeparams.$tag") or die "Problems generating the temp exeparams.$tag file\n";
    $exetcontent=$exe_condition;
    $exetcontent=~s/=/=$FS/g;
    $exetcontent=~s/$sep_param/\n/g;
    print($EXETOFILE $exetcontent);
    close($EXETOFILE) or die "Problems closing the temp exeparams.$tag file\n";
    
    #Generating name of the temporary output file for this exe-condition 
    $tempofile=$output_file;
    $tempofile=~s/\/[^\/]+$//g;

    for (my $i=0; $i<scalar @sample_bams; ++$i)
    {
        for (my $j=$i+1; $j<scalar @sample_bams; ++$j)
        {
            my $nameA=sprintf("S%0*d",length scalar @sample_bams, $i);
            my $nameB=sprintf("S%0*d", length scalar @sample_bams, $j);
            my $thistempofile="$tempofile/${casename}/${casename}_${nameA}_${nameB}.$tag.csv";
            $job_id=submit_job("$qsub $deps $mithe_time -J MITHE.${casename}_${nameA}_${nameB}.$tag $helper_sh $helper_pl -a $nameA -b $nameB -e exeparams.$tag -f $offile --NABfilt_cond_inputfile $onfile --covaltB_cond_inputfile $ocfile --popAF_cond_inputfile $opfile -o $thistempofile --output_vcf $output_vcf --output_list $output_list --comp $comp --n_cores $n_cores");
            print("\tThe filtering and analysis of the vcf files is being conducted with the job_id $job_id\n");
        }
    }
} 

$deps=dependencies_string();

for (my $i=0; $i<scalar @sample_bams; ++$i)
{
    for (my $j=$i+1; $j<scalar @sample_bams; ++$j)
    {
        my $nameA=sprintf("S%0*d",length scalar @sample_bams, $i);
        my $nameB=sprintf("S%0*d", length scalar @sample_bams, $j);
        my $thistempofile="$tempofile/${casename}/${casename}_${nameA}_${nameB}.csv";

        #TODO: This is not completely optimized, these jobs should depend only on a combination of samples, across all exe_params. Due to the order I am submitting the jobs prior to this, it would probably not save much time.
        $job_id=submit_job("$qsub_noparallel $deps $reduce_sh $thistempofile");
        print("The reduction of results of parallelized jobs for samples $i and $j, will be written in the file $tempofile/$casename/${casename}_${nameA}_${nameB}.csv by the job_id $job_id\n");   
    }
}

#I may need to add a final reduce per case here
my $finaljobids="";

unless ($output_vcf==0 || $comp != 2)
{
    $deps=dependencies_string();
    $job_id=submit_job("$qsub $deps $annotation_sh $n_cores"); 
    $finaljobids=$job_id;
    print("Variant annotation with annovar submitted with job_id $job_id\n");
    $job_id=submit_job("$qsub $deps $tstv_sh $output_dir");
    $finaljobids="$finaljobids$sep_dep$job_id";
    print("Ts/Tv statistics are being calculated in the job_id $job_id\n");

}
print("Jobs submitted!\n$finaljobids");

write_vcfrefname_dict("vcfdict.csv");

exit;

###################################################################################
###FUNCTIONS
###################################################################################

sub submit_job
{
    my $refhash=\%job_ids;
    my $command;

    if (scalar @_ == 2)
    {
        ($command,$refhash)=@_;
    }
    else
    {
        ($command)=@_;
    }

    my $job_id="";
    $job_id=`$command | $sed`;
    chomp($job_id);
    while ($job_id eq "")
    {
	    print("WARNING: Job submission has failed, trying again\n");
        sleep($sleep);
        $job_id=`$command | $sed`;
   		chomp($job_id);
    }
    ${$refhash}{$job_id}=1;
    return $job_id;
}

sub dependencies_string
{
    my $deps=compile_dependencies(@_);
    if ($deps ne "")
    {
        return "$dep_prefix$sep_dep$deps";
    }

    return $deps;
}

sub compile_dependencies
{
    my $option="";
    my $refhash=\%job_ids;
    if (scalar(@_) ==2)
    {
        ($refhash,$option)=@_;
    }
    elsif (scalar(@_) ==1)
    {
        ($refhash)=@_;
    }
    my $dependencies="";
    my $keep=0;
    if ($option =~/keep/i)
    {
        $keep=1;    
    }
    for my $key (keys %{$refhash})
    {
        $dependencies="$dependencies$key$sep_dep";
        $keep == 0 and delete(${$refhash}{$key});
    }
    chop($dependencies);
    return $dependencies;
}

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

#Parse vcf
###########################################
sub parse_vcf
{
	my ($vcf1_file)=@_;
	open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
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
			$key=$vcf1[$i];
			$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it!
            chomp($key);
			$hash{$key}=1;
			#print("DEBUG: New variant being hashed $key\n");
		}
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
                		#print("DEBUG: Name $name\n");
				$flag=1;
			}
			$key=$vcf1[$i];
			$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it!
            		chomp($key);
			$hash{$key}=1;
			#print("DEBUG: New variant being hashed $key\n");
		}
	}
	return (\%hash,$name);

}
#Compare two vcf files and generate statistics
##############################################
sub vcf_compare2
{
	my ($vcf1_file,$vcf2_file,$ref_statistics)=@_;
	
	open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
	open(my $VCF2,$vcf2_file) or die "The file $vcf2_file is not located in the output folder. Please, check if the variant caller has successfully finished";
	my @vcf1=<$VCF1>;
	my @vcf2=<$VCF2>;
	close($VCF1);
	close($VCF2);
	clean_vcfcontent(\@vcf1);
	clean_vcfcontent(\@vcf2);
	my %variants1=%{variants_to_hash(\@vcf1)};
	my %variants2=%{variants_to_hash(\@vcf2)};
	my %commonvariants;
	my $n1=scalar(keys %variants1);
	my $n2=scalar(keys %variants2);
	
	if ($n1<$n2) ##This could be improved in terms on speed following the strategy of the next subroutine. I'm not using this one any more, so I won't update it... at least so far.
	{
		foreach my $variant (keys %variants1)
		{
			if (exists $variants2{$variant})
			{
				$commonvariants{$variant}=1;
				#print("DEBUG: New common variant $variant");
			}
		}
	}
	else
	{
		foreach my $variant (keys %variants2)
                {
                        if (exists $variants1{$variant})
                        {
                                $commonvariants{$variant}=1;
				#print("DEBUG: New common variant $variant");
                        }
                }
	}
	@{ $ref_statistics }=($n1,$n2,scalar (keys %commonvariants));
	return \%commonvariants;

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
				last; ##No more pending comparisons
			}
	}	
    
    if ($n1<$n2)
	{
		print("WARNING: The number of variants in the unfiltered sample is smaller than the one of the filtered\n");		
	}
	
    my $n_common=scalar(keys %commonvariants);
    if($n2==0)
    {
        @{ $ref_statistics }=(0,$n_common); ##Stats= proportion of selected reads in the reference, number of selected variants

    }
    else
    {
	    @{ $ref_statistics }=($n_common/$n2,$n_common); ##Stats= proportion of selected reads in the reference, number of selected variants
    }
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
        #print("DEBUG: Variant $variant_to_remove ");
		if ($tag1==0 and exists($variants1{$variant_to_remove}))
		{
			delete($variants1{$variant_to_remove});
            #print("deleted in vcf1");
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
            #print("DEBUG: No more variants in vcf1\n");
		}	
		if($tag2==0 and scalar(keys %variants2)==0)
		{
			$tag2=1;
            #print("DEBUG: No more variants in vcf2\n");
		}
		
		if($tag1==1 and $tag2==1)
		{
            #print("DEBUG: No more variants\n");
			last;
		}
		
	}
	my $n_selvariants=scalar keys %variants1;
	my $n_filtvariants=scalar keys %variants2;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_statistics }=(0,$n_selvariants); 
    }
    else
    {   	
	    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
	}
    return (\%variants1,\%variants2);
}

#Combine the variants of two samples. Calculate the union of the two samples
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
    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
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
		$key=${$array}[$i];
		$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
		$hash{$key}=1;
		#print("DEBUG: New variant being hashed $key\n");
	}
	return \%hash;
}
# ############################################################
sub write_variant_list
{
    my ($ref_hash,$filename)=@_;
    open(my $FILE, ">$filename");
    print($FILE "#CHROM,POS\n");
    foreach my $variant (keys %{$ref_hash})
    {
        print($FILE join($OFS,split(",",$variant)),"\n");
    }
    close $FILE;
}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################
#($ref_common_variantsAfilt,"Afilt$sep_param${condition}.vcf","A.vcf","##Comment for the header");
sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment)=@_;
    open(my $OFILE, ">$filename");
    open(my $IFILE, "$vcf");
    my @icontent= <$IFILE>;
    close($IFILE);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;
     
    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for(my$i=0;$i< scalar @icontent; ++$i)
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
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            chomp($key);
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
    my $realfilename=$output_dir . "/" . $filter . $n_name . ".vcf";
    $dictrealnamevcf{$filename}=$realfilename;
    return $realfilename;
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
        print($FILE "$file$OFS$output_dir/$file\n");
    }
    close($FILE);
}
