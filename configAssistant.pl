#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

#ERRORS
my $SRC_ERR="ERROR: the specified base directory of MITHE does not contain MITHE sources\n";
my $FILE_ERR="ERROR: the file %s cannot be found\n";
my $DIR_ERR="ERROR: the directory %s cannot be found\n";
my $REQ_ERR="ERROR: this option is required and cannot be left blank\n";
my $MOD_ERR="ERROR: the module %s cannot be loaded\n";

#IO conf
my $MAX_HLEVEL=3;

print("\nMITHE Configurator assistant\n---------------------------\n");
my $dirname = abs_path(dirname(__FILE__));

print("Please, refer to the README file before running this assistant. Dependencies and supporting data should be resolved and prepared before running this script.\n\nConfiguring MITHE components and supporting data...\n");

#MITHE_HOME
unless (promptyn("\tSetting base directory of MITHE installation.\n\tAutomatically detected as \'$dirname\'. Is this correct?"))
{
    $dirname=prompt("\n\t\tInput the base directory of the MITHE installation: ");
}

if (! -f "$dirname/MITHE_loop.sh")
{
    die $SRC_ERR;
}

open(my $OUTFILE, ">$dirname/config.txt") or die "ERROR: impossible to generate the output file $dirname/config.txt\n";

writecomment("Main path and supporting data info",2);

writeparam("MITHE_HOME",$dirname);

#MITHE_INT
#I am automatically putting this one together

if (! -f "$dirname/internal/MITHE.pl")
{
    die $SRC_ERR;
}
writeparam("MITHE_INT","$dirname/internal/");

#MITHE_HUMANDB_DIR
promptDirAndSetVar("\tInput the base directory of annovar's human reference genome. This is the result of running annotate_variation.pl --downdb refGene DIRECTORY --build hg19, with an user-specified DIRECTORY location","","MITHE_HUMANDB_DIR");

#MITHE_HUMAN_GENOME
promptFileAndSetVar("\tInput the location of the human reference genome in fasta format","","MITHE_HUMAN_GENOME");

#MITHE_GNOMAD
promptFileAndSetVar("\tInput the location of the GNOMAD pre-processed population allele frequency information (i.e., the .vcf.bgz file). See MITHE manual for instruction son how to create this file","","MITHE_GNOMAD");

writecomment("Workload manager environment",1);

# Environment/SLURM related variables
print("Done\n\nSetting up the HPC environment in relation with your workload manager. In this section, defaults are provided for SLURM\n");

#MITHE_SUBMIT_CMD
promptAndSetVarDefault("\tIndicate the command to submit a job to your workload scheduler.","MITHE_SUBMIT_CMD","sbatch");

#MITHE_SUBMIT_SED
promptAndSetVarDefault("\tIndicate a command that, when executed piping in the output of your submit command, will extract the job id of that job.","MITHE_SUBMIT_SED",'sed "s/Submitted batch job \(.*\)/\1/"');

#MITHE_SUBMIT_PAR
promptAndSetVarDefault("\tIndicate the argument to add to the submit command to indicate the partition/queue specified when submitting a job.","MITHE_SUBMIT_PAR",'--partition=');

#MITHE_SUBMIT_MUL
promptAndSetVarDefault("\tIndicate the arguments needed to add to the submit command to submit a multi-threaded job. The number of threads is specified differently, and will be directly appended to this variable.","MITHE_SUBMIT_MUL",'-N 1 -n 1 -c ');

#MITHE_SUBMIT_DEP
promptAndSetVarDefault("\tIndicate the arguments needed to add to the submit command to indicate a dependency.","MITHE_SUBMIT_DEP",'--dependency=afterok');

#MITHE_SUBMIT_SEP
promptAndSetVarDefault("\tIndicate the character used to separate job-ids to indicate multiple dependencies.","MITHE_SUBMIT_SEP",':');

#MITHE_NCPUS_VAR
promptAndSetVarDefault("\tIndicate the environment variable that will indicate the number of CPU threads available for multi-threaded jobs (within the job).","MITHE_NCPUS_VAR",'SLURM_JOB_CPUS_PER_NODE');

#MITHE_MAX_TIME
promptAndSetVarDefault("\tIndicate the argument needed to add the maximum time limit to submit long jobs. This is not required, but important to be used for optimization runs.","MITHE_MAX_TIME",'-t 4-00');

#MITHE_MAX_MEM
promptAndSetVarDefault("\tIndicate the argument needed to add the maximum RAM usage to submit heaby jobs. This is not required, but important to be used for runs with vcf output. No condition should require more than 16G.","MITHE_MAX_MEM",'--mem=16G');

# Modules

writecomment("Software environment/modules",1);

my @modules=("PERL", "ANNOVAR", "PLATYPUS", "VCFTOOLS", "GNUPARALLEL", "BEDOPS", "SAMTOOLS");

my @modulesJ=("GATK", "SNPSIFT");
my @modulesJVars=("MITHE_GATKJAR", "MITHE_SNPSIFTJAR");

print("Done.\n\nSetting up the software environment. For each program, you will have the possibility of setting a module name to load and some code to execute before using it\n");

#MITHE_MOD_ISA
my $isAvail=`module is-avail 2>&1 | grep -c ERROR`;
if ($isAvail > 0)
{
    print "\tModule environments not compatible with is-avail detected\n";
    writeparam("MITHE_MOD_ISA",0);
}
else
{
    print "\tModule environments compatible with is-avail detected\n";
    writeparam("MITHE_MOD_ISA",1);
}

foreach my $module (@modules)
{
    print "\n\tProgram: $module\n";
    writecomment("Program $module",2);
    promptAndSetVarModule($module,"MITHE_MODULE_${module}");
    promptAndSetVar("\t\tInsert the code to execute before running $module\n","MITHE_EXE_${module}");
}

for (my $imodule=0; $imodule < scalar @modulesJ; ++$imodule)
{
    my $module=$modulesJ[$imodule];
    writecomment("Program $module",2);
    print "\n\tProgram: $module\n";
    promptAndSetVarModule($module,"MITHE_MODULE_${module}");
    promptAndSetVarReqVar($module,"$modulesJVars[$imodule] with ${module}\'s jar","MITHE_EXE_${module}");
}

print("\nDone. Configuration finished. Enjoy using MITHE\n");

close($OUTFILE);

#SUBROUTINES
#IO
sub prompt
{
    my ($query) = @_;
    local $| = 1; # activate autoflush to immediately show the prompt
    print $query;
    chomp(my $answer = <STDIN>);
    return $answer;
}

sub promptyn
{
    my ($query) = @_;
    my $answer=prompt($query." [y/n]: ");
    return lc($answer) eq 'y';
}

sub writeparam
{
    my ($var,$value)=@_;
    if ($value eq '')
    {
        print($OUTFILE "#export $var='$value'\n");
    }
    else
    {
        print($OUTFILE "export $var='$value'\n");
    }
}

sub writecomment
{
    my ($value,$level)=@_;
    my $head='#' x ($MAX_HLEVEL-$level);
    print($OUTFILE "\n$head$value\n");
}

sub promptFile
{
    my ($query,$file) = @_;
    my $answer=prompt($query);
    my @cdirs=($answer);
    $file ne '' and push(@cdirs,$file);
    $file=join("/",@cdirs);
    if (! -f $file)
    {
        die sprintf($FILE_ERR,$file);
    }
    return abs_path($file);
}

sub promptDir
{
    my ($query,$dir) = @_;
    my $answer=prompt($query);
    my @cdirs=($answer);
    $dir ne '' and push(@cdirs,$dir);
    $dir=join("/",@cdirs);
    if (! -d $dir)
    {
        die sprintf($DIR_ERR,$dir);
    }
    return abs_path($dir);
}

sub promptFileAndSetVar
{
    my($query,$file,$var) = @_;
    my $answer = promptFile("$query: ",$file);
    writeparam($var,$answer);
    return $answer;
}

sub promptDirAndSetVar
{
    my($query,$file,$var) = @_;
    my $answer = promptDir("$query: ",$file);
    writeparam($var,$answer);
    return $answer;
}

sub promptAndSetVarDefault
{
    my($query,$var,$default) = @_;
    my $answer = prompt("$query Default: \'$default\' : ");
    if ($answer eq "")
    {
        writeparam($var,$default);
        return $default;
    }
    else
    {
        writeparam($var,$answer); 
        return $answer;
    }
}

sub promptAndSetVar
{
    my($query,$var) = @_;
    my $answer = prompt($query);
    writeparam($var,$answer);
    return $answer;
}

sub promptAndSetVarModule
{
    my($module,$var) = @_;
    my $answer = promptAndSetVar("\t\tInsert the module name to load $module, or leave blank if no module needs to be loaded: ",$var);
    my $test=`module load $answer 2>&1 | grep -c ERROR`;
    $test > 0 ? die sprintf($MOD_ERR,$var) : print ("\t\t$var properly set\n");
    return $answer;
}

sub promptAndSetVarReqVar
{
    my($module,$reqvar,$var) = @_;
    return promptAndSetVar("\t\tInsert the code to execute before running $module.\n\t\t***WARNING***: this requires setting the variable $reqvar. : ",$var);
}

sub promptAndSetRequired
{
    my($query,$var,$default) = @_;
    my $answer = prompt($query);
    if ($answer eq "")
    {
        writeparam($var,$default);
        return $default;
    }
    else
    {
        writeparam($var,$answer); 
        return $answer;
    }
}
