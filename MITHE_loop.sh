#!/bin/bash

##Functions
function getFilePath {
    file=$1
    dir=$2
    if [[ $(echo $file | grep -c "^/") -eq 1 ]]
    then
        echo $(readlink -f $file)
    else
        echo $(readlink -f $dir/$file)
    fi
}

usage="MITHE Usage:\n-----------\n$0 directory torun_file exe_params filtering_params NAB_params covB_params popAF_params n_cores output_vcf output_list comprehensiveness filterINDELS [queue]\nThis script executes MITHE_control.pl for each case in a directory with its name. Then it integrates all the information in a file named results.csv and results_basictstv.csv.\nManifest file structure: output N_file S1_file S2_file ... Sn_file"

queue=""

echo -e "MITHE 0.1\n--------"

##Problems with number of arguments or the user is asking for help

if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
    echo -e $usage
    exit 0
fi

if (! ([[ $# -eq 12 ]] || [[ $# -eq 13 ]]))
then
    echo -e "ERROR: this script needs 12 or 13 arguments, but only $# have been provided.\n\n$usage"
    exit 1
fi

if [[ ! -d $1 ]]
then
    mkdir -p $1
fi

if [[ ! -s $2 ]]
then
    echo -e "ERROR: the manifest file cannot be located. Please, check the information below on how to call this program properly\n\n$usage"
    exit 1
fi

paramfilenames=("exe_params" "filtering_params" "NAB_params" "NAB2_params" "covB_params" "popAF_params")

for ifile in {3..7}
do
    if [[ ! -s ${!ifile} ]]
    then
        thefilename=$(( $ifile - 2 ))
        echo -e "ERROR: the parameter file for ${paramfilenames[!thefilename]}, was specified as ${!ifile} and cannot be located. Please, check the information below on how to call this program properly\n\n$usage"
        exit 1
    fi
done

basedir=$(dirname $0)

if [[ ! -s $basedir/config.txt ]]
then
    echo -e "ERROR: config.txt could not be sourced from $basedir/config.txt. Please, make sure you generate the config.txt file for your system following the installation tutorial in the README\n"
    exit 1
fi

source $basedir/config.txt

##Parsing ARGV
dir=$(readlink -f $1)
torun=$(readlink -f $2)
exe_params=$(readlink -f $3)
filtering_params=$(readlink -f $4)
NAB_params=$(readlink -f $5)
covB_params=$(readlink -f $6)
popAF_params=$(readlink -f $7)
n_cores=$8
output_vcf=${9}
output_list=${10}
comp=${11}
filterINDELS=${12}
queue=${13}

reldir=$(dirname $2)

#Currently, I am not using the dependencies, but I will need to when I fix the summarization step

dependency=""

while read line
do
    manifest=($line)
    output=${manifest[0]}
    control=${manifest[1]}
    samples=("${manifest[@]:1,2}")
    abscontrol=$(getFilePath $control $reldir)
    abssamplestext=""
    
    echo -e "Launching job for $output"
    if [[ ! -f $abscontrol ]]
    then
        echo -e "\tERROR: cannot open the file $abscontrol as the control sample for $output, according to the manifest file. MITHE will stop submitting jobs. You may want to stop the jobs submitted and empty the output directory before re-running this program. List of jobs already submitted:"
        echo $(echo $dependency | sed "s/${MITHE_SUBMIT_SEP}/ /g")
        exit 1
    fi

    for ifile in $(seq 0 $(( ${#samples[@]} - 1 )))
    do
        thisabsfile=$(getFilePath ${samples[$ifile]} $reldir)
        abssamplestext="$abssamplestext $thisabsfile"
        if [[ ! -f $thisabsfile ]]
        then
            echo -e "\tERROR: cannot open the file $thisabsfile as the $(( $ifile + 1 )) sample for $output, according to the manifest file. MITHE will stop submitting jobs. You may want to stop the jobs submitted and empty the output directory before re-running this program. List of jobs already submitted:"
            echo $(echo $dependency | sed "s/${MITHE_SUBMIT_SEP}/ /g")
            exit 1
        fi
    done

    abssamples=($abssamplestext)
    
    if [[ $queue == "" ]]
    then
        id=$($MITHE_INT/perl.sh $MITHE_INT/MITHE_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $abscontrol --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp --filterINDELS $filterINDELS ${abssamples[@]} | tee $dir/${output}.out | tail -n 1)
        #echo $MITHE_INT/perl.sh $MITHE_INT/MITHE_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $abscontrol --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp --filterINDELS $filterINDELS ${abssamples[@]}
    else
        id=$($MITHE_INT/perl.sh $MITHE_INT/MITHE_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $abscontrol --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp --queue $queue --filterINDELS $filterINDELS ${abssamples[@]} | tee $dir/${output}.out | tail -n 1) 
        #echo $MITHE_INT/perl.sh $MITHE_INT/MITHE_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $abscontrol --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp --queue $queue --filterINDELS $filterINDELS ${abssamples[@]}
    fi
    dependency="${dependency}${MITHE_SUBMIT_SEP}${id}"    
done < $torun

tstv=1

if [[ $comp -ne 2 ]] || [[ ${9} == 0 ]]
then
    tstv=0
fi

if [[ $queue == "" ]]
then
    $MITHE_SUBMIT_CMD ${MITHE_SUBMIT_DEP}$dependency $MITHE_MAX_MEM $MITHE_INT/summarizeResults.sh $dir $torun $tstv
else
    $MITHE_SUBMIT_CMD ${MITHE_SUBMIT_PAR}$queue ${MITHE_SUBMIT_DEP}$dependency $MITHE_MAX_MEM $MITHE_INT/summarizeResults.sh $dir $torun $tstv
fi
