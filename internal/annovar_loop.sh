#!/bin/bash
###The memory may generate problems

if [[ -n "$MITHE_MODULE_GNUPARALLEL" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_GNUPARALLEL") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_GNUPARALLEL" || module load "$MITHE_MODULE_GNUPARALLEL"
        else
            echo "ERROR: The module $MITHE_MODULE_GNUPARALLEL is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_GNUPARALLEL"
    fi
fi

if [[ -n "$MITHE_EXE_GNUPARALLEL" ]]
then
    eval "$MITHE_EXE_GNUPARALLEL"
fi

if [[ -n "$MITHE_MODULE_ANNOVAR" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_ANNOVAR") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_ANNOVAR" || module load "$MITHE_MODULE_ANNOVAR"
        else
            echo "ERROR: The module $MITHE_MODULE_ANNOVAR is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_ANNOVAR"
    fi 
fi

if [[ -n "$MITHE_EXE_ANNOVAR" ]]
then
    eval "$MITHE_EXE_ANNOVAR"
fi


if [[ -z "$MITHE_HUMANDB_DIR" ]]
then
    echo "ERROR: The environment variable MITHE_HUMANDB_DIR needs to point to the human genome for annovar to work\n";
    exit 1
else
    if [[ $(ls $MITHE_HUMANDB_DIR/*.fa | wc -l) -ne 1 ]]
    then
        echo "ERROR: The environment variable MITHE_HUMANDB_DIR is not properly set, there is no *.fa genome file there\n"
        exit 1
    fi
fi


n_cores=$1
#I have to get the number of cores from the arguments in $@ instead of from SLURM_JOB_CPUS_PER_NODE since if I increase the required memory I get more nodes that the ones I ask for!!!!!

annotate ()
{
    outputfile=$(basename $1)

    if [[ -f ${outputfile}.annotated.log ]]
    then
        echo "This file had already been anotated"
    else
        convert2annovar.pl --format vcf4 $1 > ${outputfile}.inputann
        annotate_variation.pl --geneanno --buildver hg19 --outfile ${outputfile}.annotated ${outputfile}.inputann $MITHE_HUMANDB_DIR
        rm -f ${outputfile}.inputann
    fi
}

export -f annotate

rdir=$PWD
for i in */vcfdict.csv
do
    dir=$(dirname $i)
    files=$(awk -v dir=$dir 'BEGIN{FS=","}{if($1~/^.*filt.*/){print $2}}' $i)
    cd $dir
    parallel --delay "0.2" -j $n_cores --joblog annovar.log annotate {1} ::: $files
    cd $rdir
done

