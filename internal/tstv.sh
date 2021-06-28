#!/bin/bash

if [[ -n "$MITHE_MODULE_VCFTOOLS" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_VCFTOOLS") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_VCFTOOLS" || module load "$MITHE_MODULE_VCFTOOLS"
        else
            echo "ERROR: The module $MITHE_MODULE_VCFTOOLS is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_VCFTOOLS"
    fi
fi

if [[ -n "$MITHE_EXE_VCFTOOLS" ]]
then
    eval "$MITHE_EXE_VCFTOOLS"
fi

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

if [[ $# -ne 1 ]] || [[ ! -d $1 ]]
then
    echo "$0 input_folder"
    exit
fi

dir=$1

cd $dir

fun() {
    
    name=$(echo $1 | sed "s/.vcf//")
    file=$2
    realname=$(echo $2 | sed "s/.vcf//")
    ratio=$(vcftools --vcf $file --TsTv-summary --out $realname 2>&1 | sed -n -e "/Ts\/Tv ratio/p" | sed "s/Ts\/Tv ratio: //")
    echo "$name,$ratio"
}

export -f fun

rdir=$PWD
for i in */vcfdict.csv
do
    thisdir=$(dirname $i)
    cd $thisdir
    echo "File,tstv" > tstv.csv
    parallel --colsep ',' -j ${!MITHE_NCPUS_VAR} --joblog tstv.log fun {1} {2} :::: vcfdict.csv >> tstv.csv
    cd $rdir
done
