#!/bin/bash

if [[ -n "$MITHE_MODULE_PLATYPUS" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_PLATYPUS") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_PLATYPUS" || module load "$MITHE_MODULE_PLATYPUS"
        else
            echo "ERROR: The module $MITHE_MODULE_PLATYPUS is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_PLATYPUS"
    fi
fi

if [[ -n "$MITHE_EXE_PLATYPUS" ]]
then
    eval "$MITHE_EXE_PLATYPUS"
fi

if [[ -z "$MITHE_HUMAN_GENOME" ]]
then
    echo "ERROR: The environment variable MITHE_HUMAN_GENOME needs to point to the fasta file with the human reference genome\n";
    exit 1
else
    if [[ ! -f "$MITHE_HUMAN_GENOME" ]]
    then
        echo "ERROR: The environment variable MITHE_HUMAN_GENOME is not set properly\n"
        exit 1
    fi
fi

bamfiles=$1
output=$2
logfilename=$3
filterINDELS=$4

if [[ ! -s ${bamfiles}.bai ]]
then
    if [[ ! -s $bamfiles ]]
    then
        echo -e "ERROR: the file $bamfiles cannot be found\n"
    else
        echo -e "WARNING: the file $bamfiles is not indexed. MITHE will try to use samtools to index it before running platypus\n"
        
        ##SAMTOOLS ENVIRONMENT
        if [[ -n "$MITHE_MODULE_SAMTOOLS" ]]
        then
            if [[ "$MITHE_MOD_ISA" -eq 1 ]]
            then
                if [[ $(module is-avail "$MITHE_MODULE_SAMTOOLS") -eq 1 ]]
                then
                    module is-loaded "$MITHE_MODULE_SAMTOOLS" || module load "$MITHE_MODULE_SAMTOOLS"
                else
                    echo "ERROR: The module $MITHE_MODULE_SAMTOOLS is not available"
                    exit 1
                fi
            else
                module load "$MITHE_MODULE_SAMTOOLS"
            fi
        fi
        
        if [[ -n "$MITHE_EXE_SAMTOOLS" ]]
        then
            eval "$MITHE_EXE_SAMTOOLS"
        fi
        ###
    
        samtools index $bamfiles
    fi
fi

if [[ $# -ge 4 ]]
then
    shift 4
    platypus callVariants --nCPU=${!MITHE_NCPUS_VAR} --bamFiles=$bamfiles --refFile=${MITHE_HUMAN_GENOME} --output=$output --logFileName=$logfilename $@ 
else
    echo "Error, the number of arguments for this script is not appropriate"
fi

mv $output ${output}_bkp_multisnv
$MITHE_INT/perl.sh $MITHE_INT/separateMultipleSNVPlatypus.pl -i ${output}_bkp_multisnv -o $output -f $filterINDELS
