#!/bin/bash

usage="$0 bed outputname Nbam"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $3 ]]
then
    echo -e $usage
    exit 1
fi

bed=$1
out=$2
bam1=$3

if [[ -n "$MITHE_MODULE_GATK" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_GATK") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_GATK" || module load "$MITHE_MODULE_GATK"
        else
            echo "ERROR: The module $MITHE_MODULE_GATK is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_GATK"
    fi
fi

if [[ -n "$MITHE_EXE_GATK" ]]
then
    eval "$MITHE_EXE_GATK"
fi

if [[ -n "$MITHE_MODULE_BEDOPS" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_BEDOPS") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_BEDOPS" || module load "$MITHE_MODULE_BEDOPS"
        else
            echo "ERROR: The module $MITHE_MODULE_BEDOPS is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_BEDOPS"
    fi
fi

if [[ -n "$MITHE_EXE_BEDOPS" ]]
then
    eval "$MITHE_EXE_BEDOPS"
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


if [[ ! -f $out ]]
then
    name_out=$(echo $out | sed "s/.tsv//g")
    java -Xms512m -Xmx6G -jar $MITHE_GATKJAR -T UnifiedGenotyper -R $MITHE_HUMAN_GENOME -I $bam1 -o "$name_out.vcf" --intervals $bed --output_mode EMIT_ALL_SITES -glm BOTH -dcov 10000 > "$name_out.log" 2>&1
    cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
else
    echo "The file $out is present and will be reused"
fi
