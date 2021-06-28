#!/bin/bash

if [[ -n "$MITHE_MODULE_SNPSIFT" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_SNPSIFT") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_SNPSIFT" || module load "$MITHE_MODULE_SNPSIFT"
        else
            echo "ERROR: The module $MITHE_MODULE_SNPSIFT is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_SNPSIFT"
    fi
fi

if [[ -n "$MITHE_EXE_SNPSIFT" ]]
then
    eval "$MITHE_EXE_SNPSIFT"
fi

$MITHE_INT/perl.sh $MITHE_INT/vcf_filtering.pl $@ 
