#!/bin/bash

if [[ -n "$MITHE_MODULE_PERL" ]]
then
    if [[ "$MITHE_MOD_ISA" -eq 1 ]]
    then
        if [[ $(module is-avail "$MITHE_MODULE_PERL") -eq 1 ]]
        then
            module is-loaded "$MITHE_MODULE_PERL" || module load "$MITHE_MODULE_PERL"
        else
            echo "ERROR: The module $MITHE_MODULE_PERL is not available"
            exit 1
        fi
    else
        module load "$MITHE_MODULE_PERL"
    fi
fi

if [[ -n "$MITHE_EXE_PERL" ]]
then
    eval "$MITHE_EXE_PERL"
fi

perl $@
