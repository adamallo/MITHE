#!/bin/bash

usage="$0 outputname\n"
sufix=".csv"


if [[ $# -ne 1 ]]
then
    echo -e $usage
    exit 1
fi

input=$1
dir=$(dirname $input)
name=$(basename $input | sed "s/$sufix$//" | sed "s/\.[^.]*$//")
pref=$dir/$name
output="$pref$sufix"
withheader=0

#cat outs
for i in ${pref}.[0-9]*$sufix
do
    if [[ $withheader -eq 0 ]]
    then
        cat $i > $output
        withheader=1
    else
        tail -n+2 $i >> $output
    fi
    rm -f $i
done

cat $pref/vcfdict* | sort | uniq > $pref/vcfdict.csv
cat $pref/listdict* | sort | uniq > $pref/listdict.csv
rm -f $pref/vcfdict.*.*
rm -f $pref/listdict.*.*
