#!/bin/bash

usage="\n$0 main_results_directory filelist \nThis script parses the output of MITHE_control.pl and generates a list of variants for each patient taking into consideration the information of all their samples. It generates a VCF file per sample, a sample file with information of number of variants in each sample, variants file with the presence of each in the different samples and classification of clonal, subclonal, and private, and a CSV file with the genotype of each sample for each variant. FASTA and CloneFinder files are quick and dirty solutions at this point and need to be improved\n\nOptions:\n--------\n\t--min_coverage: minimum number of reads to recall a variant. Below that threshold, a variant in a sample will be considered ?/?\n\t--min_alternative_reads: minimum number of reads to call a non-called variant present in another sample.\n\t--filt_file: file with a list of variants to keep, the rest are not taken into consideration. If this file is not used, all of them are taken into consideration\n\nInfo:\n\tRecall process: if a variant has not been called in the sample is considered ?/? if it does not have enough coverage (min_coverage). If it does have enough coverage, it will be considered 0/0 if the number of alternative reads is n<min_alternative_reads. Otherwise, the variant will be recalled. In the non-recalling output, the latter would be ?/? (never call something that has not been called, unless it is a reference 0/0)\n"

if [[ $# -lt 2 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
    echo -e $usage
    exit 1
fi

dir=$1
torun=$2

shift 2

basedir=$(dirname $0)

if [[ ! -s $basedir/config.txt ]]
then
    echo -e "ERROR: config.txt could not be sourced from $basedir/config.txt. Please, make sure you generate the config.txt file for your system following the installation tutorial in the README\n"
    exit 1
fi

source $basedir/config.txt

while read line
do
    manifest=($line)
    output=${manifest[0]}
    control=${manifest[1]}
    samples=("${manifest[@]:1,2}")
    nsamples=${#samples[@]}

   
    echo "Submitting job to obtain the final variant information from patient $output:"
    echo $MITHE_SUBMIT_CMD $MITHE_MAX_MEM $MITHE_INT/perl.sh $MITHE_INT/getVariants.pl -d $dir/$output -n $nsamples -o "variants" $@
    
done < $torun
