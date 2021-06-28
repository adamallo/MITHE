#!/bin/bash
#SBATCH --mem 16000

usage="\n$0 directory torun_file tstv\n\ntorun_file structure: Manifest file structure: output N_file S1_file S2_file ... Sn_file\n-------------------------------------------------\n
\n
This script postprocess the output of MITHE_control.pl for each sample in a directory with its name, integrating all the information in a file named results.csv. Results results_basictstv.csv is generated only if tstv==1.\n"

if [[ $# -ne 3 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    tstv=${3}
fi

cat $torun | while read line
do
    manifest=($line)
    output=${manifest[0]}
    control=${manifest[1]}
    samples=("${manifest[@]:1,2}")
    
    outputHeader=1

    for ifile in $(seq 0 $(( ${#samples[@]} - 1 )))
    do
        for jfile in $(seq $(( $ifile + 1 )) $(( ${#samples[@]} - 1 )))
        do
            comparison=$dir/${output}/${output}_S${ifile}_S${jfile}
            comparisonName=${output}_S${ifile}_S${jfile}

            if [[ ! -d $comparison ]]
            then
                echo "ERROR: folder $comparison not found!"
                exit 1
            fi
            cat "${comparison}.csv" | awk 'BEGIN{OFS=",";FS=","}{out="";for (i=2;i<=NF;++i){out=out $i OFS};out=substr(out,0,length(out)-1);print out}' > $comparison/csvtomerge.temp

            if [[ $tstv -eq 1 ]]
            then
                $MITHE_INT/perl.sh $MITHE_INT/tstv_pretabulate_heter.pl $comparison/tstv.csv $comparison/tstv.tomerge
                $MITHE_INT/perl.sh $MITHE_INT/merge_csv.pl $comparison/csvtomerge.temp $comparison/tstv.tomerge $comparison/${comparisonName}_withtstv.temp
                rm -f $comparison/csvtomerge.temp
                $MITHE_INT/perl.sh $MITHE_INT/tabulate_results.pl -i $comparison/${comparisonName}_withtstv.temp -o $comparison/${comparisonName}_withtstv_tab.temp
                cat $comparison/${comparisonName}_withtstv_tab.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Comparison\",\$0} else {print \"$comparisonName\",\$0}}" > ${comparison}/${comparisonName}_withtstv.csv
                rm -f $comparison/${comparisonName}_withtstv.temp $comparison/${comparisonName}_withtstv_tab.temp
        
                #BasicTsTv postprocessing and general tstv data gathering
                a=$(cat $comparison/tstv.csv | sed -n "/^A,/p" | sed "s/^A,//")
                b=$(cat $comparison/tstv.csv | sed -n "/^B,/p" | sed "s/^B,//")
                n=$(cat $comparison/tstv.csv | sed -n "/^N,/p" | sed "s/^N,//")
            else 
                $MITHE_INT/perl.sh $MITHE_INT/tabulate_results.pl -i $comparison/csvtomerge.temp -o $comparison/csvtomerge_tab.temp
                cat $comparison/csvtomerge_tab.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Comparison\",\$0} else {print \"$comparisonName\",\$0}}" > ${comparison}/${comparisonName}_withtstv.csv
                rm -f $comparison/csvtomerge.temp $comparison/csvtomerge_tab.temp
            fi
            if [[ $outputHeader -ne 1 ]]
            then
                tail -n +2 ${comparison}/${comparisonName}_withtstv.csv >> ${dir}/${output}_results.csv
            else
                cat ${comparison}/${comparisonName}_withtstv.csv > ${dir}/${output}_results.csv
                if [[ $tstv -eq 1 ]]
                then
                    echo "Comparison,TsTv_A,TsTv_B,TsTv_N" > ${dir}/${output}_results_basictstv.csv
                fi
                outputHeader=0
            fi
            if [[ $tstv -eq 1 ]]
            then
                echo "$comparisonName,$a,$b,$n" >> ${dir}/${output}_results_basictstv.csv
            else
                rm -f ${comparison}/${comparisonName}_withtstv.csv
            fi
        done
    done
done
