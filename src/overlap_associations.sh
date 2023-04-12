#!/bin/bash

cut -f 1-2 $i-corr-lasso-filter.tsv > *-cpg-taxa

## Cut -f 1-2 to get features and taxa
for filename in *-cpg-taxa; do 
    base=$(basename $filename -cpg-taxa)
    awk 'NR==FNR{huty[$1]=$2;next}NR>FNR{
        if($1 in huty)
        {
            print $0 "\t" huty[$1]
            } 
            else 
            {
                print $0 "\t" "None"
        }
    }' cpg-meta ${base}-cpg-taxa > ${base}.vlookup.tsv
    
    cut -f 2,3 ${base}.vlookup.tsv > ${base}-m-gene
    cut -f 2 ${base}.vlookup.tsv > ${base}-m
    i=1
    y=1
    cat ${base}-m | while read pro
        do
        sed -i "$i s/;/;${pro}\t/g" ${base}-m-gene
        i=$(($i+$y))
        done
done

## find correlations
for filename in *-m-gene; do 
    base=$(basename $filename -m-gene)
    sed -i "s/;/\n/g" ${base}-m-gene
    cut -f 2 ${base}-m-gene > ${base}-g
    cut -f 1 ${base}-m-gene > ${base}-m
    paste ${base}-g ${base}-m > ${base}-g-m
    sed -i "s/\t/--/g" ${base}-g-m
done

## find overlap
for filename in *-gene-taxa; do 
    base=$(basename $filename -gene-taxa)
    grep -Ff ${base}-gene-taxa  ${base}-g-m > ${base}.overlap
    sort ${base}.overlap | uniq > ${base}.overlap.uniq
    wc -l ${base}.overlap.uniq >> wcl
done