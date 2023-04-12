#!/bin/bash
# This script can be used to screen candidate associations

## check same correlations test,spear and lasso results are needed。projects_list.txt includes 32 cancer types
cat projects_list.txt | while read pro
do
	#sort
	sed -i "s/^y\t/111-gene\t/g" ${pro}.lasso.projection.test.result.tsv
	sort -k 1 ${pro}.lasso.projection.test.result.tsv > ${pro}.sort.lasso.projection.test.result.tsv
	cut -f 1 ${pro}.lasso.projection.test.result.tsv > ${pro}-gene
	sort ${pro}-gene | uniq > ${pro}-gene-uniq
	rm -rf ${pro}-gene
	sed -i "/111-gene/d" ${pro}-gene-uniq
done


## check and merge features, spear and lasso results are needed
for filename in *.lasso.projection.test.result.tsv
do
	base=$(basename $filename .lasso.projection.test.result.tsv)

	cut -f 1-2 ${base}.lasso.projection.test.result.tsv > ${base}-lasso-name
	cut -f 1-2 ${base}.spear.correlation.test.result.tsv > ${base}-spear-name

	## make sure same feature in every row
	diff ${base}-lasso-name ${base}-spear-name >> diff.txt
	rm -rf *-name
	paste ${base}.spear.correlation.test.result.tsv ${base}.lasso.projection.test.result.tsv > ${base}-corr-lasso.tsv
done

## screen associations
cat projects_list.txt | while read i
 do
	# lasso
	awk '$7<0.1 && $14<0.1 && $20=="TRUE"' $i-corr-lasso.tsv > $i-corr-lasso-filter.tsv

	# spearman
	awk '$7<0.1' $i-corr-lasso.tsv > $i-corr-filter.tsv

	# spearman and lasso
	awk '$14<0.1 && $20=="TRUE"' $i-corr-lasso.tsv > $i-lasso-filter.tsv
done

## add cancer name
for filename in *-corr-lasso.tsv; do 
	base=$(basename $filename -corr-lasso.tsv)
	sed -i "s/k__/${base}\tk__/g" ${base}-corr-lasso.tsv
done


