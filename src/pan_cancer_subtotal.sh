#!/bin/bash
## find pan-cancer associations
for filename in *-corr-lasso-filter.tsv; do 
	base=$(basename $filename -corr-lasso-filter.tsv)
	cut -f1-2 ${base}-corr-lasso-filter.tsv > ${base}-name
	sed -i "s/\t/--/g" ${base}-name
done

cat *-gene-taxa > all-gene-taxa
sort all-gene-taxa > sort-all-gene-taxa

cat sort-all-gene-taxa|awk -F '\t' '{
	x[$1]++;
	} 
	END
		{
		for(i in x) print(i ":" x[i])
		}' > Sub_total