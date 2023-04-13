## GI cancer list: COAD，CHOL，ESCA，LIHC，PAAD，READ，STAD

cat GI_cancer_list | while read pro
do
    grep "	$pro	" protein-all >> protein-GI
done
cut -f 1,3 protein-GI > protein-GI-m-g
sed -i "s/\t/-/g" protein-GI-m-g

cat protein-GI-m-g|awk '{
    x[$1]++;
    } 
    END
        {
        for(i in x) print(i ":" x[i])
        }' > protein-GI-pancancer

wc -l protein-GI-pancancer > wcl