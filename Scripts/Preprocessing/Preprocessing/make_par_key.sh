#bedtools intersect -a /oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf -b Output/preprocessing_v8/par_regions.bed -wb | awk -F'"' '{print $2"\t"$NF}' | awk -F"\t" '{print $1"\t"$NF}' > Output/preprocessing_v8/par_table.txt
bedtools intersect -a /oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf -b Output/preprocessing_v8/evolutionary_strata.txt -wb | awk -F'"' '{print $2"\t"$NF}' | awk -F"\t" '{print $1"\t"$NF}' > Output/preprocessing_v8/strata_table.txt

