#!/bin/bash

gtf=$1 #GTF FILE
out_a=$2 #Out file: autosomal
out_x=$3 #out file: x
#gtf=/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf

less -S $gtf | awk -F"[\t;]" '{if($1~/chr[0-9]{1,2}/){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' > ${out_a}
less -S $gtf | awk -F"[\t;]" '{if($1=="chrX"){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' > ${out_x}
