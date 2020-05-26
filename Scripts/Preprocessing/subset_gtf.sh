#!/bin/bash

gtf=$1 #input: general gencode gtf file
out_a=$2 #output: autosomal only file
out_x=$3 #output: x only file


#take only autosomal genes
#take only lincRNA and protein coding genes
#outputs to out_a
less -S $gtf | awk -F"[\t;]" '{if($1~/chr[0-9]{1,2}/){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' > ${out_a}

#take only x chromsome genes
#take only lincRNA and protein coding genes
#outputs to out_x
less -S $gtf | awk -F"[\t;]" '{if($1=="chrX"){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' > ${out_x}
