#!/bin/bash

set -o nounset -o errexit -o pipefail

gtf=${1}
out_a=${2}
out_x=${3}
bed_a=${4}
bed_x=${5}
outdir=${6}
#gtf=/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf


#print gtf files which are 2 cols: 1=gene_name, 2=[protein_coding/lincRNA] sorted by gene name for aut and x separately
less -S $gtf | awk -F"[\t;]" '{if($1~/chr[0-9]{1,2}/){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' | sort -k 1 > ${out_a}
less -S $gtf | awk -F"[\t;]" '{if($1=="chrX"){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' | sort -k 1 > ${out_x}


# gene bed file
cat ${gtf} | awk 'BEGIN{OFS="\t"}{
	if (substr($1,1,1)=="#") {next};
	if ($3!="gene") {next};
	name=substr($10,2,length($10)-3);
	print $1,$4,$5,name
}' > ${outdir}/gtf.bed #${gtexprefix}.bed

# gene bed file with 10kb added on either side
#cat ${gtexprefix}.bed | awk 'BEGIN{OFS="\t"}{
cat ${outdir}/gtf.bed | awk 'BEGIN{OFS="\t"}{
	low = $2-10000;
	if (low < 0) {low = 0};
	print $1,low,$3+10000,$4
}' > ${outdir}/gtf_padded10kb.bed


rm $outdir/gtf.bed

echo "sort and join"
sort -k 4 ${outdir}/gtf_padded10kb.bed > tmp.txt && mv tmp.txt ${outdir}/gtf_padded10kb.bed
join --nocheck-order -j1 1 -j2 4 ${out_a} ${outdir}/gtf_padded10kb.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > ${bed_a} #${outdir}/gtf_aut_padded10kb.bed
join --nocheck-order -j1 1 -j2 4 ${out_x} ${outdir}/gtf_padded10kb.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > ${bed_x} #${outdir}/gtf_x_padded10kb.bed
echo "DONE"
