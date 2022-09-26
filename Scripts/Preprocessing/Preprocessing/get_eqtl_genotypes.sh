#!/bin/bash

### Generate list of top eQTLs for each gene in each tissue
egenes_path=$1*.v8.egenes.txt.gz
GTEX_WGSv8="$2"
out_bed_path="$3"
out_vcf_path="$4"
out_cis_eqtl_path="$5"
GTEX_WGSv8_sorted="$6"

echo "$egenes_path"
echo $GTEX_WGSv8
echo $out_bed_path
echo $out_vcf_path
echo $out_cis_eqtl_path
echo $GTEX_WGSv8_sorted

#egenes_path: $GTEXv8/eqtl/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz
#GTEX_WGSv8
#out_bed_path: $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_cis_eQTLs.bed
#out_vcf_path: $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf
#out_cis_eqtl_path: $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs
echo "0"
#sort -V -k1,1 -k2,2 
##zcat $egenes_path | \
##        cut -f14-17 | grep -v variant_pos | sed 's/^X/23/g' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' | sort -V -k1,1 -k2,2n | \
##        sed 's/^23/X/g' | uniq > "$out_bed_path"

echo "out_bed_path"
ls $out_bed_path

wait
echo "1"


#if [[ 1 -eq 0 ]]
#then
### Extract these sites from the GTEx v8 VCF using bedtools
#tmp comment out
#zcat "$GTEX_WGSv8" | head -4000 | grep '#' > "$out_vcf_path"

wait


echo "1.5"
#zcat $GTEX_WGSv8  | \
#        cut -f14-17 | grep -v variant_pos | sed 's/^X/23/g' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' | sort -k 1,1 -k2,2n | \
#        sed 's/^23/X/g' | uniq > $GTEX_WGSv8_sorted

#bedtools intersect -a "$out_bed_path" -b "$GTEX_WGSv8_sorted" -loj -sorted | \
bedtools intersect -a "$out_bed_path" -b "$GTEX_WGSv8" -loj -sorted | \
        grep -v '\-1' | cut -f6- | sort -k1,1 -k2,2n | uniq >>  "$out_vcf_path"

wait
echo "2"
bgzip -f $out_vcf_path

wait
echo "3"
tabix -p vcf "$out_vcf_path.gz"

wait
### Convert the cis-eQTL genotypes in VCF format to the number of alternate alleles using VCFTools
echo "4"
vcftools --gzvcf "$out_vcf_path.gz" --out "$out_cis_eqtl_path" --012 --maf 0.01

#Rscript process_gtex_v8_cis_eqtl_genotypes.R

echo "DONE"
#fi
