#!/bin/bash

# Author: Emily Tsang

echo "running vcf2bedfiles_processVCF.sh"

# uses vcftools to extract relevant information from vcf files
# helper script for vcf2bedfiles

minmaf=0.000001 # would need to change this if using a huge dataset
maxmaf=0.25

# take as input:
# * vcf file
# * file with individuals to include
# * output directory

if [ $# -ne 3 ]; then
    echo "usage: vcf2bedfiles_helper_provessVCF.sh <vcf> <individuals file> <output prefix>"
    exit
fi

vcf=$1
indincl=$2 # this should be only european american individuals
outprefix=$3

# individuals to include from the allele frequency calculation
# generate file name from the indincl file name
keeplist=${indincl%_ids.txt}_VCFids.txt # as above, but IDs as they appear in the snp/indel VCF

# get their ids as shown in the vcf file
zcat $vcf | awk '{if(substr($1,1,2)!="##"){print; exit}}' | \
awk -v indincl=$indincl 'BEGIN{
    while((getline<indincl)>0){
        IDs[$1]
    }
}{
    for(i=10;i<=NF;i++){
        if($i in IDs){print $i}
    }
}' > $keeplist


echo "KEEP DONE"
date

if [ ! -f ${outprefix}_SNPs ]
then
	# SNPs
	vcftools --gzvcf $vcf --out ${outprefix}_SNPs --maf $minmaf --max-maf $maxmaf --remove-filtered-all --keep $keeplist --remove-indels --max-missing-count 10 --freq &
	vcftools --gzvcf $vcf --out ${outprefix}_SNPs --maf $minmaf --max-maf $maxmaf --remove-filtered-all --keep $keeplist --remove-indels --max-missing-count 10 --extract-FORMAT-info GT &
fi

if [ ! -f ${outprefix}_indels ]
then
	# indels
	vcftools --gzvcf $vcf --out ${outprefix}_indels --maf $minmaf --max-maf $maxmaf --remove-filtered-all --keep $keeplist --keep-only-indels --max-missing-count 10 --freq &
	vcftools --gzvcf $vcf --out ${outprefix}_indels --maf $minmaf --max-maf $maxmaf --remove-filtered-all --keep $keeplist --keep-only-indels --max-missing-count 10 --extract-FORMAT-info GT
fi
wait
