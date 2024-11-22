#!/bin/bash
#SBATCH --job-name=feature_gt
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
##SBATCH --time=06:59:00
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/featuregt_sexgtex_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/featuregt_sexgtex_%j.err"

module load R vcftools tabix bedtools python3;

#called from vcf2bedfiles_processVCFtoolsOutput.sh
i=${1}
chr=${2}
dir=${3}
af=${4}
outdir=${5}
TYPE=${6}
gt=${7}

#sample=`head -n1 $gt | cut -f$i`
sample="GTEX-1R9PO"
echo "i: $i , SAMPLE: $sample, CHR: $chr, TYPE: $TYPE"
shopt -s expand_aliases

# skip first line of file
alias skip_header="tail -n+2"

# to remove missing genotypes
alias filter_genotypes="grep -v '\./\.' | grep -v '\.|\.'"

# to turn genotypes into 0,1,2 (the number of non-reference alleles)
# 0/0 -> 0; 0/*,*/0 -> 1; */* -> 2 (where * is any number >0)
# now considers phased genotypes (i.e. 0|0 -> 0, ...)
alias clean_genotypes="sed 's%0\(/\||\)0%0%' | sed 's%[1-9]\(/\||\)0%1%' | sed 's%0\(/\||\)[1-9]%1%' | sed 's%[1-9]\(/\||\)[1-9]%2%'"

if [ "$TYPE" = "SNPs" ]; then
  cat $gt | cut -f1,2,$i | skip_header | filter_genotypes | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3}' | sort -k1,1 -k2,2n | \
  bedtools intersect -sorted -wa -wb -a stdin -b ${dir}/${chr}_$af | \
  awk 'BEGIN{OFS="\t"}{split($4,geno,"[/|]"); i=10+geno[1]; j=10+geno[2];print $1,$2,$3,$4,$8,$9,$10,$11}' | \
  clean_genotypes | \
  awk 'BEGIN{OFS="\t"}{if($5>0.25){next}; if($4==1 || $4==$6){print $1,$2,$3,$5,$4,$7,$8}}' \
  > ${outdir}/${chr}_${sample}_${TYPE}.bed
else
cat $gt | cut -f1,2,$i | skip_header | filter_genotypes | clean_genotypes | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3}' | sort -k1,1 -k2,2n | \
  bedtools intersect -sorted -wa -wb -a stdin -b ${dir}/${chr}_$af | \
  awk 'BEGIN{OFS="\t"}{if($8>0.25){next}; if($4==1 || $4==$9){print $1,$2,$3,$8,$4}}' \
  > ${outdir}/${chr}_${sample}_${TYPE}.bed
fi  

echo "DONE"
