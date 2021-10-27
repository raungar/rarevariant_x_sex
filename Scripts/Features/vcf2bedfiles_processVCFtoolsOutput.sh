#!/bin/bash

# author: Emily Tsang
# updated for v8: Nicole Ferraro, summer 2018
echo "running vcf2bedfiles_processVCFtoolsOutput.sh"

set -o nounset
shopt -s expand_aliases

# helper script for vcf2bedfiles.sh
# processes output from vcf tools to create bed files for each individual
# the bed files contain MAF and genotype in 4th and 5th columns
# only includes sites where the individual carries the minor allele and the MAF <= 0.25
# (i.e. individual is either het (1) or homozygous for the minor allele (0 or 2, depending on the site))

# required files:
# GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_*.frq
# GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_*.GT.FORMAT

# process input parameter (which is the wildcard part of the filenames above)
if [ $# -ne 3 ]; then
    echo "usage: vcf2bedfiles_helper_processVCFtoolsOutput.sh <SNPs/indels/other> <vcftools prefix>"
    exit
fi



TYPE=$1
VCFPREFIX=$2
scriptdir=$3

##############
# FILE PATHS #
##############
dir=${VCFPREFIX%/*} # keep same directory as vcf prefix
outdir=${dir}/individualsRA # some future scripts (e.g. CADD one) assume this directory, so do not change this without changing downstream
mkdir -p $outdir
freq=${VCFPREFIX}_${TYPE}.frq
gt=${VCFPREFIX}_${TYPE}.GT.FORMAT
# output
#af=${dir}/AF_${TYPE}.bed
af=AF_${TYPE}.bed
# make sure output directory for individual beds exists
if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi


##################
# SOME FUNCTIONS #
##################
# skip first line of file
alias skip_header="tail -n+2"

# to remove missing genotypes
alias filter_genotypes="grep -v '\./\.' | grep -v '\.|\.'"

# to turn genotypes into 0,1,2 (the number of non-reference alleles)
# 0/0 -> 0; 0/*,*/0 -> 1; */* -> 2 (where * is any number >0)
# now considers phased genotypes (i.e. 0|0 -> 0, ...)
alias clean_genotypes="sed 's%0\(/\||\)0%0%' | sed 's%[1-9]\(/\||\)0%1%' | sed 's%0\(/\||\)[1-9]%1%' | sed 's%[1-9]\(/\||\)[1-9]%2%'"



# PROCESS ALLELE FREQUENCY FILE
# make bed file with minor allele frequency and whether the minor allele homozygote is 0 or 2
# for SNPs: columns 6 onwards are the alleles in order presented in the vcf with "." for columns up to 10 (ATCG*) if no extra alleles
# if there are multiple non reference alleles, the non-ref AF is set to 1-ref AF
# skips chrom X #RAU : ADD AN X PART
# Doesn't remove X and af=-nan for HallLab SV because there are duplicate entries for the same coordinate and that needs to be dealt with separately
# Changing to $2 and $2+1 to match up with gnomad (September 2018)
if [ ! -f ${dir}/aut_$af ]
then
	#autosomal
	if [ "$TYPE" = "SNPs" ]; then
    	cat $freq | tr ":" "\t" | skip_header | \
		awk 'BEGIN{OFS="\t"}{
               	af=$6;hz=0; 
               	if(af=="-nan" || $1=="chrX"){next}; 
               	if(af>0.5){af=1-$6;hz=2}; 
               	printf "%s\t%i\t%i\t%f\t%i",$1,$2,$2+1,af,hz; 
               	for(i=5;i<=NF-1;i=i+2){printf "\t%s",$i};
               	for(j=NF+1;j<=13;j=j+2){printf "\t."};
               	printf "\n"}' | \
		sort -k1,1 -k2,2n > ${dir}/aut_$af
	fi
	if [ "$TYPE" = "indels" ]; then
    	cat $freq | tr ":" "\t" | skip_header | \
		awk 'BEGIN{OFS="\t"}{
               	af=$6;hz=0; 
               	if(af=="-nan" || $1=="chrX"){next}; 
               	if(af>0.5){af=1-$6;hz=2};
               	print $1,$2,$2+1,af,hz}' | \
		sort -k1,1 -k2,2n > ${dir}/aut_$af
	fi

	# remove duplicate positions, if applicable.
	# this is to deal with the new gtex vcfs that have multi-allelic sites over multiple lines.
	# we are ignoring these sites.
	cut -f1-3 ${dir}/aut_$af | uniq -d > ${dir}/aut_${af}.dup
	bedtools subtract -f 1 -r -a ${dir}/aut_$af -b ${dir}/aut_${af}.dup > ${dir}/aut_${af}.dedupped
	mv ${dir}/aut_${af}.dedupped ${dir}/aut_$af
fi	

#x (repeat above)
if [ ! -f ${dir}/x_$af ]
then
	if [ "$TYPE" = "SNPs" ]; then
    	cat $freq | tr ":" "\t" | skip_header | \
		awk 'BEGIN{OFS="\t"}{
               	af=$6;hz=0; 
               	if(af=="-nan" || $1!="chrX"){next}; 
               	if(af>0.5){af=1-$6;hz=2}; 
               	printf "%s\t%i\t%i\t%f\t%i",$1,$2,$2+1,af,hz; 
               	for(i=5;i<=NF-1;i=i+2){printf "\t%s",$i};
               	for(j=NF+1;j<=13;j=j+2){printf "\t."};
               	printf "\n"}' | \
		sort -k1,1 -k2,2n > ${dir}/x_$af
	fi
	if [ "$TYPE" = "indels" ]; then
    	cat $freq | tr ":" "\t" | skip_header | \
		awk 'BEGIN{OFS="\t"}{
               	af=$6;hz=0; 
               	if(af=="-nan" || $1!="chrX"){next}; 
               	if(af>0.5){af=1-$6;hz=2};
               	print $1,$2,$2+1,af,hz}' | \
		sort -k1,1 -k2,2n > ${dir}/x_$af
	fi
	
	
	# remove duplicate positions, if applicable.
	# this is to deal with the new gtex vcfs that have multi-allelic sites over multiple lines.
	# we are ignoring these sites.
	cut -f1-3 ${dir}/x_$af | uniq -d > ${dir}/x_${af}.dup
	bedtools subtract -f 1 -r -a ${dir}/x_$af -b ${dir}/x_${af}.dup > ${dir}/x_${af}.dedupped
	mv ${dir}/x_${af}.dedupped ${dir}/x_$af
fi
				
# PROCESS GENOTYPE FILE
# do each column separately
# extract individual information from sample name
# get bed file of genotypes (in 0,1,2 format) and intersect with AF
# filter to keep heterozygotes and the minor allele homozygotes
# output bedfile with AF and genotype as columns 4 and 5

# deals with SVs entirely differently (takes care of AF filtering at this step)
# this is because there are duplicate entries for the same positions, so can't use a bed intersect
# probably could have dealt with SNPs and indels in a similar way, but keeping them as is because it already works.
ncol=`head -n1 $gt | wc -w`

# Get rid of printing "chr" beforehand (already included in v8)
processgt() {
    i=$1
    chr=$2
    sample=`head -n1 $gt | cut -f$i`
    echo "SAMPLE: "$sample", CHR: "$chr

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
}

cols=$(eval echo "{3..$ncol}")
echo "COLS: $cols"


#gt=${VCFPREFIX}_${TYPE}.GT.FORMAT
for i in $cols
do 
	sample=`head -n1 $gt | cut -f$i`
	echo "${outdir}/aut_${sample}_${TYPE}.bed"

	if [[ ! -f "${outdir}/aut_${sample}_${TYPE}.bed" && ! -f "${outdir}/aut_${sample}_${TYPE}.bed.gz" ]]
	then
		echo "submitting : ${outdir}/aut_${sample}_${TYPE}"
		sbatch ${scriptdir}/vcf2bedfiles_processVCFtoolsOutput_processgt.sh $i "aut" ${dir} ${af} ${outdir} ${TYPE} ${gt}
	fi
	if [[ ! -f "${outdir}/x_${sample}_${TYPE}.bed" && ! -f "${outdir}/x_${sample}_${TYPE}.bed.gz" ]]
	then
		echo "submitting : ${outdir}/x_${sample}_${TYPE}"
		sbatch ${scriptdir}/vcf2bedfiles_processVCFtoolsOutput_processgt.sh $i "x" ${dir} ${af} ${outdir} ${TYPE} ${gt}
	else
		echo "already have : ${outdir}/x_${sample}_${TYPE}"
	fi
done #| parallel --jobs 20 processgt {1}

