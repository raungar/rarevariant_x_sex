# author: Emily Tsang
#rachel replaced intersectBed w bedtool intersect to preoperly deal with wildcards

# script to put together some features for the GTEx WGS samples
# restricts to variants with an AF of 0.25 or less

# INPUT
# * input bed file, see format below
# * 1000 genomes allele frequency file
#   (bed file with 5 additional columns with the AF for each super population [EAS, AMR, AFR, EUR, SAS])
# * output directory

# REQUIRES
# * $CONSOLIDATED (15 cols, with chr)
# * $KG_AF_INDELS (8 cols, with chr)
# * $KG_AF_SNPS (ditto)
# * $TFDIR/*merged.bed.gz (3 cols, with chr)
# * $STATEDIR/[prom/dyadic/enh]/*.bed.gz (3 cols, with chr)

set -o nounset -o errexit -o pipefail


date

######################################

# check that there are the correct arguments
#if [ $# -ne 4 ]; then
#    echo "usage: add_features_variant_beds.sh input_file_name 1KG_AF_bed_file outdir"
#    exit
#fi
echo "add_features_variant_beds.sh"

f=${1}
kgfile=${2}
gnomadfile=${3}
outdir=${4}
STATEDIR=${5}
TFDIR=${6}
CONSOLIDATED=${7}

fname=`basename ${f}`
ind=${fname%%.*}


if [ -f ${outdir}/${ind}_features.bed.gz ]
then
	echo "already exists"
	exit "This file already exists, no need to rerun!"
fi


echo ${ind}
echo "filename: ${f}"
echo "STATEDIR: ${STATEDIR}"
echo "TFDIR: ${TFDIR}"
echo "CONSOLIDATED: ${CONSOLIDATED}"
#################################
# input bed file has
# 1* chromosome (with "chr")
# 2* Pos0
# 3* Pos1
# 4* MAF
# 5* genotype (0,1,2)
# 6* CADD raw score (SNP files only; add NA for indels)
# 7* CADD phred score (SNP files only; add NA for indels)
#
# add the following annotations
# 8*  number of TFBS overlaps (Pouya's Roadmap motifs)
# 9*  number of cell lines with promoters (from Wouter Meuleman's call for Roadmap Epigenomics)
# 10* number of cell lines with enhancers (ditto)
# 11* number of cell lines with dyadic elements (ditto)
# 12* CpG (this and all features through 26 are from CADD)
# 13* priPhCons
# 14* mamPhCons
# 15* verPhCons
# 16* priPhyloP
# 17* mamPhyloP
# 18* verPhyloP
# 19* GerpN
# 20* GerpS
# 21* GerpRS
# 22* GerpRSpval
# 23* fitCons
# 24* TFBS
# 25* TFBSPeaks
# 26* TFBSPeaksMax
# 27* 1KG allele frequency from 5 super populations (5 columns: EAS, AMR, AFR, EUR, SAS)
##################################

module load bedtools/2.27.1 

# If you want to add extra columns, that you just need a count (-c),
# put them before the consolidated/sorted.consolidated.annotations.bed.gz
# then adjust the "allend" variable in the last awk line
# otherwise, just add them at the end and adjust the awk accordingly
# also make sure to use the -loj option

# NOTE: using -loj output is slightly messed up, so need to use sed to remove two consecutive tabs (see below)

# first add two NA columns for indels that have only 5 columns (no CADD scores)
zcat $f | sort -k 1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{if(NF==5){print $0,"NA","NA"} else {print $0}}' | \
    bedtools intersect -sorted -wa -c \
    -a stdin \
    -b ${TFDIR}/*.merged.bed.gz | \
    bedtools intersect -sorted -wa -c \
    -a stdin \
    -b ${STATEDIR}/prom/*.bed.gz | \
    bedtools intersect -sorted -wa -c \
    -a stdin \
    -b ${STATEDIR}/enh/*.bed.gz | \
    bedtools intersect -sorted -wa -c \
    -a stdin \
    -b ${STATEDIR}/dyadic/*.bed.gz | \
    bedtools intersect -sorted -wa -wb -loj \
    -a stdin \
    -b ${CONSOLIDATED} | \
    sed 's/\t\t/\t/g' | \
    bedtools intersect -sorted -wa -wb -loj \
    -a stdin \
    -b ${kgfile} | \
    bedtools intersect -sorted -wa -wb -loj \
    -a stdin \
    -b ${gnomadfile} | \
    sed 's/\t\t/\t/g' | \
    awk 'BEGIN{roadmapEnd=11; consEnd=29}{
         for(i=1; i<=roadmapEnd; i++){ printf("%s\t",$i)}; 
         for(i=roadmapEnd+4; i<=consEnd; i++){ printf("%s\t",$i)}; 
         for(i=consEnd+4; i<NF; i++){ printf("%s\t",$i)}; 
         print $NF}' | \
    sed 's/\t\./\tNA/g' | \
    gzip -c \
    > ${outdir}/${ind}_features.bed.gz

echo "${fname} done"

date
