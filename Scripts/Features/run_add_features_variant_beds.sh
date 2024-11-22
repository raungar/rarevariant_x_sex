#!/bin/bash

# Author: Emily Tsang
# Updated: Nicole Ferraro, summer 2018
# Adapted for SCG: Rachel Ungar, Spring 2020

set -o nounset -o errexit -o pipefail

echo " RUNNING run_add_features_variant_beds.sh"

RAREDIR=${1}
scriptdir=${2}
KG_AF_SNPS=${3}
GNOMAD_SNPS=${4}
KG_AF_INDELS=${5}
GNOMAD_INDELS=${6}
AF_SVs=${7}
CONSOLIDATED=${8}
TFDIR=${9}
STATEDIR=${10}

echo "11: STATEDIR"
echo ${STATEDIR}
# get current working directory of the script
#scriptdir=`dirname \$(readlink -f "\$0")`

## run add_features_variant_beds.sh in batches

## arguments: takes as input the version and (optionally) the number of cores to use
#if [ $# -eq 0 ]; then
#    ncores=25
#elif [ $# -eq 1 ]; then
#    ncores=$1
#else
#    echo "usage: run_add_features_variant_beds.sh [ncores]"
#    exit
#fi

## set output directory
dir=${RAREDIR}/features_v8

outdir=${dir}/bySite

mkdir -p $outdir

#parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individualsRA/*_SNPs.bed.gz ::: ${KG_AF_SNPS} ::: ${GNOMAD_SNPS} ::: ${outdir}
#parallel --jobs $ncores --xir}/variantBeds/individualsRA/*_SNPs.bed.gz ::: ${KG_AF_SNPS} ::: ${GNOMAD_SNPS} ::: ${outdir}

sbatch_header="#!/bin/bash\n#SBATCH --job-name=add_features\n#SBATCH --cpus-per-task=4\n#SBATCH --partition=interactive\n#SBATCH --account=default\n#SBATCH --time=08:59:00\n#SBATCH --mem-per-cpu=6G\n#SBATCH --output=\"/oak/stanford/groups/smontgom/raungar/Sex/Jobs/features_sexgtex_%j.out\"\n#SBATCH --error=\"/oak/stanford/groups/smontgom/raungar/Sex/Jobs/features_sexgtex_%j.err\"\n"

if [[ -f ${dir}/variantBeds/individualsRA/*bed ]]
then
	for file in `ls ${dir}/variantBeds/individualsRA/*bed`
	do
		gzip "$file"
	done
fi

counter=0
for this_file in `ls ${dir}/variantBeds/individualsRA/x_*indels.bed.gz`
do
counter=$counter+1
echo date
echo $this_file 
file=`echo $this_file | awk -F"_indels.bed.gz" '{print $1}'`
ID=`echo $file | awk -F"/" '{print $NF}'`
echo "SBATCH1"
echo -e $sbatch_header > ${scriptdir}/submit_add_variant_beds1.sh

#echo "${scriptdir}/add_features_variant_beds.sh  ${dir}/variantBeds/individualsRA/${file}_SNPs.bed.gz \\" >> ${scriptdir}/submit_add_variant_beds1.sh
echo "${scriptdir}/add_features_variant_beds.sh  ${file}_SNPs.bed.gz \\" >> ${scriptdir}/submit_add_variant_beds1.sh
echo "${KG_AF_SNPS} ${GNOMAD_SNPS} ${outdir} \\" >> ${scriptdir}/submit_add_variant_beds1.sh
echo "${STATEDIR} ${TFDIR} \\" >> ${scriptdir}/submit_add_variant_beds1.sh
echo "${CONSOLIDATED}" >> ${scriptdir}/submit_add_variant_beds1.sh
if [ ! -f  ${outdir}/${ID}_SNPs_features.bed.gz ]
then
	echo "submitting ${ID} snps"
	###sbatch ${scriptdir}/submit_add_variant_beds1.sh
fi
#parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/${file}_indels.bed.gz ::: ${KG_AF_INDELS} ::: ${GNOMAD_INDELS} ::: ${outdir}
#parallel --jobs $ncores --xapply ::: ${KG_AF_INDELS} ::: ${GNOMAD_INDELS} ::: ${outdir}
echo "SBATCH2"
echo -e	$sbatch_header > ${scriptdir}/submit_add_variant_beds2.sh
echo "${scriptdir}/add_features_variant_beds.sh ${file}_indels.bed.gz \\" >> ${scriptdir}/submit_add_variant_beds2.sh
echo "${KG_AF_INDELS} ${GNOMAD_INDELS} ${outdir} \\" >> ${scriptdir}/submit_add_variant_beds2.sh
echo "${STATEDIR} ${TFDIR} \\" >> ${scriptdir}/submit_add_variant_beds2.sh
echo "${CONSOLIDATED}" >> ${scriptdir}/submit_add_variant_beds2.sh


if [ ! -f  ${outdir}/${ID}_indels_features.bed.gz ]
then
	echo "submitting ${ID} indels" 
	###sbatch ${scriptdir}/submit_add_variant_beds2.sh
fi
#parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/${file}_HallLabSV.bed.gz ::: ${AF_SVs} ::: ${outdir}
#parallel --jobs::: ${AF_SVs} ::: ${outdir}
#echo "SBATCH3"
#echo -e	$sbatch_header > ${scriptdir}/submit_add_variant_beds3.sh
#echo  "${scriptdir}/add_features_variant_beds.sh ${dir}/variantBeds/individuals/${file}_HallLabSV.bed.gz \\"  >> ${scriptdir}/submit_add_variant_beds3.sh
#echo  "${scriptdir}/add_features_variant_beds.sh ${file}_HallLabSV.bed.gz \\"  >> ${scriptdir}/submit_add_variant_beds3.sh
#echo "${AF_SVs} ${outdir} ${STATEDIR} \\" >> ${scriptdir}/submit_add_variant_beds3.sh
#echo "${TFDIR} \\" >> ${scriptdir}/submit_add_variant_beds3.sh
#echo "${STATEDIR} ${CONSOLIDATED}" >> ${scriptdir}/submit_add_variant_beds3.sh
#sbatch ${scriptdir}/submit_add_variant_beds3.sh

if [[ $counter == 20 ]] 
then
	break
fi

done

wait
