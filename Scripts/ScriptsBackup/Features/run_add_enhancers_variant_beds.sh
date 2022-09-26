#!/bin/bash

# Author: Emily Tsang
# Updated: Nicole Ferraro, summer 2018
# Adapted for SCG: Rachel Ungar, Spring 2020

set -o nounset -o errexit -o pipefail

echo "running run_add_enhancers_variant_beds.sh"

RAREDIR=${1}
scriptdir=${2}
GNOMAD_SNPS=${3}
GNOMAD_INDELS=${4}
AF_SVs=${5}

# get current working directory of the script


## set output directory
dir=${RAREDIR}/features_v8

outdir=${dir}/bySite
AF_SVS=${AF_SVs}
mkdir -p $outdir

echo "header"
#																																								
sbatch_header="#!/bin/bash\n#SBATCH --job-name=add_features\n#SBATCH --cpus-per-task=2\n#SBATCH --partition=interactive\n#SBATCH --account=default\n#SBATCH --time=03:30:00\n#SBATCH --mem-per-cpu=2G\n#SBATCH --output='/oak/stanford/groups/smontgom/raungar/Sex/Jobs/feature_sexgtex_%j.out'\n#SBATCH --error='/oak/stanford/groups/smontgom/raungar/Sex/Jobs/feautre_sexgtex_%j.err'\n"

#sbatch ${scriptdir}/submit_add_variant_beds1.sh

echo "SBATCH1"
echo -e $sbatch_header > ${scriptdir}/submit_add_enhac_variant_beds1.sh
echo "sh ${scriptdir}/add_predEnhancers_variant_beds.sh ${dir}/variantBeds/individuals/*_SNPs.bed.gz ${GNOMAD_SNPS} ${outdir}" >> ${scriptdir}/submit_add_enhac_variant_beds1.sh
sbatch ${scriptdir}/submit_add_enhac_variant_beds1.sh
#parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_SNPs.bed.gz ::: ${GNOMAD_SNPS} ::: ${outdir}


echo "SBATCH2"
echo -e $sbatch_header > ${scriptdir}/submit_add_enhac_variant_beds2.sh
echo "sh ${scriptdir}/add_predEnhancers_variant_beds.sh ${dir}/variantBeds/individuals/*_indels.bed.gz ${GNOMAD_INDELS} ${outdir}" >> ${scriptdir}/submit_add_enhac_variant_beds2.sh
sbatch ${scriptdir}/submit_add_enhac_variant_beds2.sh
#parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_indels.bed.gz ::: ${GNOMAD_INDELS} ::: ${outdir}

outdir=${dir}/bySite/HallLabSV
echo "SBATCH3"
echo -e $sbatch_header > ${scriptdir}/submit_add_predEnhancers_variant_beds3.sh
echo "sh ${scriptdir}/add_predEnhancers_variant_beds.sh ${dir}/variantBeds/individuals/HallLabSV_hg38/*_HallLabSV.bed.gz ${AF_SVS} ${outdir}" >> ${scriptdir}/submit_add_predEnhancers_variant_beds3.sh
sbatch ${scriptdir}/submit_add_predEnhancers_variant_beds3.sh
#parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/HallLabSV_hg38/*_HallLabSV.bed.gz ::: ${AF_SVS} ::: ${outdir}



wait
