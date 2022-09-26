#!/bin/bash
#SBATCH --job-name=peer
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
##SBATCH --time=06:59:00
#SBATCH --time=02:59:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_%j.err"

set -o nounset -o errexit -o pipefail

module load legacy/.base R/3.3.1
srun="srun -n1 -N1 --exclusive"



## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

peerdir=${1} #${RAREDIR}/preprocessing_v8/PEER_v8
scriptdir=${2} #`dirname \$(readlink -f "\$0")`
RAREDIR=${3}
pcs=${4}
md=${5}
factors_type=${6} #factors or factors prefix for when transform by sex
logfile=${7} #log file change
incl_sex=${8}  # T or F
sex_contin=${9}
traitsFileName=${10}

echo "${traitsFileName}"
#for traitsFileName in `ls ${peerdir}/*log2.ztrans*.txt`
#do
    prefix=${traitsFileName%.log2.ztrans*txt}
    sex=`echo $traitsFileName | awk -F"\\." '{print $(NF-1)}'`
    #prefix=${traitsFileName%.tpm.log2.ztrans.txt}
    tissue=`basename "$prefix"`
    nsamples=$(cat $traitsFileName | wc -l) # this is actually n samples + 1
    if [ $nsamples -le 150 ]; then
        maxFactorsN=15
    elif [ $nsamples -le 249 ]; then
        maxFactorsN=30
    elif [ $nsamples -le 349 ]; then
        maxFactorsN=45
    else
        maxFactorsN=60
    fi
    maxIterations=10000
    boundTol=0.001
    varTol=0.00001
    e_pa=0.1
    e_pb=10
    a_pa=0.001
    a_pb=0.1
    outdir=${prefix}_Factors"$maxFactorsN"_$sex
    indir=${prefix}_Factors"$maxFactorsN"_$sex
    echo $outdir
    mkdir -p $outdir

    if [[ "$incl_sex" == "T" ]]
    then
        sex="${sex}.sex"
    fi
    
    if [[ "$factors_type" == *"regress"* ]]
    then
        sex="${sex}_regress"
    fi

    outfile=${prefix}.${sex}.peer.ztrans.txt

   echo $outfile

    # computing residuals
    echo "${scriptdir}/calculate_PEER_residuals.R $traitsFileName ${pcs} \
            ${indir}/factors.tsv  \
        ${prefix}.${sex}.peer.ztrans.txt &> ${outdir}/log.residuals.txt  "

	#3 and 4 removed
    Rscript ${scriptdir}/calculate_PEER_residuals.R $traitsFileName ${pcs} \
            ${indir} \
        ${outfile} ${md} ${factors_type} ${incl_sex} ${sex_contin} &> ${outdir}/${logfile}

    #break
#done
