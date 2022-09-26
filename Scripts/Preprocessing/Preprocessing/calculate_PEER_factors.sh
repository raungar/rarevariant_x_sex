#!/bin/bash
#SBATCH --job-name=peer_factors
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
##SBATCH --time=06:59:00
#SBATCH --time=03:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_factors_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_factors_%j.err"

set -o nounset -o errexit -o pipefail

#####FOR SCG
module load legacy/.base R/3.3.1
srun="srun -n1 -N1 --exclusive"

peerdir=$1 #${RAREDIR}/preprocessing_v8/PEER_v8
scriptdir=$2 #`dirname \$(readlink -f "\$0")`
RAREDIR=$3
traitsFileName=$4
## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

echo "beginning"
#for traitsFileName in `ls ${peerdir}/*log2.ztrans*.txt`
#do

    prefix=${traitsFileName%.log2.ztrans*txt}
    #prefix=${traitsFileName%.reads.txt}
    #sex=`echo $traitsFileName | awk -F"\\." '{print $(NF-1)}'  | awk -F"\\." '{print $1}'`
    sex=`echo $traitsFileName | awk -F"\\." '{print $(NF-1)}'`
    #prefix=${traitsFileName%.tpm.log2.ztrans.txt}
    tissue=`basename "$prefix"`
    echo "sex is $sex and tissue is $tissue"
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

    ## actual calculation of peer factors
    echo "Rscript calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" > ${outdir}/log.txt
     Rscript ${scriptdir}/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue >> ${outdir}/log.txt 2>&1
#  break
#done
