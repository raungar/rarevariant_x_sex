#!/bin/bash
#SBATCH --job-name=peer
#SBATCH --cpus-per-task=4
#SBATCH --partition=interactive
#SBATCH --account=default
##SBATCH --time=06:59:00
#SBATCH --time=02:59:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/peer_%j.err"

set -o nounset -o errexit -o pipefail
echo $1
echo $2
echo $3
echo $4
echo $5
echo $6 #md
#####FOR SCG
#module load R peer

module load legacy/.base R/3.3.1
srun="srun -n1 -N1 --exclusive"



## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples



    traitsFileName=$5


    peerdir=$1 #${RAREDIR}/preprocessing_v8/PEER_v8
    scriptdir=$2 #`dirname \$(readlink -f "\$0")`
    gtex_v8_eqtl_dir=$3 #${GTEXv8}/eqtl/GTEx_Analysis_v8_eQTL
    RAREDIR=$4
    pcs=$6
    md=$7
    factors_type=$8




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


    which R
    ## actual calculation of peer factors
    ###echo "Rscript calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" > ${outdir}/log.txt
    ##Rscript ${scriptdir}/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
    ##        $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue >> ${outdir}/log.txt 2>&1
    
    # computing residuals
    echo "${scriptdir}/calculate_PEER_residuals.R $traitsFileName ${pcs} \
            ${indir}/factors.tsv ${gtex_v8_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
        $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
        ${prefix}.${sex}.peer.v8ciseQTLs.ztrans.txt &> ${outdir}/log.residuals.txt  "

    Rscript ${scriptdir}/calculate_PEER_residuals.R $traitsFileName ${pcs} \
            ${indir} ${gtex_v8_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
        $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
        ${prefix}.${sex}.peer.v8ciseQTLs.ztrans.txt ${md} ${factors_type} &> ${outdir}/log.residuals.txt  

