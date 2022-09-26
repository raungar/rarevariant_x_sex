set -o nounset -o errexit -o pipefail


#####FOR SCG
module load R
srun="srun -n1 -N1 --exclusive"



## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

peerdir=$1 #${RAREDIR}/preprocessing_v8/PEER_v8
scriptdir=$2 #`dirname \$(readlink -f "\$0")`
gtex_v8_eqtl_dir=$3 #${GTEXv8}/eqtl/GTEx_Analysis_v8_eQTL
RAREDIR=$4
pcs=$5
md=$6
#parallel="parallel -N 1 --delay .2 -j 10 --joblog parallel_joblog --resume"
#parallel --jobs 10 runPeer ::: ${peerdir}/*.log2.ztrans*.txt
#$srun sh $scriptdir/calculate_PEER.sh $1 $2 $3 $4 ::: ${peerdir}/Adipose_Subcutaneous.log2.ztrans*.txt
#$srun sh $scriptdir/calculate_PEER.sh $1 $2 $3 $4 ::: ${peerdir}/Adipose_Subcutaneous.log2.ztrans*.txt
#$srun ./oak/stanford/groups/smontgom/raungar/Sex/Scripts/Preprocessing/runPeer.sh ::: ${peerdir}/Adipose_Subcutaneous.log2.ztrans*.txt

for file in `ls ${peerdir}/*log2.ztrans*.txt`
do
	sbatch $scriptdir/calculate_PEER.sh $1 $2 $3 $4 $file $5 $6
	echo "$scriptdir/calculate_PEER.sh $1 $2 $3 $4 $file $5 $6"
	#break
done


#echo "DONE!"
echo "submitted"
