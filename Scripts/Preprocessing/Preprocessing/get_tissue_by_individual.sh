#!/bin/bash


dir=$1
# Go though GTEx data to figure out which individuals have data for which tissues.
# Create a file that has 2 columns:
# one column with the individual id and one column with the tissue name.
# Also generate files with the lists of unique tissue names and individual IDs.

set -o nounset -o errexit -o pipefail

#dir=${RAREDIR}/preprocessing_v8

#group=("m" "f" "both" "both.sex" "both.regress" "both_half" "both_half.sex" "both_half.regress" )
group=("m" "f" "both_half" "both")
#group=("m" "f" "both_half" "both_half.sex_regress" )

for group in ${group[@]}
do
	out=${dir}/gtex_2017-06-05_tissue_by_ind_${group}.txt 
	tissues=${dir}/gtex_2017-06-05_tissues_all_normalized_samples_${group}.txt
	inds=${dir}/gtex_2017-06-05_individuals_all_normalized_samples_${group}.txt

	echo ${out}
	echo ${tissues}
	echo ${inds}

	# put header
	echo -e "Tissue\tId" > $out


	for f in ${dir}/PEER_v8/*.${group}.peer.ztrans.txt
	do
	    fname=`basename $f`
	    tissue=${fname%.${group}.peer.ztrans.txt}
	    head -n1 $f | awk -v tissue=$tissue '{for(i=2; i<=NF; i++){print tissue"\t"$i}}' >> $out
	done

	# get the lists of individuals and tissues separately
	cut -f1 $out | sort | uniq | grep -P -v '^Tissue' > $tissues
	cut -f2 $out | sort | uniq | grep -P -v '^Id' > $inds
done
