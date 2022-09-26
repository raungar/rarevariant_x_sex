#!/bin/bash


#infile=$1
#outfile=$2

infile="Output/sexdeg_v8/all_genes_sexDEGs_rvs_BREAST.txt.gz"

num_sexdegs=0
num_nonsexdegs=0
num_rvs_sexdegs=0
num_rvs_nonsexdegs=0
awk -F"\t" -v num_sexdegs=$num_sexdegs -v num_nonsexdegs=$num_nonsexdegs -v num_rvs_sexdegs=$num_rvs_sexdegs -v num_rvs_nonsexdegs=$num_rvs_nonsexdegs '{
	#if beta is not zero means that it is NOT a sexDEG, skip header
	if ($2 !="beta"){
		if ($2 == 0){
			num_nonsexdegs=num_nonsexdegs+1
			if ($4 == 1){
				num_rvs_nonsexdegs=num_rvs_nonsexdegs+1

			}
		} else {
			num_sexdegs=num_sexdegs+1
			if ($4 == 1){
				num_rvs_sexdegs=num_rvs_sexdegs+1
			}
		}

	}
}'

echo "FINALLY"

echo "$num_rvs_nonsexdegs over $num_nonsexdegs"
echo "--------------------------------------"
echo "$num_rvs_sexdegs over $num_sexdegs"

