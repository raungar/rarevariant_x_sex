#!/bin/bash


infile=$1
outfile_aut=$2
outfile_x=$3

#infile="Output/sexdeg_v8/all_genes_sexDEGs_rvs_BREAST.txt.gz"
#infile="Output/sexdeg_v8/all_genes_sexDEGs_rvs_BREAST.txt.gz"

num_sexdegs=0
num_nonsexdegs=0
num_rvs_sexdegs=0
num_rvs_nonsexdegs=0
vals=`less $infile | awk -F"\t" -v num_sexdegs=$num_sexdegs -v num_nonsexdegs=$num_nonsexdegs -v num_rvs_sexdegs=$num_rvs_sexdegs -v num_rvs_nonsexdegs=$num_rvs_nonsexdegs '{
	#if beta is not zero means that it is NOT a sexDEG, skip header
	if ($1 == "chrX"){
		if ($3 !="beta"){
			if ($3 == 0){
				num_nonsexdegs=num_nonsexdegs+1
				if ($5 != 1){
					num_rvs_nonsexdegs=num_rvs_nonsexdegs+1

				}
			} else {
				num_sexdegs=num_sexdegs+1
				if ($5 != 1){
					num_rvs_sexdegs=num_rvs_sexdegs+1
				}
			}
		print num_rvs_nonsexdegs,num_nonsexdegs,num_rvs_sexdegs,num_sexdegs
		}
	}
}' | tail -1`

echo "num_rvs_nonsexdegs num_nonsexdegs num_rvs_sexdegs num_sexdegs" > $outfile_x
echo "$vals" >> $outfile_x
#echo "$num_rvs_sexdegs over $num_sexdegs"
#echo "--------------------------------------"
#echo "$num_rvs_nonsexdegs over $num_nonsexdegs"


num_sexdegs=0
num_nonsexdegs=0
num_rvs_sexdegs=0
num_rvs_nonsexdegs=0
vals=`less $infile | awk -F"\t" -v num_sexdegs=$num_sexdegs -v num_nonsexdegs=$num_nonsexdegs -v num_rvs_sexdegs=$num_rvs_sexdegs -v num_rvs_nonsexdegs=$num_rvs_nonsexdegs '{
	#if beta is not zero means that it is NOT a sexDEG, skip header
	if ($1 != "chrX"){
		if ($3 !="beta"){
			if ($3 == 0){
				num_nonsexdegs=num_nonsexdegs+1
				if ($5 != 1){
					num_rvs_nonsexdegs=num_rvs_nonsexdegs+1

				}
			} else {
				num_sexdegs=num_sexdegs+1
				if ($5 != 1){
					num_rvs_sexdegs=num_rvs_sexdegs+1
				}
			}
		print num_rvs_nonsexdegs,num_nonsexdegs,num_rvs_sexdegs,num_sexdegs
		}
	}
}' | tail -1`

echo "num_rvs_nonsexdegs num_nonsexdegs num_rvs_sexdegs num_sexdegs" > $outfile_aut
echo "$vals" >> $outfile_aut
#echo "$num_rvs_sexdegs over $num_sexdegs"
#echo "--------------------------------------"
#echo "$num_rvs_nonsexdegs over $num_nonsexdegs"

