#!/bin/bash

#												
dir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8"
file_b="outliers_zthresh3_nphen5_noglobal_medz_varAnnot_aut_both_half.regress.txt "
file_m="outliers_zthresh3_nphen5_noglobal_medz_varAnnot_aut_m.txt"
file_f="outliers_zthresh3_nphen5_noglobal_medz_varAnnot_aut_f.txt"

if [ ! -f Output/Analysis/outlier_aut.txt ]
then
	grep outlier $dir/$file_b | awk -F"\t" '{print $1"\t"$2"\t"$5"\tboth"}' > Output/Analysis/outlier_aut.txt
	grep outlier $dir/$file_m | awk -F"\t" '{print $1"\t"$2"\t"$5"\tfemale"}' >> Output/Analysis/outlier_aut.txt
	grep outlier $dir/$file_f | awk -F"\t" '{print $1"\t"$2"\t"$5"\tmale"}' >> Output/Analysis/outlier_aut.txt
fi

awk -F"\t" '{print $1"\t"$2}' Output/Analysis/outlier_aut.txt | sort | uniq -c | awk '{if($1==1){print $2"\t"$3}}' > Output/Analysis/oulier_aut_only1.txt


while IFS= read -r line
do
	#mygrep=`grep "$line" $dir/$file_b`
	#if [ -z "$mygrep" ]
	if grep -q "$line" $dir/$file_b
	then
		grep "$line" $dir/$file_b >> Output/Analysis/oulier_aut_only1_grepped.txt
	fi
	
done < Output/Analysis/oulier_aut_only1.txt
