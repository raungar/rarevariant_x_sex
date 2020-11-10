#!/bin/bash

rv_file=$1
sexdegs_file=$2
outfile=$3

#rv_file="Output/features_v8/collapsed_maf_both_x.tsv.gz"
#sexdegs_file="Output/sexdeg_v8/all_genes_sexDEGs_BREAST.txt"

#read in all genes into a dictionary, default beta (sexdeg) is 0
declare -A genes_dic
while read line;
do
	ensg=`echo $line | awk '{print $1}' `
	beta=`echo $line | awk '{print $2}' `
	genes_dic[$ensg]="$beta,0,1,NA" # this will be beta, # inds, min MAF, sex [na/m/f/b]
	#genes_dic[$ensg]=["$beta","NA",1,"NA"]
	#echo "$beta,NA,1,NA"
done < $sexdegs_file


i=1
gzcat $rv_file | while read l
do
	if [[ $i == 1 ]]
	then
		i=$(($i+1))
		continue
	fi

	#echo $l
	ensg=`echo $l | awk -F '{print $3}'` #get teh gene name
	#ind= `echo $l | awk -F '{print $5}'` #ind seen
	sex=`echo $l | awk -F '{print $6}'` #either male, female, or both
	maf=`echo $l | awk -F '{print $10}'` 

	gene_match=genes_dic["${ensg}"]
	gene_match_beta=`awk -F"," '{print $1}'`
	gene_match_ind=`awk -F"," '{print $2}'`
	gene_match_mafmin=`awk -F"," '{print $3}'`
	gene_match_sex=`awk -F"," '{print $4}'`
	#if gene_match_sex is NA , then no rv has been found near it yet, so skip ahead
	if [[ $gene_match_sex != "NA" ]]
	then
		#record this smallest maf
		if [[ $gene_match_mafmin < $maf ]]
		then
			maf=$gene_match_mafmin
		fi

		if [[ $gene_match_sex != $sex ]]
		then
			sex="both"
		fi
	fi
	
	ind=$(($gene_match_ind+1))
	genes_dic["${ensg}"]="$gene_match_beta,$ind,$sex,$maf"

done #< $rv_file

#finally, print the dictionary such that
#key is ensg, val is beta value
#print so tab delimeted ensg/beta
for key in "${!genes_dic[@]}"
do 
	echo $key","${genes_dic["${key}"]} >> "${outfile}"
done

