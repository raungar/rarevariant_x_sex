#!/bin/bash

rv_file=$1
sexdegs_file=$2
outfile=$3

#rv_file="Output/features_v8/collapsed_maf_both_x.tsv.gz"
#sexdegs_file="Output/sexdeg_v8/all_genes_sexDEGs_BREAST.txt"

#read in all genes into a dictionary, default beta (sexdeg) is 0
declare -A genes_dic

j=1
while read line;
do
	j=$(($j+1))
	ensg=`echo $line | awk '{print $1}' `
	beta=`echo $line | awk '{print $2}' `
	genes_dic[$ensg]="$beta,0,1,NA" # this will be beta, # inds, min MAF, sex [na/m/f/b]
	#genes_dic[$ensg]=["$beta","NA",1,"NA"]
	#echo "$beta,NA,1,NA"


done < $sexdegs_file




#echo "${!genes_dic[@]}"
#echo "above is ALLLLLL"


i=1
zcat $rv_file | while read l
do
	if [[ $i == 1 ]]
	then
		i=$(($i+1))
		continue
	fi

	#echo $l
	ensg=`echo $l | awk '{print $3}'` #get teh gene name
	#ind= `echo $l | awk'{print $5}'` #ind seen
	sex=`echo $l | awk '{print $6}'` #either male, female, or both
	maf=`echo $l | awk '{print $10}'` 
	########ensg="ENSG00000280233.1"
	##### FIX HERE ---- OUTPUT SEEMS TO BE REVERTING TO LINE FOR SOME REASON
	##### PAST THIS POINT VARIABLES SEEM TO CHANGE.....
	gene_match=`echo ${genes_dic[$ensg]}`
	gene_match_beta=`echo $gene_match | awk -F"," '{print $1}'`
	gene_match_ind=`echo $gene_match | awk -F"," '{print $2}'`
	gene_match_mafmin=`echo $gene_match | awk -F"," '{print $3}'`
	gene_match_sex=`echo $gene_match | awk -F"," '{print $4}'`
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
	echo "$gene_match_beta,$ind,$sex,$maf"
	genes_dic[$ensg]="$gene_match_beta,$ind,$sex,$maf"
done #< $rv_file

#finally, print the dictionary such that
#key is ensg, val is beta value
#print so tab delimeted ensg/beta
for key in "${!genes_dic[@]}"
do 
	echo $key","${genes_dic["${key}"]} >> "${outfile}"
done

