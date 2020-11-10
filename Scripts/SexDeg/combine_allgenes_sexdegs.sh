#!/bin/bash

all_genes_file=$1
sexdegs_file=$2
outfile=$3

#all_genes_file="Output/sexdeg_v8/all_genes.txt"
#sexdegs_file="Files/SexDEGs/sexDEGS-BREAST.csv"


#read in all genes into a dictionary, default beta (sexdeg) is 0
declare -A genes_dic
while read line;
do
	genes_dic["${line}"]=0

done < $all_genes_file


#now actually read the sexdeg file. for each gene assign the sexdeg value in the dictionary
i=1
while read l
do
	if [[ $i == 1 ]]
	then
		i=$(($i+1))
		continue
	fi

	#echo $l
	ensg=`echo $l | awk -F"," '{print $1}'` #get teh gene name
	beta=`echo $l | awk -F"," '{print $4}'` #get the beta (effect size of sexDEG)
	genes_dic["${ensg}"]=$beta
	#echo ${genes_dic["${ensg}"]}

done < $sexdegs_file


#finally, print the dictionary such that
#key is ensg, val is beta value
#print so tab delimeted ensg/beta
for key in "${!genes_dic[@]}"
do 
	echo -e $key"\t"${genes_dic["${key}"]} >> "${outfile}"
done

