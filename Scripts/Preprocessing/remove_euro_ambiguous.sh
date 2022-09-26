#!/bin/bash


euro_ids="Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
euro_amb="Output/preprocessing_v8/euro_ambiguous.txt"
outfile="Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids_notambiguous.txt"

while read l
do
	in_amb=`grep $l $euro_amb`
	if [[ $in_amb == "" ]]
	then
		echo $l >> $outfile
	fi

done<$euro_ids
