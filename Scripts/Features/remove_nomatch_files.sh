#!/bin/bash

i=0
for file in `ls Output/features_v8/bySiteXFIX/x_GTEX-X*`
do
	file_len=`zcat $file |  awk -F"\t" '{print $35}' | sort | uniq -c | wc -l`
	echo "$file has $file_len columns"
	if [[ $file_len -eq 1 ]]
	then
		echo "need to remove this ain't good"
		rm $file
		echo "remoevd $file"
		echo $file >> removed_files.txt


	fi

	i=$(($i+1))
done
