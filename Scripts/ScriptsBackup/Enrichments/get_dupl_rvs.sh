#!/bin/bash



rvs_all=$1 #Output/enrichments_v8/all_rvs_inds_types.txt.gz
temp_out=$2 #Output/enrichments_v8/temp_out.txt
out=$3 #Output/enrichments_v8/all_rvs_inds_types_repeated.txt

if [ ! -f $temp_out ]
then
	"Creating temporary file of duplicated RVs"
	less $rvs_all | awk '{print $2"\t"$3}' | sort | uniq -c | awk '$1 > 1 { print $2"\t"$3}' > $temp_out
fi

while read line
do
	zgrep "$line" $rvs_all >> $out
done<$temp_out

gzip $out
