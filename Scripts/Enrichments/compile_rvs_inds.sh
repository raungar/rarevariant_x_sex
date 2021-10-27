#!/bin/bash


file_dir="Output/features_v8/bySite"
pattern="x*gz"
outfile=""

echo -e chr"\t"start"\t"end"\t"maf"\t"sample"\t"vartype > ${outfile}

for f in $file_dir/$pattern
do
	echo $f
	gtex_id=`echo $f | awk -F"/" '{print $NF}' | grep -o -P 'GTEX.{0,6}'`
	var_type=`echo $f | awk -F"/" '{print $NF}' | egrep -o 'indel|SNP|SV'`
	zcat $f | awk -F"\t" -v gtex_id=$gtex_id -v var_type=$var_type '{print $1"\t"$2"\t"$3"\t"$4"\t"gtex_id"\t"var_type}' \
	>> ${outfile}
done

