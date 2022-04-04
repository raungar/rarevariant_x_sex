#!/bin/bash

gtf=$1
out_prefix=$2
filt_name=$3

echo "gtf input is: $gtf"
echo "out_prefix input is: $out_prefix"
echo "only grabbing: $filt_name"


#check if files exist
if ls "$out_prefix*gtf" 1> /dev/null 2>&1
then
	rm "$out_prefix*gtf"
fi
echo "any existing gtf files removed."


while read line
do
	echo "reading: $line"
	type=`echo $line | awk '{print $3}'`	
	echo $type
	if [[ $type != "$filt_name" ]]
	then
		continue
	fi
	chr_num=`echo $line | awk '{print $1}'`
	echo "$out_prefix$chr_num.presort.gtf"
	echo $line | sed 's/\s/\t/g' | awk -F"\t" '{for(i=1; i<=20; i++){printf $i"\t"} ;print""}'  >> $out_prefix$chr_num.presort.gtf

done < $gtf

####this is necessary to sort after alll is split
echo "sort gtf files"
for file in `ls $out_prefix*presort.gtf`
do
	file_chrnum=`echo $file | awk -F"/" '{print $NF}' | awk -F"." '{print $1}'`
	echo $file_chrnum
	cat $file | awk -F "\t" '{OFS="\t"}{printf $1"\t"$4"\t"$5; for(i=6; i<=20;i++){printf "\t"$i} print ""}'  | sort -k 1,1 -k2,2n >  $out_prefix$file_chrnum.gtf
	rm $file
done

