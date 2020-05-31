#!/bin/bash
#SBATCH --job-name=correct_sex
#SBATCH --cpus-per-task=12
#SBATCH --partition=interactive
#SBATCH --account=default
##SBATCH --time=06:59:00
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/correct_sexgtex_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/correct_sexgtex_%j.err"


echo $4
tissue_dir=`echo $4 | awk -F"/" '{print $(NF-1)}'`
Rscript $3 --dir_peer $1 --metadata_file $2 --file_to_correct $4 --outfile $1/$tissue_dir/$5
echo "$1/$tissue_dir/$5"
