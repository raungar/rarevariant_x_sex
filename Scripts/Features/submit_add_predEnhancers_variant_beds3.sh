#!/bin/bash
#SBATCH --job-name=add_features
#SBATCH --cpus-per-task=2
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=03:30:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output='/oak/stanford/groups/smontgom/raungar/Sex/Jobs/feature_sexgtex_%j.out'
#SBATCH --error='/oak/stanford/groups/smontgom/raungar/Sex/Jobs/feautre_sexgtex_%j.err'

sh /oak/stanford/groups/smontgom/raungar/Sex/Scripts/Features/add_predEnhancers_variant_beds.sh /oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/variantBeds/individuals/HallLabSV_hg38/*_HallLabSV.bed.gz /oak/stanford/groups/smontgom/raungar/Sex/Files/AF_HallLabSV.bed /oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/bySite/HallLabSV
