#!/bin/bash
#SBATCH --job-name=add_features
#SBATCH --cpus-per-task=4
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=08:59:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/features_sexgtex_%j.out"
#SBATCH --error="/oak/stanford/groups/smontgom/raungar/Sex/Jobs/features_sexgtex_%j.err"

Scripts/Features/add_features_variant_beds.sh /oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/variantBeds/individualsRA/aut_GTEX-111CU_HallLabSV.bed.gz \
/oak/stanford/groups/smontgom/raungar/Sex/Files/AF_HallLabSV.bed /oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/bySite /oak/stanford/groups/smontgom/raungar/Sex/Files/annotations/epigenomicsRoadmap \
/oak/stanford/groups/smontgom/raungar/Sex/Files/annotations/TFBS_Pouya \
/oak/stanford/groups/smontgom/raungar/Sex/Files/annotations/epigenomicsRoadmap /oak/stanford/groups/smontgom/raungar/Sex/Files/sorted.consolidated.annotations.bed.gz
