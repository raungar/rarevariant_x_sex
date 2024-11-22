#!/bin/bash

vcf="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
output_vcf="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/chrX_GQ.vcf"
bcftools view $vcf --regions chrX > $output_vcf
