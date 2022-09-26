#!/usr/bin/env

## Merge outliers with variant annotations 

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)


parser = ArgumentParser()
parser$add_argument('--method', help = 'Method') #medz or splicing or ase
parser$add_argument('--outfile', help = 'outfile')
parser$add_argument('--variant_file', help='variant file') #/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv
parser$add_argument('--infile', help = 'input file for medz, splicing, or ase')
args = parser$parse_args()

method = args$method
outfile = as.character(args$outfile)
infile = as.character(args$infile)
this_variant_file<- as.character(args$variant_file)

  medz_outliers = fread(infile)
  colnames(medz_outliers)[1:2] = c('ind','ensg')
  variant_file = fread(paste0("zcat -f ",this_variant_file),header=T,fill=TRUE)
  variant_file$ensg = as.character(variant_file$ensg)
  print("VARIANT FILE")
  head(variant_file)
  print("MEDZ OUTLIERS")
  head(medz_outliers)

  medz_variants = merge(medz_outliers,variant_file,by=c('ind','ensg'),all.x=T)
  write.table(medz_variants,file=outfile,sep='\t',quote=F,row.names=F)

