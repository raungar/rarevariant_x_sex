#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(argparse)

## Script for combining outlier calls and features then computing enrichments.



parser = ArgumentParser()
parser$add_argument('--dir.suffix', help = 'Suffix of directory name (e.g., v7 or multi_omics).')
parser$add_argument('--outliers.file', help = 'Outlier file.')
parser$add_argument('--outfile', help = 'output file.')
parser$add_argument('--window', help = 'Window name.')
parser$add_argument('--z_thresh', default = 3, help = 'Z-score threshold.')
parser$add_argument('--output.suffix', default = '', help = 'suffix for output file')
parser$add_argument('--scriptdir', default = NULL, help = 'script directory')
parser$add_argument('--RAREDIR', default = NULL, help = 'RAREDIR')
args = parser$parse_args()

baseDir = as.character(args$RAREDIR)
outfile = as.character(args$outfile)
scriptdir_path<-as.character(args$scriptdir)

source(paste0(scriptdir_path,"/enrichment_functions.R"))


#iodir = paste0(baseDir, '/data_', args$dir.suffix, '/outliers/')
iodir = paste0(baseDir,'/enrichments_v8/')
if(!dir.exists(iodir)){dir.create(iodir)}
featdir = paste0(baseDir, '/features_', args$dir.suffix, '/byGene/', args$window, '/')
if(!dir.exists(featdir)){dir.create(featdir)}

## Get files with outliers
z_thresh = as.numeric(args$z_thresh)

## First get medz outliers
#medz_outliers = fread(paste0(baseDir, '/data_v8/outliers/', args$outliers.file)) %>% mutate(Method = 'MEDZ')
medz_outliers = fread(args$outliers.file) %>% mutate(Method = 'MEDZ')

all_logits_top = get.all.enrich(medz_outliers, featdir, counts = TRUE)

## Set factor levels and labels for plotting
maf.levels = unique(all_logits_top$Maf)
maf.labels = sub('MAF', '', maf.levels)
maf.order = order(as.numeric(str_split_fixed(maf.labels, '-', 2)[,1])) # get order by first number
maf.levels = maf.levels[maf.order]
maf.labels = maf.labels[maf.order]
all_logits_top$Maf = factor(all_logits_top$Maf, levels = maf.levels, labels = maf.labels)

#save(all_logits_top, file = paste0(iodir, 'enrichments_', args$window, '_Z', z_thresh, '_', args$output.suffix, '.RData'))
save(all_logits_top, file = outfile)

