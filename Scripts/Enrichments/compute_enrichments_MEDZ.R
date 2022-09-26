#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(argparse)

## Script for combining outlier calls and features then computing enrichments.



parser = ArgumentParser()
parser$add_argument('--featuresdir', help = 'd')
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
featdir<-as.character(args$featuresdir)

source(paste0(scriptdir_path,"/enrichment_functions.R"))

## Produce all enrichments statistics for a given set of outliers (across all mafs and variant types)
## Input: outlier/control data frame, path to features, whether they are count features, whaether to scale features
## Output: Single data frame with all the enrichments from logistic regression
get.all.enrich <- function(outliers, feature_dir, counts = FALSE, scale = TRUE, suffix = NULL) {
    if (is.null(suffix)) {
        if (counts) {
            feature_dir = paste0(feature_dir, '/counts/')
            suffix = '_counts_bygene.txt'
            scale = FALSE
        } else {
            suffix = '_features_bygene.txt'
        }
    }

    ## get mafs and vartypes - assumes the same vartypes for each MAF
    mafs = list.dirs(feature_dir, full.names = FALSE, recursive = FALSE)
    mafs = c('MAF0-1','MAF1-5','MAF5-10','MAF10-25')
    if (!counts) {
#        mafs = 'MAF0-1'
        mafs = mafs[!(mafs == 'counts')] # removes counts dir, if there
    }
    vartypes = unique(list.dirs(paste0(feature_dir, '/', mafs), full.names = FALSE, recursive = FALSE))
    #vartypes = vartypes[which(vartypes != 'HallLabSV')] # removing for v8
    if (!counts) {
        vartypes = vartypes[which(vartypes != 'HallLabSV')] # Remove SVs for features
    }
    vartypes = c('SNPs','HallLabSV') # adding to combine enrichments for SNPs+indels
    #vartypes = c('SNPs')
    ## for all combinations of mafs and vartypes, merge features with outliers
    ## then collapse into a single data frame
    if (!counts) {
	print("NEED TO IMPLEMENT COMPUTE.ENRICHMENTS")
	logit_all=-1
     # logit_all = do.call(rbind,
     #                     mapply(compute.enrichments,
      #                           rep(mafs, length(vartypes)), rep(vartypes, each = length(mafs)),
      #                           MoreArgs = list(outliers, feature_dir, suffix, scale, counts), SIMPLIFY = FALSE))
    } else {
      logit_all = do.call(rbind,
                          mapply(compute.relative.risk,
                                 rep(mafs, length(vartypes)), rep(vartypes, each = length(mafs)),
                                 MoreArgs = list(outliers, feature_dir, suffix, scale, counts), SIMPLIFY = FALSE))
    }
    return(logit_all)
}

#iodir = paste0(baseDir, '/data_', args$dir.suffix, '/outliers/')
iodir = paste0(baseDir,'/enrichments_v8/')
if(!dir.exists(iodir)){dir.create(iodir)}
#featdir = paste0(baseDir, '/features_', args$dir.suffix, '/byGene/', args$window, '/')
#if(!dir.exists(featdir)){dir.create(featdir)}

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

