#!/usr/bin/env Rscript

## Script to call outliers given a file with Z-scores.
## The two first columns have the gene and tissue or phenotype, and the gene columne is named "Gene".
## Subsequent columns have Z-score data for each individual and the individual ID is the column name.

print("IN R SCRIPT")
## Get path to the GOATs data directories
#dir = Sys.getenv('RAREDIR')

## Load required packages
require(data.table)
require(optparse)
require(reshape2)
require(dplyr)
require(foreach)
require(doMC)
require(robustbase)

## Register parallel backend
registerDoMC(cores = 3)

#--------------- FUNCTIONS

#### Function to pick top outliers per gene
## Input: Data frame of outlier test statistics and the number of phenotypes required per test
## Output: Assignments of outliers and controls for genes with outliers
## Performs BF correction on the gene levels and Bh across genes
pick.outliers <- function(test.stats, nphen, zthresh,medz){
    test.stats = test.stats %>% filter(Df >= nphen) %>%
        group_by(Gene) %>%
        mutate(N = n())
        out = test.stats %>% arrange(desc(abs(MedZ))) %>%
            mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
            ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y)
    return(out)
}

#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, nphen, zthresh,metric) {
    data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
    colnames(data.melted)[3:4] = c('Ind', 'Z')
    head(data.melted)
    test.stats = data.melted %>% group_by(Ind, Gene) %>%
        summarise(MedZ = median(Z, na.rm = T),
                  Df = sum(!is.na(Z)))
    head(test.stats)
    outliers = pick.outliers(test.stats, nphen, zthresh, medz = T)
    return(outliers)
}

#### Function to write the outlier information to file.
## Input: data frame with data to write and output filename
## Output: None, writes directly to file.
write.outliers <- function(outliers, filename) {
    write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
}

#--------------- MAIN

##-- Read command line arguments
option_list = list(
	make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
	make_option(c('--GLOBAL'), type = 'character', default = NA, help = 'global outlier file'),
	make_option(c('--N.PHEN'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
	make_option(c("--outfile"),type="character", default=NULL, help= "output file path"),
	make_option(c("--outdir"),type="character", default=NULL, help= "output directory"),
	make_option(c('--ZTHRESH'), type = 'numeric', default = 2, help = 'threshold for abs(MedZ) outliers')
	)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

nphen = opt$N.PHEN
zthresh = opt$ZTHRESH
globalFile = opt$GLOBAL
outfile=opt$outfile
zscore=opt$Z.SCORES
outdir=opt$outdir

print(paste0("ZSCORES: ", zscore))
##-- Analysis

if (!(dir.exists(outdir))){
	dir.create(outdir)
}

## Read in the normalized data
data = as.data.frame(fread(paste0('zcat -f ', zscore)))


##-- Call outliers using median Z-score
# Filtering global outliers
if (!is.na(globalFile)) {
    goutliers = fread(globalFile, header=F)
    rinds = which(colnames(data) %in% goutliers$V1)
    data = data[,-1*rinds]
}


## Median Z-score
print('MEDZ')
outliers.medz = call.outliers(data, nphen, zthresh,metric = 'medz')
#write.outliers(outliers.medz, paste0(dir, prefix, '.medz.txt'))
write.outliers(outliers.medz, outfile)
