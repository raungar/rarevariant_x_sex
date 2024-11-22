#!/usr/bin/env Rscript


rm(list = ls())

## R script to process VCFtools output containing the allele counts for each V8 individual at the V8 cis-eQTL calls

require(ggplot2)
require(data.table)
require(dplyr)
require(reshape2)
require(plyr)
library(optparse)

#--------------- FUNCTIONS
get_opt_parser<-function(){
        option_list = list(
                 make_option(c("-d", "--RAREDIR"), type="character", default=NULL, help="RAREDIR directory"),
                 make_option(c("-i", "--inds"), type="character", default=NULL, help="inds file"),
                 make_option(c("-p", "--pos"), type="character", default=NULL, help="pos file"),
                 make_option(c("-g", "--genotypes"), type="character", default=NULL, help="genotypes files"),
                 make_option(c("-o", "--out"), type="character", default=NULL,help="output file")
        )
        opt_parser = OptionParser(option_list=option_list)
        return(opt_parser)

}

#--------------- MAIN
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)

baseDir = as.character(args$RAREDIR) #Sys.getenv('RAREDIR')
inds_str = as.character(args$inds)
pos_str = as.character(args$pos)
genos_str = as.character(args$genotypes)
out_str=as.character(args$out)

dir = paste0(baseDir, '/preprocessing_v8')


# Load the list of individuals, the site position info, and the allele counts
inds = read.table(paste( inds_str, sep = ''), header = F, stringsAsFactors = F)[, 1]
pos = read.table(paste( pos_str, sep = ''), header = F)
genos = t(as.data.frame(fread(paste(genos_str, sep = ''), na.strings = c('-1'))))[-1, ]
#inds = read.table(paste(dir, '/', inds_str, sep = ''), header = F, stringsAsFactors = F)[, 1]
#pos = read.table(paste(dir, '/', pos_str, sep = ''), header = F)
#genos = t(as.data.frame(fread(paste(dir, '/', genos_str, sep = ''), na.strings = c('-1'))))[-1, ]

# Combine the position and genotype info into a single data frame
# Add the individual IDs as headers
out = cbind(pos, genos)
colnames(out) = c('Chrom', 'Pos', inds)

# Remove duplicated positions
out = out %>% mutate(ID = paste(Chrom, '_', Pos, sep = ''))
ids = as.data.frame(table(out$ID))
ids.to.keep = ids %>% filter(Freq == 1)
out = out %>% filter(ID %in% ids.to.keep$Var1) %>% select(-ID)

# Write out the combined data frame
#write.table(out, paste(dir, out_str, sep = ''), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(out, out_str, sep = '\t', col.names = T, row.names = F, quote = F)
