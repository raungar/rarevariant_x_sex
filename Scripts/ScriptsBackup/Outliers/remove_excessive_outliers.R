library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
library(tidyverse)
library(optparse)


option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
   make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
  make_option(c("--max_outliers"), type = 'numeric', default = NULL, help = "min cadd score")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
outfile <- as.character(opt$outfile)
max_outliers <- as.numeric(opt$max_outliers)

# infile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8/Current/outliers_noglobal_medz_zthresh3_nphen3_x_f.txt.gz"
#infile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/Current/outliersSkin-Sun-Exposed-Lower-leg_noglobal_medz_zthresh3_nphen1_x_both.txt.gz"
 #max_outliers=3


outliers<-fread(infile)
#ind     ensg
colnames(outliers)[c(1,2)]<-c("Ind","Gene")
outliers_red<-(outliers) %>% group_by(Gene, Y) %>% mutate(n_genes=n()) %>% 
  mutate(filter_me=ifelse(n_genes>=max_outliers & Y=="outlier","filter","ok")) %>% 
  dplyr::filter(filter_me=="ok") %>% select(-n_genes,-filter_me)
fwrite(outliers_red,  file=outfile,quote=F,sep="\t")

