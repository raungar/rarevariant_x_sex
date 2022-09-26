#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
                make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF"),
                make_option(c("--cadd_min"), type = 'numeric', default = NULL, help = "min cadd score")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
outfile<-as.character(opt$outfile)
cadd_min <- as.numeric(opt$cadd_min)
sex <- as.character(opt$sex)
min_maf<-as.numeric(opt$min_maf)
max_maf<-as.numeric(opt$max_maf)
# #
# infile<-"/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Outliers/outliers_noglobal_varAnnot_medz_zthresh2_nphen3_aut_f_beta0.111_cadd15_aut.txt.gz"
# # infile="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8/OutliersAndRVs/all_outliers_noglobal_medz_varAnnot_zthresh3_nphen5_f_CADDtypesALL_linc_prot.txt.gz"
# gtf_code_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/x_proteincoding_lncrna.gtf"
#my_zscore<-2
# min_maf<-0.1
# nphen=3
# max_maf<-0.2
# sex="f"
# cadd_min=15
# 
# if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","AbsBetaSex","MaxOutlierBetaSex","MinOutlierBetaSex","AbsOutlierBetaSex",
                      "chr","chrNum","start","end","vartype","sex",
                     "gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","numrv")

exp_data$cadd_phred<-as.numeric(exp_data$cadd_phred)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0

print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
  this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
  this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
  #print(paste0("MAF=",this_maf,", cadd=",this_cadd_phred))
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_maf)){"common"}
  else if (this_maf>=min_maf & this_maf<max_maf){
    if(this_cadd_phred>=cadd_min){"rare"}
    else if(this_cadd_phred<cadd_min){"common"}
    else (stop("ERROR: VARIANT NOT MAKING SENSE (at cadd level)"))
  }
  else if (this_maf>=max_maf){"common"}
  else if (this_maf<min_maf){"common"}
  else (stop("ERROR: VARIANT NOT MAKING SENSE (at maf level)"))
})
exp_data$has_variant<-has_variant
write.csv(exp_data,  file=out_exp,quote=F,sep="\t")

