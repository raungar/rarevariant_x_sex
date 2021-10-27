#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)

# 
# option_list = list(
#   make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
#   make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
#   make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
#   make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
#   make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF"),
#   make_option(c("--beta_min"), type = 'numeric', default = NULL, help = "beta_min"),
#   make_option(c("--tissue"), type = 'character', default = NULL, help = "tissue")
# )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
outfile <- as.character(opt$outfile)
tissue <- as.character(opt$tissue)
zscore <- as.numeric(opt$zscore)
min_maf<-as.numeric(opt$min_maf)
max_maf<-as.numeric(opt$max_maf)
beta_min<-as.numeric(opt$beta_min)

# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/SexDeg/aut_beta0.111_SKINS_z3_nphen5_f,_linc_prot.txt.gz"
# #infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen4_noglobal_medz_varAnnot_x_m.txt"
# zscore<-3
# min_maf<-0
# max_maf<-0.01
# beta_min<-0.111
z=3
nphen=5
  
for(this_sex in sex)
  ### expression
  print('Reading expression')
  print(paste0("reading in: ",infile))
  exp_data = fread(infile,data.table=F)
  
  print("relative risk")
  ### Relative risk
  risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
  sexdeg_types<-c("all","male_biased","female_biased")
  beta_min=0.111
  for(this_sexdeg in sexdeg_types){
    if(this_sexdeg=="all"){
      print("all sexDEGs")
      exp_data_sexdeg=exp_data %>%dplyr::filter(abs(beta)>=0.111)
      exp_data_nonsexdeg=exp_data %>%dplyr::filter(abs(beta)<0.111)
      
    }else if(this_sexdeg=="male_biased"){
      print("male sexDEGs")
      exp_data_sexdeg<-exp_data %>% dplyr::filter(beta<=(-1*beta_min))
      exp_data_nonsexdeg<-exp_data%>% dplyr::filter(beta>(-1*beta_min))
    }else if(this_sexdeg=="female_biased"){
      print("female sexDEGs")
      exp_data_sexdeg<-exp_data%>% dplyr::filter(beta>=beta_min)
      exp_data_nonsexdeg<-exp_data%>% dplyr::filter(beta<beta_min)
    }
    exp_data_sexdeg_outlier=exp_data_sexdeg%>%dplyr::filter(abs(MedZ)>=3)%>%nrow()
    exp_data_sexdeg_nonoutlier=exp_data_sexdeg%>%dplyr::filter(abs(MedZ)<3)%>%nrow()
    exp_data_nonsexdegoutlier=exp_data_nonsexdeg%>%dplyr::filter(abs(MedZ)>=3)%>%nrow()
    exp_data_nonsexdeg_nonoutlier=exp_data_nonsexdeg%>%dplyr::filter(abs(MedZ)<3)%>%nrow()
    exptable<-rbind(c(exp_data_sexdeg_outlier,exp_data_sexdeg_nonoutlier),c(exp_data_nonsexdegoutlier,exp_data_nonsexdeg_nonoutlier))
    err = epitab(exptable, method = 'riskratio')
    risks = rbind(risks,data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                    sexdeg=this_sexdeg,
                                   num_sexdegs=nrow(exp_data_sexdeg),
                                   num_outliers=exp_data_nonsexdegoutlier+exp_data_sexdeg_outlier),
                                  sex=this_sexdeg,z=z,nphen=nphen
                  )
    
    
  }



