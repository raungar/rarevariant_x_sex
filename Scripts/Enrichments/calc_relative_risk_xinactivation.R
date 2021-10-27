#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--xinact_file"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--max_maf"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--gtf_code_file"), type = 'character', default = NULL, help = "gtf_code_file preprocessing"),
                make_option(c("--min_maf"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--nphen"), type = 'numeric', default = NULL, help = "min tissues"),

                make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),             
                make_option(c("--out_rdata_relative"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
                make_option(c("--cadd_min"), type = 'numeric', default = NULL, help = "min cadd score")
        )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
xinact_file <- as.character(opt$xinact_file)
out_rdata_relative <- as.character(opt$out_rdata_relative)
zscore <- as.numeric(opt$zscore)
cadd_min <- as.numeric(opt$cadd_min)
min_maf<-as.numeric(opt$min_maf)
max_maf <- as.numeric(opt$max_maf)
sex <- as.character(opt$sex)
nphen <- as.numeric(opt$nphen)

gtf_code_file<-as.character(opt$gtf_code_file)
# zscore<-2
 # infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/x_outlier_noglobal_medz_varAnnot_zthresh3_nphen5_m_CADDtypesGQ5BlacklistRemovedALL_linc_prot.txt.gz"
#  #out_rdata_relative<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_z3_x_f.xci.RData"
#  xinact_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact.tsv"
# zscore<-3
# red_row="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/collapsed_outliers_rvs_zthresh3_nphen5_x_m.RData"


if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

   
### expression
print('Reading expression')
print(paste0("reading in: ",infile))

exp_data = fread(infile,data.table=F)
#Df is number of tisseus
colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","chr","start","end","vartype","sex",
                      "gtex_maf","gnomad_maf","use_maf","genetype","numrv","cadd_raw","cadd_phred")

exp_data$gene_id_red<-sapply(strsplit(exp_data$ensg,"\\."), "[[",1)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
  this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
  this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_maf)){"common"}
  else if (this_maf>=min_maf & this_maf<max_maf){
    if(this_cadd_phred>=cadd_min){"rare"}
    else if(this_cadd_phred<cadd_min){"common"}
    else (stop("ERROR: VARIANT NOT MAKING SENSE"))
    }
  else if (this_maf>=max_maf){"common"}
  else (stop("ERROR: VARIANT NOT MAKING SENSE"))
})
exp_data$has_variant<-has_variant
gtf_code<-fread(gtf_code_file,header=F)
genetype_dic=gtf_code$V2
names(genetype_dic)<-gtf_code$V1
exp_data$variant_cat<-genetype_dic[exp_data$ensg] 

xinact<-fread(xinact_file,data.table=F)
xinact$GeneID_red<-sapply(strsplit(xinact$`Gene ID`,"\\."), "[[",1)

exp_data_xci<-merge(exp_data,xinact[c("Combined XCI status", "GeneID_red")], by.x="gene_id_red",by.y="GeneID_red")
colnames(exp_data_xci)[colnames(exp_data_xci)=="Combined XCI status"]<-"XCI_STATUS"

print("filter")
### get relative risk per category

exp_outliers = exp_data_xci %>% dplyr::filter(abs(MedZ) >= zscore)  %>% dplyr::filter(as.numeric(Df)>=nphen)
exp_controls = dplyr::filter(exp_data_xci, (abs(MedZ) < zscore) | (abs(MedZ) >= zscore & as.numeric(Df)<nphen))







print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
xci_status = unique(exp_data_xci$XCI_STATUS)
for (xci in na.omit(xci_status)) {
   print(xci)

    exp_nn = nrow(exp_controls %>% dplyr::filter(XCI_STATUS == xci) %>% dplyr::filter( has_variant != "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_ny = nrow(exp_controls %>% dplyr::filter(XCI_STATUS == xci) %>% dplyr::filter(has_variant == "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_yn = nrow(exp_outliers %>% dplyr::filter(XCI_STATUS == xci) %>% dplyr::filter(has_variant != "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_yy = nrow(exp_outliers %>% dplyr::filter(XCI_STATUS == xci) %>% dplyr::filter(has_variant == "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
   print(paste0("building table with: ",exp_nn," and ",exp_ny," and ",exp_yn," and ",exp_yy))

   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
   err = epitab(exptable, method = 'riskratio')
   risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   XCI_STATUS = xci,
                                   exp_nn=exp_nn,
                                   exp_ny=exp_ny,
                                   exp_yn=exp_yn,
                                   exp_yy=exp_yy,
                                   num_outliers=nrow(exp_outliers%>% dplyr::filter(XCI_STATUS == xci)),
                                   sex=sex,z=zscore,nphen=nphen,cadd=cadd_min,
                                   Type = 'XCI'))   
}
 
risks = risks %>% arrange(by=Risk) 
risks$XCI_STATUS = factor(risks$XCI_STATUS, levels=unique(risks$XCI_STATUS))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))



write.csv(risks,  file=out_rdata_relative,quote=F,sep="\t")
print(paste0("COMPLETE: ",out_rdata_relative))
# save(all_coefs, file=out_rdata_continuous)

