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
                make_option(c("--maf_max"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--gtf_code_file"), type = 'character', default = NULL, help = "gtf_code_file preprocessing"),
                make_option(c("--maf_min"), type = 'character', default = NULL, help = "path of x inactivation file"),
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
my_zscore=zscore
cadd_min <- as.numeric(opt$cadd_min)
min_maf<-as.numeric(opt$maf_min)
max_maf <- as.numeric(opt$maf_max)
sex <- as.character(opt$sex)
nphen <- as.numeric(opt$nphen)

gtf_code_file<-as.character(opt$gtf_code_file)
# zscore<-2
 # infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/x_outlier_noglobal_medz_varAnnot_zthresh3_nphen5_m_CADDtypesGQ5BlacklistRemovedALL_linc_prot.txt.gz"
#  #out_rdata_relative<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_z3_x_f.xci.RData"
# xinact_file<-"/Volumes/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact.tsv"
# zscore<-3
# red_row="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/collapsed_outliers_rvs_zthresh3_nphen5_x_m.RData"


if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

   
### expression
print('Reading expression')
print(paste0("reading in: ",infile))

exp_data = fread(infile,data.table=F)
#Df is number of tisseus
if("Tissue" %in% colnames(exp_data)){
  colnames(exp_data)<-c("ind","ensg","MedZ","Df","Y","Tissue","chr","chrNum","start","end","vartype","sex","var_location",
                        "gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","geno","veptype","vepconsq","numrv")
}else if(ncol(exp_data)<22){
  colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","chr","chrNum","start","end","vartype","sex","var_location",
                        "gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","geno","numrv")
  
}else{
  colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","chr","chrNum","start","end","vartype","sex","var_location",
                        "gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","geno","veptype","vepconsq","numrv")
}
exp_data$cadd_phred<-as.numeric(exp_data$cadd_phred)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0
# exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

print("filter for MAF")
#choose rare/common
has_variant<-apply(exp_data,1,function(x){
  this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
  # print(this_maf)
  this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
  this_geno<-as.numeric(x[which(colnames(exp_data)=="geno")])
  # print(this_cadd_phred)
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_maf)){
    "common"
  } else if (this_maf>=min_maf & this_maf<max_maf){
    # print("HI")
    if(this_cadd_phred>=cadd_min & this_geno>=1 & !is.na(this_geno)){"rare"}
    else if(this_cadd_phred<cadd_min){"common"}
    else if(this_cadd_phred>=cadd_min & (this_geno<1 | is.na(this_geno))){"common"}
    else (stop("ERROR: VARIANT NOT MAKING SENSE"))
  } else if (this_maf>=max_maf){
    "common"
  }else if(this_maf<min_maf){
    "common"  
  }else (stop("ERROR: VARIANT NOT MAKING SENSE"))
})
exp_data$has_variant<-has_variant


gtf_code<-fread(gtf_code_file,header=F)
genetype_dic=gtf_code$V2
names(genetype_dic)<-gtf_code$V1
exp_data$variant_cat<-genetype_dic[exp_data$ensg] 

xinact<-fread(xinact_file,data.table=F)
xinact$GeneID_red<-sapply(strsplit(xinact$`Gene ID`,"\\."), "[[",1)

exp_data_xci<-merge(exp_data,xinact[c("Combined XCI status", "GeneID_red")], by.x="ensg",by.y="GeneID_red")
colnames(exp_data_xci)[colnames(exp_data_xci)=="Combined XCI status"]<-"CATEGORY"
print(head(exp_data_xci))
print("filter")
### get relative risk per category
### get relative risk per category
exp_controls = exp_data_xci %>% dplyr::filter(Y=="control")
exp_outliers = exp_data_xci %>% dplyr::filter(Y=="outlier")
exp_outliers_over = exp_data_xci %>% dplyr::filter(Y=="outlier") %>% dplyr::filter(MedZ>0)#was just expdata
exp_controls_over = exp_data_xci %>% dplyr::filter(Y=="control") %>% dplyr::filter(MedZ>0)#was just expdata
exp_outliers_under = exp_data_xci %>% dplyr::filter(Y=="outlier") %>% dplyr::filter(MedZ<0)#was just expdata
exp_controls_under = exp_data_xci %>% dplyr::filter(Y=="control") %>% dplyr::filter(MedZ<0)#was just expdata

print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Cat = character(),
                   exp_nn=numeric(),exp_ny=numeric(),exp_yn=numeric(),exp_yy=numeric(),
                   sex=character(),z=numeric(),nphen=numeric(),Type=character())
vcats = na.omit(unique(exp_data$variant_cat))
exp_types=c("all","over","under")
veptypes=c(na.omit(unique(unlist(strsplit(exp_data$veptype,",")))),"all")
varlocations=as.character(na.omit(c("all",unique(exp_data$var_location))))
xci_status = c(unique(exp_data_xci$CATEGORY),"all")

for (this_xci in xci_status) {

  
  for(this_var_location in varlocations){
    for(this_veptype in veptypes){
      for(exp_type in exp_types){
        print(paste0("HERE FOR ", this_veptype, " and ",this_var_location, " and ",exp_type, " and ",this_xci))
        if(this_xci!="all"){
          if(exp_type=="all"){
            this_exp_outliers<-exp_outliers%>% dplyr::filter(CATEGORY == this_xci)
            this_exp_controls<-exp_controls%>% dplyr::filter(CATEGORY == this_xci)
            
          }else if(exp_type=="over"){
            this_exp_outliers<-exp_outliers_over%>% dplyr::filter(CATEGORY == this_xci)
            this_exp_controls<-exp_controls_over%>% dplyr::filter(CATEGORY == this_xci)
            
          }else if(exp_type=="under"){
            this_exp_outliers<-exp_outliers_under%>% dplyr::filter(CATEGORY == this_xci)
            this_exp_controls<-exp_controls_under%>% dplyr::filter(CATEGORY == this_xci)
          }else{stop("ERROR: INVALID EXP TYPE")}
        }else{
          if(exp_type=="all"){
            this_exp_outliers<-this_exp_outliers
            this_exp_controls<-this_exp_controls
            
          }else if(exp_type=="over"){
            this_exp_outliers<-exp_outliers_over
            this_exp_controls<-exp_controls_over
            
          }else if(exp_type=="under"){
            this_exp_outliers<-exp_outliers_under
            this_exp_controls<-exp_controls_under
          }else{stop("ERROR: INVALID EXP TYPE")}
        }
        if(this_var_location=="all" & this_veptype=="all"){
          this_exp_outliers<-this_exp_outliers
          this_exp_controls<-this_exp_controls
        }else if(this_var_location!="all" & this_veptype=="all"){
          this_exp_outliers<-this_exp_outliers%>%dplyr::filter(var_location==this_var_location)
          this_exp_controls<-this_exp_controls%>%dplyr::filter(var_location==this_var_location)
        }else if(this_var_location=="all" & this_veptype!="all"){
          #dplyr::filter(df, !grepl("RTB",TrackingPixel))
          # this_exp_outliers<-this_exp_outliers%>% filter(map_lgl(str_split(veptype, ","),
          #                                                        ~ any(str_detect(.x, this_veptype))))        
          this_exp_outliers<-this_exp_outliers%>% dplyr::filter(grepl(this_veptype,veptype))
          this_exp_controls<-this_exp_controls%>%  dplyr::filter(grepl(this_veptype,veptype))
        }else{
          this_exp_outliers<-this_exp_outliers%>%  dplyr::filter(grepl(this_veptype,veptype))%>%
            dplyr::filter(var_location==this_var_location) 
          this_exp_controls<-this_exp_controls%>% dplyr::filter(grepl(this_veptype,veptype)) %>%dplyr::filter(var_location==this_var_location)
        }
        if(nrow(this_exp_outliers)==0 | nrow(this_exp_controls)==0){print("no outliers or controls"); next}
        print(head(this_exp_outliers))
        
        for (vcat in (vcats)) {
          print(vcat)
          print(cadd_min)
          exp_nn = nrow(this_exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
          exp_ny = nrow(this_exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
          exp_yn = nrow(this_exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
          exp_yy = nrow(this_exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare")) # %>% dplyr::filter(cadd_phred>cadd_min))
          print(paste0("the following: ",exp_nn," and ",exp_ny," and ",exp_yn," and ",exp_yy))
          if(exp_nn==0 & exp_ny==0){print("NO NON OUTLIERS")}
          else if(exp_nn==1 | exp_ny == 1){print("don't swap for no reason.")}
          else if(exp_nn==0){exp_nn=1;exp_ny=exp_ny-1}
          else if(exp_ny==0){exp_ny=1;exp_nn=exp_nn-1}
          else{print("No zeros, no worries")}
          if(exp_yn==0 & exp_yy==0){print("NO OUTLIERS")}
          else if(exp_yn==1 | exp_yy == 1){print("don't swap for no reason.")}
          else if(exp_yn==0){exp_yn=1;exp_yy=exp_yy-1}
          # else if(exp_yy==0){exp_yy=1;exp_yn=exp_yn-1}
          else{print("No zeros, no worries")}
          #exp_nn = nrow(this_exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare"))
          #exp_ny = nrow(this_exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant != "rare"))
          #exp_yn = nrow(exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant== "rare"))
          #exp_yy = nrow(exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant != "rare"))
          exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
          print(exptable)
          err = epitab(exptable, method = 'riskratio')
          risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                          Lower = err$tab[2,6],
                                          Upper = err$tab[2,7],
                                          Pval = err$tab[2,8],
                                          Cat = this_xci,
                                          exp_nn=exp_nn,
                                          exp_ny=exp_ny,
                                          exp_yn=exp_yn,
                                          exp_yy=exp_yy,
                                          var_location=this_var_location,veptype=this_veptype,
                                          num_outliers=nrow(this_exp_outliers%>% dplyr::filter(variant_cat == vcat)),
                                          sex=sex,z=my_zscore,nphen=nphen,cadd=cadd_min,
                                          exp_type=exp_type,Type = 'XCI'))   
        }
        
        exp_nn_all=nrow(dplyr::filter(this_exp_controls,has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_ny_all=nrow(dplyr::filter(this_exp_controls,has_variant == "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_yn_all=nrow(dplyr::filter(this_exp_outliers,has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_yy_all=nrow(dplyr::filter(this_exp_outliers,has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        print(paste0("the following: ",exp_nn_all," and ",exp_ny_all," and ",exp_yn_all," and ",exp_yy_all))
        if(exp_nn_all==0 & exp_ny_all==0){print("NO NON OUTLIERS")}
        else if(exp_nn_all==1 | exp_ny_all == 1){print("don't swap for no reason.")}
        else if(exp_nn_all==0){exp_nn_all=1;exp_ny_all=exp_ny_all-1}
        else if(exp_ny_all==0){exp_ny_all=1;exp_nn_all=exp_nn_all-1}
        else{print("No zeros, no worries")}
        if(exp_yn_all==0 & exp_yy_all==0){print("NO OUTLIERS")}
        else if(exp_yn_all==1 | exp_yy_all == 1){print("don't swap for no reason.")}
        else if(exp_yn_all==0){exp_yn_all=1;exp_yy_all=exp_yy_all-1}
        else if(exp_yy_all==0){exp_yy_all=1;exp_yn_all=exp_yn_all-1}
        # else{print("No zeros, no worries")}
        #exp_nn_all = nrow(this_exp_outliers %>% dplyr::filter( has_variant == "rare"))
        #exp_ny_all = nrow(this_exp_outliers %>% dplyr::filter(has_variant != "rare"))
        #exp_yn_all = nrow(exp_controls %>% dplyr::filter(has_variant== "rare"))
        #exp_yy_all = nrow(exp_controls  %>% dplyr::filter(has_variant != "rare"))
        exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
        colnames(exptable_all)<-c("non-outliers","outliers")
        rownames(exptable_all)<-c("common","rare")
        err_all = epitab(exptable_all, method = 'riskratio')
        risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                        Lower = err_all$tab[2,6],
                                        Upper = err_all$tab[2,7],
                                        Pval = err_all$tab[2,8],
                                        Cat = this_xci,
                                        exp_nn=exp_nn_all,
                                        exp_ny=exp_ny_all,
                                        exp_yn=exp_yn_all,
                                        exp_yy=exp_yy_all,
                                        var_location=this_var_location,veptype=this_veptype,
                                        num_outliers=nrow(this_exp_outliers),
                                        sex=sex,z=my_zscore,nphen=nphen,cadd=cadd_min,
                                        exp_type=exp_type,Type = 'XCI'))  
        
      }
    }
  }
}
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))



write.csv(risks,  file=out_rdata_relative,quote=F,sep="\t")
print(paste0("COMPLETE: ",out_rdata_relative))
# save(all_coefs, file=out_rdata_continuous)

