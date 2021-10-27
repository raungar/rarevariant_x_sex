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
                make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
                make_option(c("--nphen"), type = 'numeric', default = NULL, help = "min tissues"),
                make_option(c("--cadd"), type = 'numeric', default = NULL, help = "min cadd"),
                make_option(c("--chr"), type = 'character', default = NULL, help = "chr (x, aut, or 7 etc)"),
                make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--gtf_code_file"), type = 'character', default = NULL, help = "gtf_code_file preprocessing"),
                make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
                make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF"),
                make_option(c("--beta_min"), type = 'numeric', default = NULL, help = "beta_min"),
                make_option(c("--tissue"), type = 'character', default = NULL, help = "tissue")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
outfile <- as.character(opt$outfile)
my_chr <- as.character(opt$chr)
my_cadd <- as.numeric(opt$cadd)
zscore <- as.numeric(opt$zscore)
nphen <- as.numeric(opt$nphen)
my_tissue <- as.numeric(opt$tissue)
sex <- as.character(opt$sex)
min_maf<-as.numeric(opt$min_maf)
max_maf<-as.numeric(opt$max_maf)
beta_min<-as.numeric(opt$beta_min)
gtf_code_file<-as.character(opt$gtf_code_file)


print(gtf_code_file)
#infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Combined/aut_beta0.111_SKINS_z3_nphen5_f_linc_prot.txt.gz"
# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/aut_beta0.111_NERVET_m_linc_prot.txt.gz"
# #infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen4_noglobal_medz_varAnnot_x_m.txt"
#gtf_code_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/both_proteincoding_lncrna.gtf"
#infile<-"/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/x_beta0.111_BREAST_m_cadd15_linc_prot.txt.gz"
# zscore<-3
# min_maf<-0
# max_maf<-0.01
# beta_min<-0.111
# chr="x"
# cadd=15
#tissue="BREAST"

### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
# colnames(exp_data)<-c("ind","ensg","z","chr.x","beta", "chr.y","pos","vartype","sex",
#                       "ref","alt","bothvars", "gtex_maf","gnomad_maf_both","gnomad_maf_m","gnomad_maf_f","genetype","mafdiff","numrv")

gtf_code<-fread(gtf_code_file,header=F)
genetype_dic=gtf_code$V2
names(genetype_dic)<-gtf_code$V1


exp_data$cadd_phred<-as.numeric(exp_data$cadd_phred)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0
# exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))
print(paste0("HELLO MY CADD IS ",my_cadd))
print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
  this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
  this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_maf)){"common"}
  else if (this_maf>=min_maf & this_maf<max_maf){
    if(this_cadd_phred>=my_cadd){"rare"}
    else if(this_cadd_phred<my_cadd){"common"}
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

print("filter for z score")
print(head(exp_data))
exp_outliers = exp_data %>% dplyr::filter(z>=zscore)
exp_controls = exp_data %>% dplyr::filter(z<zscore)




print("relative risk")



print("relative risk  all sexdegs")
### Relative risk
##rr=Noutliers-sexDEG-rv/NsexDEG-outlier
#    -------------------------
#     N-nonoutliers-nonsexDEG-rv/NsexDEG-nonoutlier
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Cat = character(),
                   sexdeg=character(),exp_nn=numeric(),exp_ny=numeric(),exp_yn=numeric(),exp_yy=numeric(),cadd=numeric(),chr=character(),
                   tissue=character(),num_outliers=character(),sex=character(),z=numeric(),nphen=numeric(),Type=character())
vcats = na.omit(unique(exp_data$variant_cat))
sexdeg_types<-c("all","male_biased","female_biased")
for (this_sexdeg in sexdeg_types){
   if(this_sexdeg=="all"){
      print("all sexDEGs")
      this_controlset<-exp_controls %>% dplyr::filter(abs(beta)>=beta_min)
      this_outlierset<-exp_outliers%>% dplyr::filter(abs(beta)>=beta_min)
   }else if(this_sexdeg=="male_biased"){
      print("male sexDEGs")
      this_controlset<-exp_controls %>% dplyr::filter(beta<=(-1*beta_min))
      this_outlierset<-exp_outliers%>% dplyr::filter(beta<=(-1*beta_min))
   }else if(this_sexdeg=="female_biased"){
      print("female sexDEGs")
      this_controlset<-exp_controls%>% dplyr::filter(beta>=beta_min)
      this_outlierset<-exp_outliers%>% dplyr::filter(beta>=beta_min)
   }else{stop("invalid sexdeg type")}
    for (vcat in (vcats)) {
         print(vcat)
         
         exp_nn = nrow(this_controlset %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare") )
         exp_ny = nrow(this_controlset %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))
         exp_yn = nrow(this_outlierset %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))
         exp_yy = nrow(this_outlierset %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare"))
         
         exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
         print(exptable)
	     if(exp_nn == 0 & exp_ny==0 & exp_yn==0 & exp_yy==0){
	       print("all zero")
         	risks = rbind(risks, data.frame(Risk = NA, #err$tab[2,5],
                                         Lower = NA, #err$tab[2,6],
                                         Upper = NA, #err$tab[2,7],
                                         Pval = NA, #err$tab[2,8],
                                         Cat = vcat,
                                         sexdeg=this_sexdeg,
                                         exp_nn=exp_nn,
                                         exp_ny=exp_ny,
                                         exp_yn=exp_yn,
                                         exp_yy=exp_yy,
                                         cadd=my_cadd,chr=my_chr,
                                         tissue=my_tissue,
                                         num_outliers=nrow(this_outlierset%>% dplyr::filter(variant_cat == vcat)),
                                         sex=sex,z=zscore,nphen=nphen,
                                         Type = 'RR_RV_ALLSexDEGs'))   
		
	 }else {
	   print("not all zeros")
	         err = epitab(exptable, method = 'riskratio')
        	 risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                         Lower = err$tab[2,6],
                                         Upper = err$tab[2,7],
                                         Pval = err$tab[2,8],
                                         Cat = vcat,
                                         sexdeg=this_sexdeg,
                                         exp_nn=exp_nn,
                                         exp_ny=exp_ny,
                                         exp_yn=exp_yn,
                                         exp_yy=exp_yy,
                                         cadd=my_cadd,chr=my_chr,
                                         tissue=my_tissue,
                                         num_outliers=nrow(this_outlierset%>% dplyr::filter(variant_cat == vcat)),
                                         sex=sex,z=zscore,nphen=nphen,
                                         Type = 'RR_RV_ALLSexDEGs'))   
	}
      }
      print("all")
      exp_nn_all=nrow(dplyr::filter(this_controlset,has_variant != "rare"))
      exp_ny_all=nrow(dplyr::filter(this_controlset,has_variant == "rare"))
      exp_yn_all=nrow(dplyr::filter(this_outlierset,has_variant != "rare"))
      exp_yy_all=nrow(dplyr::filter(this_outlierset,has_variant == "rare"))
      exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
      print(exptable)
       if(exp_nn_all == 0 & exp_ny_all==0 & exp_yn_all==0 & exp_yy_all==0){
	      risks = rbind(risks, data.frame(Risk = NA, #err_all$tab[2,5],
                                      Lower = NA, #err_all$tab[2,6],
                                      Upper = NA, #err_all$tab[2,7],
                                      Pval = NA, #err_all$tab[2,8],
                                      Cat = "all",
                                      sexdeg=this_sexdeg,
                                      exp_nn=exp_nn_all,
                                      exp_ny=exp_ny_all,
                                      exp_yn=exp_yn_all,
                                      exp_yy=exp_yy_all,
                                      cadd=my_cadd,chr=my_chr,
                                      tissue=my_tissue,
                                      num_outliers=nrow(this_outlierset),
                                      sex=sex,z=zscore,nphen=nphen,
                                      Type = 'RR_RV_ALLSexDEGs'))   

	}else {
	      err_all = epitab(exptable_all, method = 'riskratio')
	      risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                      Lower = err_all$tab[2,6],
                                      Upper = err_all$tab[2,7],
                                      Pval = err_all$tab[2,8],
                                      Cat = "all",
                                      sexdeg=this_sexdeg,
                                      exp_nn=exp_nn_all,
                                      exp_ny=exp_ny_all,
                                      exp_yn=exp_yn_all,
                                      exp_yy=exp_yy_all,
                                      cadd=my_cadd,chr=my_chr,
                                      tissue=my_tissue,
                                      num_outliers=nrow(this_outlierset),
                                      sex=sex,z=zscore,nphen=nphen,
                                      Type = 'RR_RV_ALLSexDEGs'))  
	} 
      risks = risks %>% arrange(by=Risk) 
      risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
}


#rr only for outliers
##rr=Noutliers-sexDEG-rv/NsexDEG-outlier
#    -----------------------------------------
#     Noutliers-nonsexDEG-rv/nonsexDEG-outlier
for (this_sexdeg in sexdeg_types){
  if(this_sexdeg=="all"){
    print("all sexDEGs")
    sexdeg<-exp_outliers %>% dplyr::filter(abs(beta)>=beta_min)
    non_sexdeg<-exp_outliers%>% dplyr::filter(abs(beta)<beta_min)
  }else if(this_sexdeg=="male_biased"){
    print("male sexDEGs")
    sexdeg<-exp_outliers %>% dplyr::filter(beta<=(-1*beta_min))
    non_sexdeg<-exp_outliers%>% dplyr::filter(beta>(-1*beta_min))
  }else if(this_sexdeg=="female_biased"){
    print("female sexDEGs")
    sexdeg<-exp_outliers%>% dplyr::filter(beta>=beta_min)
    non_sexdeg<-exp_outliers%>% dplyr::filter(beta<beta_min)
  }else{stop("invalid sexdeg type")}
  for (vcat in (vcats)) {
    print(vcat)
    exp_nn = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))
    exp_ny = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare"))
    exp_yn = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant != "rare"))
    exp_yy = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare") )
    exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
    print(exptable)
      if(exp_nn == 0 & exp_ny==0 & exp_yn==0 & exp_yy==0){
	    risks = rbind(risks, data.frame(Risk = NA , # err$tab[2,5],
                                    Lower = NA , # err$tab[2,6],
                                    Upper = NA , # err$tab[2,7],
                                    Pval =NA ,# err$tab[2,8],
                                    Cat = vcat,
                                    exp_nn=exp_nn,
                                    exp_ny=exp_ny,
                                    exp_yn=exp_yn,
                                    exp_yy=exp_yy,
                                    cadd=my_cadd,chr=my_chr,
                                    sexdeg=this_sexdeg,
                                    tissue=my_tissue,
                                    num_outliers=nrow(sexdeg),
                                    sex=sex,z=zscore,nphen=nphen,
                                    Type = 'RR_RV_ALLOutliers'))   

    }else {  
	    err = epitab(exptable, method = 'riskratio')
	    risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                    Lower = err$tab[2,6],
                                    Upper = err$tab[2,7],
                                    Pval = err$tab[2,8],
                                    Cat = vcat,
                                    sexdeg=this_sexdeg,
                                    exp_nn=exp_nn,
                                    exp_ny=exp_ny,
                                    exp_yn=exp_yn,
                                    exp_yy=exp_yy,
                                    cadd=my_cadd,chr=my_chr,
                                    tissue=my_tissue,
                                    num_outliers=nrow(sexdeg),
                                    sex=sex,z=zscore,nphen=nphen,
                                    Type = 'RR_RV_ALLOutliers'))   
	}
  }
  print("all")
  exp_nn_all=nrow(dplyr::filter(non_sexdeg,has_variant != "rare"))
  exp_ny_all=nrow(dplyr::filter(non_sexdeg,has_variant == "rare"))
  exp_yn_all=nrow(dplyr::filter(sexdeg,has_variant != "rare"))
  exp_yy_all=nrow(dplyr::filter(sexdeg,has_variant == "rare"))
  exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
  print(exptable)
  if(exp_nn_all == 0 & exp_ny_all==0 & exp_yn_all==0 & exp_yy_all==0){
		  risks = rbind(risks, data.frame(Risk = NA, #err_all$tab[2,5],
                                  Lower =  NA, #err_all$tab[2,6],
                                  Upper =  NA, # err_all$tab[2,7],
                                  Pval =  NA, #err_all$tab[2,8],
                                  Cat = "all",
                                  sexdeg=this_sexdeg,
                                  exp_nn=exp_nn_all,
                                  exp_ny=exp_ny_all,
                                  exp_yn=exp_yn_all,
                                  exp_yy=exp_yy_all,
                                  cadd=my_cadd,chr=my_chr,
                                  tissue=my_tissue,
                                  num_outliers=nrow((sexdeg)),
                                   sex=sex,z=zscore,nphen=nphen,
                                  Type = 'RR_RV_ALLOutliers'))   
	
  }else{
	  err_all = epitab(exptable_all, method = 'riskratio')
	  risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                  Lower = err_all$tab[2,6],
                                  Upper = err_all$tab[2,7],
                                  Pval = err_all$tab[2,8],
                                  Cat = "all",
                                  sexdeg=this_sexdeg,
                                  exp_nn=exp_nn_all,
                                  exp_ny=exp_ny_all,
                                  exp_yn=exp_yn_all,
                                  exp_yy=exp_yy_all,
                                  cadd=my_cadd,chr=my_chr,
                                  tissue=my_tissue,
                                  num_outliers=nrow((sexdeg)),
                                   sex=sex,z=zscore,nphen=nphen,
                                  Type = 'RR_RV_ALLOutliers'))   
	} 
  risks = risks %>% arrange(by=Risk) 
  risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
}

# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))
# 
# write.csv(risks,  file=outfile_sexdeg,quote=F,sep="\t")
# print(paste0("COMPLETE: ",outfile_sexdeg))
#
print("relative risk in outliers no RV")

### Relative risk
##rr=Noutliers-sexDEG/NsexDEG
#    -------------------------
#     Noutliers-nonsexDEG/nonsexDEG
sexdeg_types<-c("all","male_biased","female_biased")
# beta_min=0.111
vcats = na.omit(unique(exp_data$variant_cat))

for(this_sexdeg in sexdeg_types){
  if(this_sexdeg=="all"){
    print("all sexDEGs")
    exp_data_sexdeg=exp_data %>%dplyr::filter(abs(beta)>=beta_min)
    exp_data_nonsexdeg=exp_data %>%dplyr::filter(abs(beta)<beta_min)
    
  }else if(this_sexdeg=="male_biased"){
    print("male sexDEGs")
    exp_data_sexdeg<-exp_data %>% dplyr::filter(beta<=(-1*beta_min))
    exp_data_nonsexdeg<-exp_data%>% dplyr::filter(beta>(-1*beta_min))
  }else if(this_sexdeg=="female_biased"){
    print("female sexDEGs")
    exp_data_sexdeg<-exp_data%>% dplyr::filter(beta>=beta_min)
    exp_data_nonsexdeg<-exp_data%>% dplyr::filter(beta<beta_min)
  }  
    for (vcat in (vcats)) {
      print(vcat)
      exp_data_sexdeg_outlier=exp_data_sexdeg%>% dplyr::filter(variant_cat == vcat)%>%dplyr::filter(abs(z)>=zscore)%>%nrow()
      exp_data_sexdeg_nonoutlier=exp_data_sexdeg%>% dplyr::filter(variant_cat == vcat)%>%dplyr::filter(abs(z)<zscore)%>%nrow()
      exp_data_nonsexdegoutlier=exp_data_nonsexdeg%>% dplyr::filter(variant_cat == vcat)%>%dplyr::filter(abs(z)>=zscore)%>%nrow()
      exp_data_nonsexdeg_nonoutlier=exp_data_nonsexdeg%>% dplyr::filter(variant_cat == vcat)%>%dplyr::filter(abs(z)<zscore)%>%nrow()
      exptable<-rbind(c(exp_data_sexdeg_outlier,exp_data_sexdeg_nonoutlier),c(exp_data_nonsexdegoutlier,exp_data_nonsexdeg_nonoutlier))
      if(exp_data_sexdeg_outlier == 0 & exp_data_sexdeg_nonoutlier==0 & exp_data_nonsexdegoutlier==0 & exp_data_nonsexdeg_nonoutlier==0){
        risks = rbind(risks,data.frame(Risk = NA, 
                                       Lower = NA, 
                                       Upper = NA,
                                       Pval = NA,
                                       Cat=vcat,
                                       sexdeg=this_sexdeg,
                                       exp_nn=exp_data_nonsexdeg_nonoutlier,
                                       exp_ny=exp_data_nonsexdegoutlier,
                                       exp_yn=exp_data_sexdeg_nonoutlier,
                                       exp_yy=exp_data_sexdeg_outlier,
                                       cadd=my_cadd,chr=my_chr,
                                       tissue=my_tissue,
                                       num_outliers=exp_data_nonsexdegoutlier+exp_data_sexdeg_outlier,
                                       sex=sex,z=zscore,nphen=nphen,
                                       Type = 'RR_ALLOutliers_noRV'))
      } else{
        err = epitab(exptable, method = 'riskratio')
        risks = rbind(risks,data.frame(Risk = err$tab[2,5],
                                       Lower = err$tab[2,6],
                                       Upper = err$tab[2,7],
                                       Pval = err$tab[2,8],
                                       Cat=vcat,
                                       sexdeg=this_sexdeg,
                                       exp_nn=exp_data_nonsexdeg_nonoutlier,
                                       exp_ny=exp_data_nonsexdegoutlier,
                                       exp_yn=exp_data_sexdeg_nonoutlier,
                                       exp_yy=exp_data_sexdeg_outlier,
                                       cadd=my_cadd,chr=my_chr,
                                       tissue=my_tissue,
                                       num_outliers=exp_data_nonsexdegoutlier+exp_data_sexdeg_outlier,
                                       sex=sex,z=zscore,nphen=nphen,
                                       Type = 'RR_ALLOutliers_noRV')
         )
      }
    }
  exp_data_sexdeg_outlier=exp_data_sexdeg%>%dplyr::filter(abs(z)>=zscore)%>%nrow()
  exp_data_sexdeg_nonoutlier=exp_data_sexdeg%>%dplyr::filter(abs(z)<zscore)%>%nrow()
  exp_data_nonsexdegoutlier=exp_data_nonsexdeg%>%dplyr::filter(abs(z)>=zscore)%>%nrow()
  exp_data_nonsexdeg_nonoutlier=exp_data_nonsexdeg%>%dplyr::filter(abs(z)<zscore)%>%nrow()
  exptable<-rbind(c(exp_data_sexdeg_outlier,exp_data_sexdeg_nonoutlier),c(exp_data_nonsexdegoutlier,exp_data_nonsexdeg_nonoutlier))
  if(exp_data_sexdeg_outlier == 0 & exp_data_sexdeg_nonoutlier==0 & exp_data_nonsexdegoutlier==0 & exp_data_nonsexdeg_nonoutlier==0){
    risks = rbind(risks,data.frame(Risk = NA, #err$tab[2,5],
                                   Lower = NA, #err$tab[2,6],
                                   Upper = NA, #err$tab[2,7],
                                   Pval = NA, #err$tab[2,8],
                                   Cat="all",
                                   sexdeg=this_sexdeg,
                                   exp_nn=exp_data_nonsexdeg_nonoutlier,
                                   exp_ny=exp_data_nonsexdegoutlier,
                                   exp_yn=exp_data_sexdeg_nonoutlier,
                                   exp_yy=exp_data_sexdeg_outlier,
                                   cadd=my_cadd,chr=my_chr,
                                   tissue=my_tissue,
                                   num_outliers=exp_data_nonsexdegoutlier+exp_data_sexdeg_outlier,
                                   sex=sex,z=zscore,nphen=nphen,
                                   Type = 'RR_ALLOutliers_noRV')
    )
  } else {
    err = epitab(exptable, method = 'riskratio')
    risks = rbind(risks,data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   Cat="all",
                                   sexdeg=this_sexdeg,
                                   exp_nn=exp_data_nonsexdeg_nonoutlier,
                                   exp_ny=exp_data_nonsexdegoutlier,
                                   exp_yn=exp_data_sexdeg_nonoutlier,
                                   exp_yy=exp_data_sexdeg_outlier,
                                   cadd=my_cadd,chr=my_chr,
                                   tissue=my_tissue,
                                   num_outliers=exp_data_nonsexdegoutlier+exp_data_sexdeg_outlier,
                                   sex=sex,z=zscore,nphen=nphen,
                                   Type = 'RR_ALLOutliers_noRV')
                  )
  }
  
}


#RV ONLY
for (this_sexdeg in sexdeg_types){
  if(this_sexdeg=="all"){
    print("all sexDEGs")
    sexdeg<-exp_data %>% dplyr::filter(abs(beta)>=beta_min)
    non_sexdeg<-exp_data%>% dplyr::filter(abs(beta)<beta_min)
  }else if(this_sexdeg=="male_biased"){
    print("male sexDEGs")
    sexdeg<-exp_data %>% dplyr::filter(beta<=(-1*beta_min))
    non_sexdeg<-exp_data%>% dplyr::filter(beta>(-1*beta_min))
  }else if(this_sexdeg=="female_biased"){
    print("female sexDEGs")
    sexdeg<-exp_data%>% dplyr::filter(beta>=beta_min)
    non_sexdeg<-exp_data%>% dplyr::filter(beta<beta_min)
  }else{stop("invalid sexdeg type")}
  for (vcat in (vcats)) {
    print(vcat)

    exp_nn = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare") )
    exp_ny = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))
    exp_yn = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))
    exp_yy = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare"))
    
    exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
    print(exptable)
    if(exp_nn == 0 & exp_ny==0 & exp_yn==0 & exp_yy==0){
      risks = rbind(risks, data.frame(Risk = NA , # err$tab[2,5],
                                      Lower = NA , # err$tab[2,6],
                                      Upper = NA , # err$tab[2,7],
                                      Pval =NA ,# err$tab[2,8],
                                      Cat = vcat,
                                      sexdeg=this_sexdeg,
                                      exp_nn=exp_nn,
                                      exp_ny=exp_ny,
                                      exp_yn=exp_yn,
                                      exp_yy=exp_yy,
                                      cadd=my_cadd,chr=my_chr,
                                      tissue=my_tissue,
                                      num_outliers=nrow(sexdeg),
                                      sex=sex,z=zscore,nphen=nphen,
                                      Type = 'RR_RV_only'))   
      
    }else {  
      err = epitab(exptable, method = 'riskratio')
      risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                      Lower = err$tab[2,6],
                                      Upper = err$tab[2,7],
                                      Pval = err$tab[2,8],
                                      Cat = vcat,
                                      sexdeg=this_sexdeg,
                                      exp_nn=exp_nn,
                                      exp_ny=exp_ny,
                                      exp_yn=exp_yn,
                                      exp_yy=exp_yy,
                                      cadd=my_cadd,chr=my_chr,
                                      tissue=my_tissue,
                                      num_outliers=nrow(sexdeg),
                                      sex=sex,z=zscore,nphen=nphen,
                                      Type = 'RR_RV_only'))   
    }
  }
  print("all")
  # exp_nn_all=nrow(dplyr::filter(sexdeg,has_variant == "rare"))
  # exp_ny_all=nrow(dplyr::filter(sexdeg,has_variant != "rare"))
  # exp_yn_all=nrow(dplyr::filter(non_sexdeg,has_variant == "rare"))
  # exp_yy_all=nrow(dplyr::filter(non_sexdeg,has_variant != "rare"))
  exp_nn_all = nrow(non_sexdeg%>% dplyr::filter( has_variant != "rare") )
  exp_ny_all = nrow(non_sexdeg  %>% dplyr::filter(has_variant == "rare"))
  exp_yn_all = nrow(sexdeg %>% dplyr::filter( has_variant != "rare"))
  exp_yy_all = nrow(sexdeg %>% dplyr::filter( has_variant == "rare"))
  exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
  print(exptable_all)
  if(exp_nn_all == 0 & exp_ny_all==0 & exp_yn_all==0 & exp_yy_all==0){
    risks = rbind(risks, data.frame(Risk = NA, #err_all$tab[2,5],
                                    Lower =  NA, #err_all$tab[2,6],
                                    Upper =  NA, # err_all$tab[2,7],
                                    Pval =  NA, #err_all$tab[2,8],
                                    Cat = "all",
                                    sexdeg=this_sexdeg,
                                    exp_nn=exp_nn_all,
                                    exp_ny=exp_ny_all,
                                    exp_yn=exp_yn_all,
                                    exp_yy=exp_yy_all,
                                    cadd=my_cadd,chr=my_chr,
                                    tissue=my_tissue,
                                    num_outliers=nrow((sexdeg)),
                                    sex=sex,z=zscore,nphen=nphen,
                                    Type = 'RR_RV_only'))   
    
  }else{
    err_all = epitab(exptable_all, method = 'riskratio')
    risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                    Lower = err_all$tab[2,6],
                                    Upper = err_all$tab[2,7],
                                    Pval = err_all$tab[2,8],
                                    Cat = "all",
                                    sexdeg=this_sexdeg,
                                    exp_nn=exp_nn_all,
                                    exp_ny=exp_ny_all,
                                    exp_yn=exp_yn_all,
                                    exp_yy=exp_yy_all,
                                    cadd=my_cadd,chr=my_chr,
                                    tissue=my_tissue,
                                    num_outliers=nrow((sexdeg)),
                                    sex=sex,z=zscore,nphen=nphen,
                                    Type = 'RR_RV_only'))   
  } 
  risks = risks %>% arrange(by=Risk) 
  risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))


}




risks$padj<-p.adjust(risks$Pval,method="BH")

write.csv(risks,  file=outfile,quote=F,sep="\t",row.names = F)
print(paste0("COMPLETE: ",outfile))
