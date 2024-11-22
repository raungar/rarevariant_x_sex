#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


# option_list = list(
#                 make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
#                 make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
#                 make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
#                 make_option(c("--nphen"), type = 'numeric', default = NULL, help = "min tissues"),
#                 make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),
#                 make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
#                 make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF"),
#                 make_option(c("--beta_min"), type = 'numeric', default = NULL, help = "beta_min"),
#                 make_option(c("--tissue"), type = 'character', default = NULL, help = "tissue")
# )
# 
# 
# opt_parser = OptionParser(option_list = option_list)
# opt = parse_args(opt_parser)
# infile <- as.character(opt$infile)
# outfile <- as.character(opt$outfile)
# tissue <- as.character(opt$tissue)
# zscore <- as.numeric(opt$zscore)
# nphen <- as.numeric(opt$nphen)
# sex <- as.character(opt$sex)
# min_maf<-as.numeric(opt$min_maf)
# max_maf<-as.numeric(opt$max_maf)
# beta_min<-as.numeric(opt$beta_min)

#infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Combined/aut_beta0.111_SKINS_z3_nphen5_f_linc_prot.txt.gz"
# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/aut_beta0.111_NERVET_m_linc_prot.txt.gz"
# #infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen4_noglobal_medz_varAnnot_x_m.txt"
# zscore<-3
min_maf<-0
max_maf<-0.01
zscore=NA
nphen=1
# beta_min<-0.111
outfile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/all_risks_RVONLY.txt"
### expression



   
   
 ##########CHANGE THIS LINE FOR WHAT U WANT THE RR TO BE....



print("relative risk  all sexdegs")
### Relative risk
##rr=Noutliers-sexDEG-rv/NsexDEG-outlier
#    -------------------------
#     N-nonoutliers-nonsexDEG-rv/NsexDEG-nonoutlier
beta_min=0.111
tissue_list<-unique(c("ADPSBQ","ADPSBQ","ADPVSC","ADPVSC","ADRNLG","ADRNLG","ARTAORT","ARTAORT","ARTCRN","ARTCRN","ARTTBL","ARTTBL","BREAST","BREAST","BRNACC","BRNACC","BRNAMY","BRNAMY","BRNCDT","BRNCDT","BRNCHA","BRNCHA","BRNCHB","BRNCHB","BRNCTXA","BRNCTXA","BRNCTXB","BRNCTXB","BRNHPP","BRNHPP","BRNHPT","BRNHPT","BRNNCC","BRNNCC","BRNPTM","BRNPTM","BRNSNG","BRNSNG","BRNSPC","BRNSPC","CLNSGM","CLNSGM","CLNTRN","CLNTRN","ESPGEJ","ESPGEJ","ESPMCS","ESPMCS","ESPMSL","ESPMSL","FIBRBLS","FIBRBLS","HRTAA","HRTAA","HRTLV","HRTLV","KDNCTX","KDNCTX","LCL","LCL","LIVER","LIVER","LUNG","LUNG","MSCLSK","MSCLSK","NERVET","NERVET","PNCREAS","PNCREAS","PTTARY","PTTARY","SKINNS","SKINNS","SKINS","SKINS","SLVRYG","SLVRYG","SNTTRM","SNTTRM","SPLEEN","SPLEEN","STMACH","STMACH","THYROID","THYROID","WHLBLD","WHLBLD"))
sex_list<-c("m","f")
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Cat = character(),
                   sexdeg=character(),tissue=character(),num_outliers=character(),sex=character(),z=numeric(),nphen=numeric(),Type=character())
for(tissue in tissue_list){
  for(sex in sex_list){
    print('Reading expression')
    infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/aut_beta0.111_NERVET_m_linc_prot.txt.gz"
    infile<-paste0("/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/aut_beta0.111_",
                   tissue,"_",sex,"_linc_prot.txt.gz")
    print(paste0("reading in: ",infile))
    exp_data = fread(infile,data.table=F)
    colnames(exp_data)<-c("ind","ensg","z","chr.x","beta", "chr.y","pos","vartype","sex",
                          "ref","alt","bothvars", "gtex_maf","gnomad_maf_both","gnomad_maf_m","gnomad_maf_f","genetype","mafdiff","numrv")
    vcats = na.omit(unique(exp_data$variant_cat))
    
    
    print("filter for MAF")
    #choose rare/common
    has_variant<-apply((exp_data),1,function(x){
      this_gnomad_both_maf<-as.numeric(x[which(colnames(exp_data)=="gnomad_maf_both")])
      #no variants found w/in 10kb of gene, so not rare
      if(is.na(this_gnomad_both_maf)){"common"}
      else if (this_gnomad_both_maf>=min_maf & this_gnomad_both_maf<max_maf){"rare"}
      else if (this_gnomad_both_maf>=max_maf){"common"}
      else (stop("ERROR: VARIANT NOT MAKING SENSE"))
    })
    print("filter for z score")
    ### get relative risk per category
    exp_data$has_variant<-has_variant
    exp_data$variant_cat<-exp_data$genetype
    
    sexdeg_types<-c("all","male_biased","female_biased")
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
          exp_nn = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant == "rare"))
          exp_ny = nrow(sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))
          exp_yn = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))
          exp_yy = nrow(non_sexdeg %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare") )
          exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
          print(exptable)
          if(exp_nn == 0 & exp_ny==0 & exp_yn==0 & exp_yy==0){
            risks = rbind(risks, data.frame(Risk = NA , # err$tab[2,5],
                                            Lower = NA , # err$tab[2,6],
                                            Upper = NA , # err$tab[2,7],
                                            Pval =NA ,# err$tab[2,8],
                                            Cat = vcat,
                                            sexdeg=this_sexdeg,
                                            tissue=tissue,
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
                                            tissue=tissue,
                                            num_outliers=nrow(sexdeg),
                                            sex=sex,z=zscore,nphen=nphen,
                                            Type = 'RR_RV_only'))   
          }
        }
        print("all")
        exp_nn_all=nrow(dplyr::filter(sexdeg,has_variant == "rare"))
        exp_ny_all=nrow(dplyr::filter(sexdeg,has_variant != "rare"))
        exp_yn_all=nrow(dplyr::filter(non_sexdeg,has_variant == "rare"))
        exp_yy_all=nrow(dplyr::filter(non_sexdeg,has_variant != "rare"))
        exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
        print(exptable_all)
        if(exp_nn_all == 0 & exp_ny_all==0 & exp_yn_all==0 & exp_yy_all==0){
          risks = rbind(risks, data.frame(Risk = NA, #err_all$tab[2,5],
                                          Lower =  NA, #err_all$tab[2,6],
                                          Upper =  NA, # err_all$tab[2,7],
                                          Pval =  NA, #err_all$tab[2,8],
                                          Cat = "all",
                                          sexdeg=this_sexdeg,
                                          tissue=tissue,
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
                                          tissue=tissue,
                                          num_outliers=nrow((sexdeg)),
                                          sex=sex,z=zscore,nphen=nphen,
                                          Type = 'RR_RV_only'))   
        } 
        risks = risks %>% arrange(by=Risk) 
        risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
      }
  }
}
risks$padj<-p.adjust(risks$Pval,method="BH")

write.csv(risks,  file=outfile,quote=F,sep="\t")
print(paste0("COMPLETE: ",outfile))
