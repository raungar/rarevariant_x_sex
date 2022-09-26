library("ggplot2")
library(tidyverse)
library("forcats")
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
# mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8redo/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_all/RR"

#groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress")
groups=c("m","f","both") #,"allboth") #, "both_half","both_half.sex","both_half.regress
# groups=c("allboth") #, "both_half","both_half.sex","both_half.regress
my_cat=c("xci","par","strata","","pq")
#my_cat=c(".xci",".par",".strata")
my_cat=""
# chr_types=c(as.character(c(1:5,7:11,15:21)),"x") #,
# chr_types=c(as.character(1:22),"x") #,
chr_types=c(as.character(2:22),"x"); filter_version=c("typesALL")
 #,
chr_types=c("x","AllAut","7sub") #,"AllAut","7") #,
# chr_types=c("x","7") #,
chr_types=c("x") #aut
nphen=c(3) #,5)
z=c(2,2.5,3)
z=2.5
CADD<-c(0,15)
CADD=0
risk_cat<-c("relative_risk") #,"absolute_risk")
#       maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
# maxtomin_dic<-c("0.05","0.01","0.001","0.0001","0")
       maxmaf<-c("0.1", "0.05","0.01", "0.001","0.0001") #,"0.001","0.0001")
 maxtomin_dic<-c("0.05","0.01","0.001","0"     ,"0")
# maxtomin_dic<-c("0","0","0","0","0","0")
names(maxtomin_dic)<-maxmaf
filter_version=c("typesSeenTwice","typesALL","typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                 "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                 "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesGQ5BlacklistRemovedALL")
outlier_types=c("outliers");nphen=3; #"outliersTOP",
outlier_types=c("outliersAdipose-Subcutaneous","outliersAdipose-Visceral-Omentum","outliersAdrenal-Gland","outliersArtery-Aorta","outliersArtery-Coronary","outliersArtery-Tibial","outliersBrain-Caudate-basal-ganglia",
                             "outliersBrain-Cerebellar-Hemisphere","outliersBrain-Cerebellum","outliersBrain-Cortex","outliersBrain-Nucleus-accumbens-basal-ganglia","outliersBreast-Mammary-Tissue",
                             "outliersCells-Cultured-fibroblasts","outliersCells-EBV-transformed-lymphocytes","outliersColon-Sigmoid","outliersColon-Transverse","outliersEsophagus-Gastroesophageal-Junction",
                             "outliersEsophagus-Mucosa","outliersEsophagus-Muscularis","outliersHeart-Atrial-Appendage","outliersHeart-Left-Ventricle","outliersLiver","outliersLung","outliersMuscle-Skeletal","outliersNerve-Tibial","outliersOvary",
                             "outliersPancreas","outliersPituitary","outliersSkin-Not-Sun-Exposed-Suprapubic","outliersSkin-Sun-Exposed-Lower-leg","outliersSmall-Intestine-Terminal-Ileum","outliersSpleen","outliersStomach","outliersThyroid",
                             "outliersUterus","outliersVagina","outliersWhole-Blood");nphen=1;
outlier_types_female<-c("outliersOvary","outliersUterus","outliersVagina")
#collapsetypes=c("collapsed","collapsedbyVARLOC") #,"collapsedbyTSS")
# outlier_types=outlier_types_female
collapsetypes=c("collapsedexon","collapsedgenestart","collapsedgeneend","collapsedintron","collapsed")
# collapsetypes=c("collapsedexon","collapsedgenestart","collapsedgeneend","collapsedintron","collapsed",
collapsetypes=c("collapsedexon","collapsedgenestart","collapsedgeneend","collapsedintron","collapsed",
                "collapsed3-prime-UTR-variant","collapsed5-prime-UTR-variant","collapsedcoding-sequence-variant",
               "collapseddownstream-gene-variant","collapsedincomplete-terminal-codon-variant","collapsedintergenic-variant",
               "collapsedintron-variant","collapsedmissense-variant","collapsedNMD-transcript-variant","collapsednon-coding-transcript-exon-variant",
               "collapsednon-coding-transcript-variant","collapsedregulatory-region-variant","collapsedsplice-acceptor-variant",
               "collapsedsplice-donor-variant","collapsedsplice-region-variant","collapsedstart-lost","collapsedstop-gained","collapsedstop-lost",
               "collapsedsynonymous-variant","collapsedTF-binding-site-variant","collapsedupstream-gene-variant")
max_outliers=3
windows<-c("100","500","1000","5000","10000")
 windows<-c("500","5000","10000")
 windows=5000
collapsetypes="collapsed"
# windows=c(1000)
# maxmaf<-c("0.1","0.001","0.0001") #,"0.001","0.0001")
# maxtomin_dic<-c("0","0","0")
# names(maxtomin_dic)<-maxmaf


all_risks<-data.frame(matrix(nrow=0,ncol=10))
colnames(all_risks)<-c("Risk","Lower","Upper" ,"Pval","Cat","Type" , "z","nphen","chr","sex")
for(this_cat in my_cat){
  for(this_collapsetype in collapsetypes){
  for(this_outlier in outlier_types){
    for(this_maxmaf in maxmaf){
      print(this_maxmaf)
  for (this_group in groups){
    for(this_chr in chr_types){
      for(this_nphen in nphen){
       for(this_z in z){
         for(this_filt in filter_version){
           for(cadd_min in CADD){
             for(this_window in windows){
          # if(this_chr=="x"){
          #      file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,
          #                  "_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_",this_group,"_CADD","typesGQ5BlacklistRemovedALL",
          #                  "_CADD",cadd_min,"_linc_prot.",this_cat,".csv")
          #      #relative_risk_x_outliers_z2.5_nphen5_x_both_min0max0.001_CADDtypesGQ5BlacklistRemovedALL_CADD15_linc_prot.csv
          #      this_filt="GQ5BlacklistRemovedALL"
          #    }else{
          #    file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,"_",this_group,"_",this_chr,
          #                "_min0max",this_maxmaf,"_CADD",this_filt,"_CADD",cadd_min,"_linc_prot.",this_cat,".txt")
          #    }
             if(this_chr=="x"){
               file=paste0(mydir,"/relative_risk_x_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                           "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
                           cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
               this_filt="GQ5BlacklistRemovedALL"
               #Output/enrichments_v8eqtl/RR/relative_risk_x_collapsed_outliers_z2.5_nphen3_aut_f_min0max0.001_CADDtypesALL_CADD0_linc_prot_maxoutliers3.csv
               #_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
               #"-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
               #this_isnull,"-",use_chrgroup,"-binary.txt.gz")
             }else if(grepl("sub",this_chr)){
               file=paste0(mydir,"/relative_risk_",this_chr,"_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                           "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",
                           cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
             } else{
               file=paste0(mydir,"/relative_risk_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                           "_",this_group,"_",this_chr,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
                           "_linc_prot_maxoutliers",max_outliers,"_window",this_window,".txt")
             }
            # risks=read.csv(file[1])
#
             risks=tryCatch({read.csv(file[1])},
                            error=function(err){return("ERROR FILE NOT FOUND")},
                            warnings=function(war){return("ERROR FILE NOT FOUND")})
             if(risks=="ERROR FILE NOT FOUND"){
               print(paste0(file[1], " WAS NOT FOUND"))
               next
             }
           #load(file)
          risks$z<-this_z
          risks$nphen<-this_nphen
          risks$chr<-this_chr
          risks$sex<-this_group
          risks$maxmaf<-this_maxmaf
          risks$filt<-this_filt
          risks$cadd<-cadd_min
          risks$collapse_method<-this_collapsetype
          risks$thiscat<-this_cat
          risks$outlierType<-this_outlier
          risks$gene_window<-this_window
          risks$real_or_null<-"real"
          risks$discrete_or_continuous<-"discrete"
          colnames(risks)[6]<-"CATEGORY"
          #assign(varname,risks)
          #print(varname)
  
          all_risks<-rbind(all_risks,risks)
           }
          # break
         }
         }
       }
        }
      }
    }
  }
  }
}
}
all_risks$outliers_tested<-all_risks$exp_yn+all_risks$exp_yy
all_risks$nonoutliers_tested<-all_risks$exp_nn+all_risks$exp_ny
all_risks$prop_outliers<-all_risks$outliers_tested/all_risks$nonoutliers_tested
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/ following this
risks_to_compare<-all_risks %>% group_by(outlierType,CATEGORY,z,nphen,cadd,chr,sex,maxmaf,real_or_null,
                                         exp_type,gene_window,var_location,veptype,collapse_method) %>% 
  dplyr::summarise(across(c(Risk,Lower,Upper),log)) %>%
  mutate(SE=abs(Lower-Upper)/(2*1.96))%>%dplyr::filter(sex=="m" | sex=="f")
risks_to_compare_side<-pivot_wider(risks_to_compare,names_from=c(sex),
                                   id_cols=c(outlierType,CATEGORY,z,nphen,cadd,chr,maxmaf,real_or_null,exp_type,
                                             gene_window,var_location,veptype,collapse_method) ,
                                   values_from=c(Risk,SE)) #id_cols=, names_from=,
risks_to_compare_side$RiskDiff<-risks_to_compare_side$Risk_f-risks_to_compare_side$Risk_m
risks_to_compare_side$SEDiff<-sqrt(risks_to_compare_side$SE_f**2+risks_to_compare_side$SE_m**2)
risks_to_compare_side$zDiff<-(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
risks_to_compare_side$p<-pnorm(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
#n is 2 for cadd * 3 for z * 3 for nphen * 3 for MAF * 24 for chr  (23 + allAut)=1296 OR 54* 2 for outliers
risks_to_compare_side$padj<-p.adjust(risks_to_compare_side$p,method = "BH") #,n=108)
all_risks_p<-merge(all_risks,risks_to_compare_side[,c("outlierType","exp_type","CATEGORY","z","nphen","cadd","chr",
                                                      "maxmaf","p","padj","gene_window","var_location","veptype","collapse_method")])
all_risks_p$gene_window <- factor(all_risks_p$gene_window, levels = c("100","500","1000","5000","10000"))
all_risks_p$outlierType <-str_replace(all_risks_p$outlierType,"outliers","")
# regresstype="incl_sex"
# null_or_not="preprocessing_v8"
# chrgroups=c("x_or_aut")
# groups=c("m","f","both_half")
# chr_types=c("x","7","aut") #,"AllAut","7") #,
# 
# #relative_risk_min0max0.0001_x_outliersTOP_cadd15_noglobal_medz-zthresh3-nphen3-m-protect_sex-preprocessing_v8-x-continuous_alpha0.5.txt.gz
# for(this_isnull in null_or_not){
#   for(this_outlier in outlier_types){
#     for(this_maxmaf in maxmaf){
#       for (this_sex in groups){
#         for (this_chrgroup in chrgroups){
#           for(this_chr in chr_types){
#           for(this_nphen in nphen){
#             for(this_z in z){
#               for(this_regresstype in regresstype){
#                 for(cadd_min in CADD){
#                     if(this_chr=="x" & this_chrgroup=="x_or_aut"){
#                       use_chrgroup="x"
#                     }
#                   else if(this_chr!="x" & this_chrgroup=="x_or_aut"){
#                     use_chrgroup="aut"
#                   } 
#                   else{
#                     use_chrgroup="both"
#                     }
#                   # file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
#                   #             "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
#                   #             this_isnull,"-",use_chrgroup,"-continuous_alpha0.5.txt.gz")
#                   # file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
#                   #             "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
#                   #             this_isnull,"-",use_chrgroup,"-binary.txt.gz")
#                   if(this_chr=="x"){
#                     file=paste0(mydir,"/relative_risk_x_",this_outlier,"_z",this_z,"_nphen",this_nphen,
#                                 "_x_",this_sex,"_min0max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
#                                 cadd_min,"_linc_prot.csv")
#                                 #_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
#                                 #"-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
#                                 #this_isnull,"-",use_chrgroup,"-binary.txt.gz")
#                   }else{
#                     file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,
#                                 "_",this_sex,"_",this_chr,"_min0max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
#                                 "_linc_prot.txt")
#                   }
# 
#                    risks=read.csv(file[1])
#                   #load(file)
#                   risks$z<-this_z
#                   risks$nphen<-this_nphen
#                   risks$chr<-this_chr
#                   if(this_sex=="both_half"){use_sex="both"}else{use_sex=this_sex}
#                   risks$sex<-use_sex
#                   risks$maxmaf<-this_maxmaf
#                   risks$filt<-this_filt
#                   risks$cadd<-cadd_min
#                   risks$thiscat<-this_cat
#                   risks$outlierType<-this_outlier
#                   colnames(risks)[6]<-"CATEGORY"
#                   risks$real_or_null<-"null"
#                   risks$discrete_or_continuous<-"continuous"
#                   #assign(varname,risks)
#                   #print(varname)
#                   # break
#                   
#                   all_risks<-rbind(all_risks,risks)
# }}}}}}}}}}



# for (this_group in groups){
#   for(this_type in types){
#     for(this_risk in risk_cat){
#       file=paste0(mydir,"/",this_risk,"_",
#                        this_type,"_",this_group,".RData")
#       load(file)
#       print(file)
#       if(exists("all_coefs")){
#         varname<-paste0("contrisk_z",sapply(strsplit(sapply(strsplit(file,"continuous_risk_"),"[[",2),".RData"),"[[",1))
#         assign(varname,all_coefs)
#         rm(all_coefs)      }
#       if(exists("risks")){
#         varname<-paste0("relrisk_",sapply(strsplit(sapply(strsplit(file,"relative_risk_"),"[[",2),".RData"),"[[",1))
#         print(varname)
#         assign(varname,risks)
#         rm(risks)
#       }
#     }
#   }
# }


 # to_plot<-all_risks%>%dplyr::filter(chr=="x")  %>% dplyr::filter(z==2.5) %>% dplyr::filter(exp_type=="all") %>% 
 #   dplyr::filter(maxmaf=="0.01") %>% dplyr::filter(gene_window==10000)
# 
#%>%dplyr::filter(chr=="1")# %>%dplyr::filter(maxmaf=="0.0001")

 
  # filter(real_or_null=="real")%>%filter(cadd=="15")%>% mutate(chr=fct_relevel(chr,"x","7","AllAut"))
# to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers")%>%
#             filter(maxmaf=="0.001")%>%filter(cadd=="15") #%>% mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% filter(real_or_null=="real")
#%>%filter(chr=="AllAut")
 ###x
#,alpha=discrete_or_continuous$
#
to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="") %>% dplyr::filter(exp_type=="all") %>%dplyr::filter(var_location=="all") %>%
  dplyr::filter(z==2.5)%>%dplyr::filter(chr=="x")%>%filter(collapse_method=="collapsed") %>%
  mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%dplyr::filter(veptype=="all") %>% dplyr::filter(gene_window==5000) 

 ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),group=sex,fill=sex)) + 
   geom_crossbar(aes(y=Risk,ymin=Lower,ymax=Upper), position=position_dodge(width=1)) +
   theme_bw(base_size=8)+
   xlab("chr")+ylab("Relative Risk")+labs(fill="Group")+
   ggtitle(paste0("Relative Risk: CADD=,",unique(to_plot$cadd), #(outlierType=",unique(to_plot$outlierType),
                  ",nphen=",unique(to_plot$nphen),",window=",unique(to_plot$gene_window),",MAF=",unique(to_plot$maxmaf),"), num outlier range [",
                  min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
   # scale_fill_manual(values=c("#FCF5A9","#97D6F2"))+
  ### scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
   # scale_fill_manual(values=c("#FCF5A9","#97D6F2"))+
   # scale_color_manual(values=c("#dbab3b","#5d8596"))+
   # scale_alpha_discrete(range=c(0.3,1))+
   guides(colour=FALSE)+
   # facet_wrap(~paste0("z_type=",exp_type)*paste0("z=",z)*paste0("sex=",sex),ncol =3) + #,scales="free_y"
    facet_wrap(~paste0("z=",z)*paste0("cadd=",cadd),ncol =3,scales="free_y") + #,scales="free_y"
   
   # facet_wrap(~paste0("maf=",maxmaf)*paste0("cadd=",cadd)*paste0("nphen=",nphen),ncol =2,scales="free_y")+
   # facet_wrap(~,scales="free_y")+
   geom_text(aes(y=Risk,#group=interaction(real_or_null,sex)
                 label=paste0("n=",num_outliers,"\nRR=",round(Risk,1),"\np=",round(Pval,2)),
                 "\n[",exp_nn,"|",exp_ny,"]","\n[",exp_yn,"|",exp_yy,"]"), 
             position = position_dodge(width = 1), size=2.5,vjust=-1.25) +
   theme(axis.text.x = element_text(angle = 45,  hjust=1))+
   
   geom_hline(yintercept=1,color="red")
 
 
 ###line for mafs
 to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>%dplyr::filter(var_location=="all") %>%# %>% dplyr::filter(exp_type=="all")
   dplyr::filter(z==2.5)%>%filter(collapse_method=="collapsed") %>%dplyr::filter(num_outliers>1)%>% #dplyr::filter(sex=="allboth")%>%  #%>%dplyr::filter(chr=="x")
   mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%dplyr::filter(cadd==15) %>% dplyr::filter(gene_window==5000)  %>% dplyr::filter(exp_type=="all")%>% 
   dplyr::filter(veptype=="all")%>%mutate(is_x=ifelse(chr=="x","x","aut")) # %>% dplyr::filter(maxmaf>0.001)
 #%>%filter(maxmaf==0.001) # %>% filter(chr=="x")  %>% dplyr::filter(nphen==3)%>% dplyr::filter(sex=="both")
 ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=interaction(chr,sex),color=sex,shape=is_x,
                  label=num_outliers)) + #  shape=sex,
      theme(axis.text.x = element_text(angle = 45,  hjust=1))+
   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
   geom_line(position=position_dodge(width=0.5))+
   theme_bw(base_size=15)+
   geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
   xlab("MAF")+
   ylab("Relative Risk")+ #xlim(c(0,0.1))+
   labs(fill="Group")+
   ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),",CADD=",unique(to_plot$cadd),",sex=",unique(to_plot$sex),
                  ",nphen=",unique(to_plot$nphen),",window=",unique(to_plot$gene_window),",z=",unique(to_plot$z),"), num outlier range [",
                  min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   #  scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
   # scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
   
   # scale_color_manual(values=c("#926fa8","#47265c","#99176e"))+
   #guides(colour=FALSE)+ #paste0("z_type=",var_location)*paste0("z=",z)*
   facet_wrap(~paste0("genewindow=",(gene_window))*paste0(exp_type, " eOutlier"),
           scales="free_y") + #,scales="free_y" ncol=2
   geom_hline(yintercept=1,color="red",linetype="dashed") +
   # geom_text(aes(y=Risk,#group=interaction(real_or_null,sex)
   #               label=paste0("n=",num_outliers,"\nRR=",round(Risk,1),"\np=",round(Pval,2),
   #                            "\n[",exp_nn,"|",exp_ny,"]","\n[",exp_yn,"|",exp_yy,"]")), 
   #           position = position_dodge(width = 0.5), size=2.5,vjust=-1.25) +
   theme(axis.text.x = element_text(angle = 45,  hjust=1))
 
 
 ### fancy version
 to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>%dplyr::filter(var_location=="all") %>% dplyr::filter(chr=="7sub")%>%
   dplyr::filter(z==2.5)%>%filter(collapse_method=="collapsed") %>%dplyr::filter(num_outliers>1)%>%# dplyr::filter(sex!="both")%>%
   mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%dplyr::filter(cadd==15) %>% dplyr::filter(gene_window==5000)  %>% dplyr::filter(exp_type=="all")%>% 
   dplyr::filter(veptype=="all")%>%mutate(is_x=ifelse(chr=="x","x","aut"))# %>%dplyr::filter(chr=="7sub")
 ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=interaction(chr,sex),color=sex,shape=chr,
                    label=num_outliers)) + #  shape=sex,
   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
   geom_line(position=position_dodge(width=0.5))+
   theme_classic(base_size=15)+
   geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
   xlab("MAF")+
   ylab("Relative Risk")+ #xlim(c(0,0.1))+
   #theme(legend.position="none")+
   ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),",CADD=",unique(to_plot$cadd),",sex=",unique(to_plot$sex),
                  ",nphen=",unique(to_plot$nphen),",window=",unique(to_plot$gene_window),",z=",unique(to_plot$z),"), num outlier range [",
                  min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   #  scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
    scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+ #sex
   # scale_color_manual(values=c("#7E6B8F","#07004D","#1B998B"))+ #chrs
   # theme(axis.text.x = element_text(angle = 45,  hjust=1))
   geom_hline(yintercept=1,color="red",linetype="dashed")
 
 
 
 
 ###switch so by group
 to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>%dplyr::filter(var_location=="all") %>%# %>% dplyr::filter(exp_type=="all")
   dplyr::filter(z==2.5)%>%filter(collapse_method=="collapsed") %>%dplyr::filter(num_outliers>1)%>%  #%>%dplyr::filter(chr=="x")
   mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%dplyr::filter(cadd==15) %>% dplyr::filter(gene_window==5000)  %>% dplyr::filter(exp_type=="all")%>% 
   dplyr::filter(veptype=="all")%>%mutate(is_x=ifelse(chr=="x","x","aut")) # %>% dplyr::filter(maxmaf>0.001)
 ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=interaction(exp_type,sex),color=exp_type,shape=is_x,
                    label=num_outliers)) + #  shape=sex,
   theme(axis.text.x = element_text(angle = 45,  hjust=1))+
   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
   geom_line(position=position_dodge(width=0.5))+
   theme_bw(base_size=15)+
   geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
   xlab("MAF")+ylab("Relative Risk")+ labs(fill="Group")+
   ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),",CADD=",unique(to_plot$cadd),",sex=",unique(to_plot$sex),
                  ",nphen=",unique(to_plot$nphen),",window=",unique(to_plot$gene_window),",z=",unique(to_plot$z),"), num outlier range [",
                  min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   facet_wrap(~paste0("genewindow=",(gene_window))*paste0("chr",chr),
              scales="free_y") + #,scales="free_y" ncol=2
   geom_hline(yintercept=1,color="red",linetype="dashed") +
   geom_text(aes(y=Risk,#group=interaction(real_or_null,sex)
                 label=paste0("n=",num_outliers,"\nRR=",round(Risk,1),"\np=",round(Pval,2),
                              "\n[",exp_nn,"|",exp_ny,"]","\n[",exp_yn,"|",exp_yy,"]")), 
             position = position_dodge(width = 0.5), size=2.5,vjust=-1.25) +
   theme(axis.text.x = element_text(angle = 45,  hjust=1))
   
 
   ###by window
   to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%
     filter(cadd==15)   %>% dplyr::filter(z==2.5)  %>% dplyr::filter(exp_type=="over") %>% # dplyr::filter(maxmaf==0.001) %>% 
    dplyr::filter(collapse_method=="collapsed") %>% dplyr::filter(var_location=="all")#%>%filter(maxmaf==0.001) # %>% filter(chr=="x")  %>% dplyr::filter(nphen==3)%>% dplyr::filter(sex=="both")
   ggplot(to_plot,aes(x=((gene_window)),y=Risk,group=paste0(maxtomin_dic[maxmaf],"-",maxmaf),color=paste0(maxtomin_dic[maxmaf],"-",maxmaf),#shape=sex,
                      label=num_outliers)) + 
     theme(axis.text.x = element_text(angle = 45,  hjust=1))+
     geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
     geom_line(position=position_dodge(width=0.5))+
     theme_bw(base_size=5)+
     geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
     xlab("MAF")+
     ylab("Relative Risk")+ #xlim(c(0,0.1))+
     labs(fill="Group")+
     ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),",CADD=",unique(to_plot$cadd),
                    ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                    min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
     # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
     # scale_color_manual(values=c("#926fa8","#47265c","#99176e"))+  
     #guides(colour=FALSE)+
     facet_wrap(~paste0("sex=",sex)*paste0("z=",z),ncol =3,scales="free_y") + #,scales="free_y"
     geom_hline(yintercept=1,color="red",linetype="dashed") #+
   geom_text(aes(x=interaction(maxmaf,chr),y=Risk,label=ifelse(Pval<0.01,
                                                               paste0("n=",num_outliers,"\nRR=",round(Risk,1)))),
             position = position_dodge(width = 1), vjust=-.25,size=2.5) +
     geom_hline(yintercept=1,color="red") #+#+
   
   
   sig_veptypes=all_risks%>%dplyr::filter(Pval<0.001) %>%pull(veptype) %>% unique()
   
   ###by tissue!!!#paste0(maxtomin_dic[maxmaf],"-",maxmaf)
   to_plot=all_risks%>%dplyr::filter(CATEGORY=="all") %>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%filter(var_location=="all")%>%dplyr::filter(veptype%in%sig_veptypes)%>%
     filter(cadd==15)   %>% dplyr::filter(z==2.5)  %>% dplyr::filter(exp_type=="all") %>%  dplyr::filter(maxmaf==0.001) %>% dplyr::filter(sex=="m")%>%
     dplyr::filter(collapse_method=="collapsed")%>%filter(gene_window==500) # %>% filter(chr=="x")  %>% dplyr::filter(nphen==3)%>% dplyr::filter(sex=="both")
   # ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=outlierType,color=outlierType,
   ggplot(to_plot,aes(x=veptype,y=Risk,group=outlierType,color=outlierType,
                      
                      label=num_outliers)) + 
     geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.8)) +
     # geom_line(position=position_dodge(width=0.5))+
     theme_bw(base_size=7)+
     # geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
      geom_point(aes(size=is_sig),position=position_dodge(width=0.8))+
     xlab("MAF")+
     ylab("Relative Risk")+ #xlim(c(0,0.1))+
     labs(fill="Group")+
     theme(axis.text.x = element_text(angle = 45,  hjust=1))+
     
     ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),",CADD=",unique(to_plot$cadd),
                    ", window: ",unique(to_plot$gene_window)," ,maxmaf=",unique(to_plot$maxmaf),"collapse_method",unique(to_plot$collapse_method),
                    ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                    min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
     # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
     # scale_color_manual(values=c("#926fa8","#47265c","#99176e"))+  
     #guides(colour=FALSE)+
     # scale_color_viridis(option="magma",discrete = T)+
     # scale_shape_manual(values=seq(0,length(unique(to_plot$outlierType))))+
     facet_wrap(~paste0(sex),ncol =3,scales="free_y") + #,scales="free_y"
     geom_hline(yintercept=1,color="red",linetype="dashed") +
     theme(legend.title = element_text(size = 5), legend.text = element_text(size = 5),legend.position = c(6, 0.7))#+
   
   geom_text(aes(x=interaction(maxmaf,chr),y=Risk,label=ifelse(Pval<0.01,
                                                               paste0("n=",num_outliers,"\nRR=",round(Risk,1)))),
             position = position_dodge(width = 1), vjust=-.25,size=2.5) +
     geom_hline(yintercept=1,color="red")
   
   
 # to_plot$chr <- factor(to_plot$chr, levels=c(as.character(1:22),"x"))
ggplot(to_plot,aes(group=sex,fill=sex,label=num_outliers)) + 
    geom_crossbar(aes(x=interaction(maxmaf,chr),y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
    
  #scale_fill_manual(values=c("#97D6F2","#FCF5A9","#B1EAA2" ))+
                       #"#A8C784","#8B3CD2","#226F3A"))+
 #scale_fill_manual(values=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"))+
    # scale_alpha_discrete(range=c(0.2,1))+

  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("chr")+
  ylab("Relative Risk")+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),"cadd=",unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+  
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  # scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+   
  #scale_fill_manual(values=c("#FCF5A9","#97D6F2"))+
  #scale_color_manual(values=c("#dbab3b","#5d8596"))+
 # guides(colour=FALSE)+ facet_wrap(~paste0("maf=",maxmaf),ncol=1)+
    guides(colour=FALSE)+
  facet_wrap(~paste0("z_type=",exp_type)*paste0("z=",z),ncol =3,scales="free_y")+
  #facet_wrap(~paste0("maf=",maxmaf)*paste0("cadd=",cadd),ncol =2,scales="free_y")+
  # guides(colour=FALSE)+ facet_wrap(~paste0("maf=",maxmaf)*paste0("cadd=",cadd),ncol =1,scales="free_y")+
    # ylim(c(0.2,2.7))+
  # geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  # geom_text(aes(group=sex,x=chr,y=Risk,label=num_outliers),  position = position_dodge(width = 1), vjust=-0.25) +
  geom_text(aes(group=sex,x=interaction(maxmaf,chr),y=Risk,label=ifelse(Pval<0.01,
                                                                        paste0("n=",num_outliers,"\nRR=",round(Risk,1),"\np=",
                                                                                           round(Pval,3),"\nsex_padj=",round(padj,2)),
                                                                        paste0("n=",num_outliers,"\np=",round(Pval,3)))),
            position = position_dodge(width = 1), vjust=-.25,size=2.5) +
  # geom_text(aes(group=sex,x=interaction(maxmaf,chr),y=Risk,label=paste0("n=",num_outliers,"\nRR=",round(Risk,1))), 
  #           position = position_dodge(width = 1), vjust=-.25,size=2.5) +
   # geom_text(aes(group=sex,x=as.character(chr),y=Risk,label=paste0("p=",round(Pval,2))),  position = position_dodge(width = 1), vjust=-1.5,size=3) +
    #geom_crossbar(aes(x=XCI_STATUS,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  #geom_crossbar(aes(x=PAR_BINARY,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  geom_hline(yintercept=1,color="red") #+#+

  #stat_compare_means(aes(group = chr,x=Cat,y=Risk), label = "p.format",size=3, method = "anova")

library(facetscales)


####outliers plotting
  to_plot<-all_risks  %>% dplyr::filter(cadd==0&veptype=="all"&gene_window==5000&var_location=="all") %>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  #%>% dplyr::filter(maxmaf=="0.01") %>% dplyr::filter(chr=="5")
  all_risks_renamed<-all_risks
  chr_rename<-c("chrX","chr7","autosomes")
  names(chr_rename)<-c("x","7sub","AllAut")
  all_risks_renamed$chr_renamed<-chr_rename[all_risks_renamed$chr]
  

                                                             
                                                             # %>% dplyr::filter(exp_type=="under" | exp_type=="over") # %>%dplyr::filter(exp_type=="all") #
  to_plot_spread<-to_plot %>% select(exp_type,chr,outlierType,maxmaf,outliers_tested,sex,z,nphen,cadd) %>% spread(exp_type,outliers_tested)
  to_plot_spread$over_prop=to_plot_spread$over/to_plot_spread$all
  to_plot_spread$under_prop=to_plot_spread$under/to_plot_spread$all
  to_plot_gather<-to_plot_spread %>% gather(exp_type,outliers_tested,all:under_prop)
  to_plot_gather_filt<-to_plot_gather %>% dplyr::filter(exp_type=="over_prop" | exp_type=="under_prop") %>% dplyr::filter(outlierType=="outliers")
  
  to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5&veptype=="all"&
                                                                                    gene_window==5000&var_location=="all")%>% 
    dplyr::filter(outlierType=="outliers")%>% #dplyr::filter(sex!="both") %>% 
    # mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% dplyr::filter(exp_type!="all")  %>%
   mutate(chr=fct_relevel(chr,"x","7sub")) %>% mutate(chr=recode(chr,x="chrX","7sub"="chr7 subsetted"))  %>% 
    dplyr::filter(exp_type=="all") %>%dplyr::filter(sex=="allboth") %>%
    dplyr::filter(chr=="autosomes") #dplyr::filter(chr=="autosomes") 
#  dplyr::filter(chr=="chrX" | chr=="chr7 subsetted") #dplyr::filter(chr=="autosomes") 
  ###HERE
  ggplot(to_plot,aes(x=as.factor(chr),y=outliers_tested, group=sex, #linetype
                    fill=sex #,linetype=exp_type #  group=interaction(sex,exp_type)
                   ))+
   #scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+  
   # scale_fill_manual(values=c("#296818","#6fb35d","#abc9a3"))+  
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
    # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme_linedraw(base_size=15)+
     scale_alpha_discrete(range=c(1,0.7))+
    geom_bar(stat="identity", position=position_dodge())+
    # scale_alpha_discrete(range=c(1,0.4,0.1))+
    ##geom_line() +
    #geom_line(stat="identity") +
    #guides(linetype = guide_legend(override.aes = list(fill = NA
    #                                                  , col = "black")))+
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
    ggtitle(paste0("Outliers at MAF ",unique(to_plot$maxmaf), " for Category ", unique(to_plot$CATEGORY), " with outlier type=",paste(unique(to_plot$outlierType),sep=","),
                   ", &+ z=",paste(unique(to_plot$z),collapse=",")," & nphen=",paste(unique(to_plot$nphen),collapse=",")))+
   # geom_text(size=3, aes(group=sex,label=round(outliers_tested,3)), 
    #          position = position_dodge(width = 1),vjust = 0.1) +
    # facet_grid_sc(rows=vars(chr),as.table=F,scales = list(y = list(`chrX`=scale_y_continuous(limits=c(0,70)),
    #                                               `chr7`=scale_y_continuous(limits=c(0,70)),
    #                                               `autosomes`=scale_y_continuous(limits=c(0,1700))
    #                                               )))+
    #facet_wrap(~chr,ncol=3) + #,scales = "free_y") +
    xlab("sex") +ylab("Number of Outliers")
  
  ###this is the stacked one
  ggplot(to_plot%>%dplyr::filter(chr=="x"),aes(x=sex,y=outliers_tested,group=sex,
                     alpha=as.factor(exp_type),
                     fill=sex,size=1))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+
    theme_classic(base_size=10)+
    scale_alpha_discrete(range=c(0.9,0.5))+
    geom_bar(stat="identity",position="stack",group="exp_type") +
    #geom_bar(stat="identity", position=position_dodge()) +
    #guides(linetype = guide_legend(override.aes = list(fill = NA
    #                                                  , col = "black")))+
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+
    ggtitle(paste0("Outliers at MAF ",unique(to_plot$maxmaf), " for Category ", unique(to_plot$CATEGORY), " with outlier type=",paste(unique(to_plot$outlierType),sep=","),
                   ", &+ z=",paste(unique(to_plot$z),collapse=",")," & nphen=",paste(unique(to_plot$nphen),collapse=",")))+
    # geom_text(size=5, aes(group=sex,label=outliers_tested),
    #           position = position_dodge(width = 1),vjust = 0.1) +
    geom_text(size=5, aes(group=exp_type,label=outliers_tested),alpha=1, position = position_stack(vjust = .5))+
    facet_wrap(~as.factor(chr),ncol=3) +
    xlab("sex") +ylab("Number of Outliers")
  ggplot(to_plot,aes(x=(sex),y=outliers_tested,group=exp_type, alpha=exp_type,fill=sex,size=1))+
     scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme_linedraw(base_size=15)+

     #scale_fill_manual(values=c("#dbab3b","#5d8596"))+
    ggtitle(paste0("Outliers at MAF ",unique(to_plot$maxmaf), " for Category ", unique(to_plot$CATEGORY), " with outlier type=",paste(unique(to_plot$outlierType),sep=","),
                   ", &+ z=",paste(unique(to_plot$z),collapse=",")," & nphen=",paste(unique(to_plot$nphen),collapse=",")))+
    geom_bar(stat="identity", position=position_dodge()) +
   #scale_alpha_discrete(range=c(1,0.8,0.2))+
     facet_wrap(~outlierType*as.factor(z),scales="free_y") +xlab("sex") +ylab("Number of Outliers")
  #facet_wrap(~as.factor(paste0("chr: ",chr)*paste0("and ",real_or_null)),scales="free_y")
    # facet_wrap(~as.factor(paste0("nphen=",nphen))) #+ylim(c(0,6700))

  ###CADD
to_plot<-all_risks %>% dplyr::filter(CATEGORY=="all") %>% dplyr::filter(outlierType=="outliersTOP") 
#to_plot$chr <- factor(to_plot$chr, levels=c(as.character(1:22),"x"))
#
#%>% dplyr::filter(z==3)
ggplot(to_plot,aes(group=cadd,alpha=sex)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("CADD")+
  ylab("Relative Risk")+
  scale_alpha_discrete(range=c(0.5,1))+
  ggtitle(paste0("RR z=",unique(to_plot$z)))+
  guides(colour=FALSE)+ labs(fill = "cadd")+ facet_wrap(~maxmaf*sex,scales="free_y")+
  geom_crossbar(aes(x=chr,y=Risk,ymin=Lower,ymax=Upper,fill=as.character(cadd)),position="dodge") +
  geom_hline(yintercept=1,color="red") #+#+



###nphen
to_plot<-all_risks %>% dplyr::filter(CATEGORY=="all") %>% dplyr::filter(outlierType=="outliersTOP")  %>% dplyr::filter(cadd=="15") 
#to_plot$chr <- factor(to_plot$chr, levels=c(as.character(1:22),"x"))
#
#%>% dplyr::filter(z==3)
ggplot(to_plot,aes(group=nphen,alpha=sex)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("CADD")+
  ylab("Relative Risk")+
  scale_alpha_discrete(range=c(0.5,1))+
  ggtitle(paste0("RR z=",unique(to_plot$z)))+
  guides(colour=FALSE)+ labs(fill = "cadd")+ facet_wrap(~maxmaf*sex,scales="free_y")+
  geom_crossbar(aes(x=chr,y=Risk,ymin=Lower,ymax=Upper,fill=as.character(cadd)),position="dodge") +
  geom_hline(yintercept=1,color="red") #+#+

relrisk_final_div_nozero <-relrisk_final_div %>% dplyr::filter(Risk !=0)
ggplot(relrisk_final_div_nozero,aes(fill=sex)) + 
  
  scale_fill_manual(values=c("#97D6F2","#FCF5A9","#B1EAA2"))+
  #                      "#A8C784","#8B3CD2","#226F3A"))+
  #scale_fill_manual(values=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Variant Type")+
  ylab("Relative Risk X/A")+
  labs(fill="Group")+
  ggtitle("Relative Risk X/A Across Variant Types")+
  #ylim(c(1,35))+
  #geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  geom_crossbar(aes(x=Cat,y=Risk,ymin=Lower,ymax=Upper),position="dodge")+
  geom_hline(yintercept=1,color="red") #+

relrisk_final_div_nozero_sex <-relrisk_final_div_sex %>% dplyr::filter(Risk !=0) %>% dplyr::filter(Risk !=Inf)
ggplot(relrisk_final_div_nozero_sex,aes(color=chr, fill=chr) ) + 
  scale_fill_manual(values=c("#AB63EC","#EC7063"))+ #,"#B1EAA2"))+
  # #scale_fill_manual(values=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"))+

   #                      "#A8C784","#8B3CD2","#226F3A"))+
  # scale_fill_manual(values=c("#97D6F2","#FCF5A9","#B1EAA2"))+
  # #                      "#A8C784","#8B3CD2","#226F3A"))+
# scale_fill_manual(values=c("#97D6F2","#FCF5A9")) +#,"#B1EAA2"))+
  #scale_color_manual(values=c("#5d8596","#dbab3b"))+
  scale_color_manual(values=c("#7126A9","#A93226"))+
  #                      "#A8C784","#8B3CD2","#226F3A"))+
  #scale_fill_manual(values=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Variant Type")+
  ylab("Relative Risk F/M")+
  labs(fill="Group")+
  ggtitle("Relative Risk F/M Across Variant Types")+
  guides(colour=FALSE)+ labs(fill = "chr")+
  #ylim(c(0,5))+
  #geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  geom_crossbar(aes(x=Cat,y=Risk,ymin=Lower,ymax=Upper),position="dodge", stat="identity")+
  geom_hline(yintercept=1,color="red") #+
  




  ##plot it all
# ggplot(all_risks %>%dplyr::filter(Cat=="inactive" & maxmaf==0.001 ),aes(fill=sex,group=filt)) + 
 ggplot(all_risks %>%dplyr::filter(filt=="typesBlacklistRemovedALL" & maxmaf=="0.01" ),aes(fill=sex,group=filt)) + 

  scale_fill_manual(values=c("#FCF5A9","#97D6F2")) + #"#B1EAA2"))+
  #                      "#A8C784","#8B3CD2","#226F3A"))+
  #scale_fill_manual(values=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Variant Type")+
  ylab("Relative Risk")+
  labs(fill="Group")+
  ggtitle("Relative Risk Strata: Z=3,NPHEN=5,MAF=0.001")+
  facet_wrap(~sex)+
  #ylim(c(1,35))+
  #geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  geom_crossbar(aes(x=Subregion,y=Risk,ymin=Lower,ymax=Upper),position="dodge")+
  geom_hline(yintercept=1,color="red") #+

 
 ###PAR
 to_plot<-all_risks %>% dplyr::filter(outlierType=="outliers") %>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% 
   dplyr::filter(cadd=="15")  %>% dplyr::filter(exp_type=="over")
 
 ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=sex,color=sex,shape=sex,label=num_outliers)) + 
   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
   geom_line(position=position_dodge(width=0.5))+
   theme_bw(base_size=8)+
   geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
   # geom_crossbar(aes(x=CATEGORY,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
   # geom_crossbar(aes(x=interaction(maxmaf,chr),y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
   theme(axis.text.x = element_text(angle = 45,  hjust=1))+
   xlab("chr")+
   ylab("Relative Risk")+
   labs(fill="Group")+
   ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),"cadd=",unique(to_plot$cadd),
                  ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                  min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
   scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+  
   guides(colour=FALSE)+
   facet_wrap(~paste0("z_type=",exp_type)*paste0('REGION: ',CATEGORY),ncol =3,scales="free_y")+
 geom_text(aes(group=sex,x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,label=ifelse(Pval<0.01,
                                                                         paste0("n=",num_outliers,"\nRR=",round(Risk,1),"\np=",
                                                                                round(Pval,3)),
                                                                         paste0("n=",num_outliers,"\np=",round(Pval,3)))),
             position = position_dodge(width = 1), vjust=-.25,size=2.5) +
   #geom_crossbar(aes(x=XCI_STATUS,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
   geom_hline(yintercept=1,color="red") #+#+
 
 
 ##outliers
 to_plot<-all_risks %>% dplyr::filter(outlierType=="outliers") %>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% 
   dplyr::filter(cadd=="15")
 ggplot(to_plot,aes(x=CATEGORY,y=prop_outliers,color=sex))+
 geom_line(position=position_dodge(width=0.5))+
   theme_bw(base_size=8)+
   scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+ 
   geom_point(position=position_dodge(width=0.5))+facet_wrap(~exp_type)
 
