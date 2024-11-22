library("ggplot2")
library(tidyverse)
library(data.table)
library("forcats")
library(ggpubr)
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR"
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
maxmaf<-c("0.1", "0.05","0.01", "0.001","0.0001") #,"0.001","0.0001")
maxtomin_dic<-c("0.05","0.01","0.001","0"     ,"0")
names(maxtomin_dic)<-maxmaf

get_all_risks<-function(mydir,maxtomin_dic){
    #groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress")
    groups=c("m","f","both","allboth") #, "both_half","both_half.sex","both_half.regress
    my_cat=c("xci","par","strata","")
    #my_cat=c(".xci",".par",".strata")
    my_cat=""
    # chr_types=c(as.character(c(1:5,7:11,15:21)),"x") #,
    # chr_types=c(as.character(1:22),"x") #,
    # chr_types=c(as.character(1:22),"x") #,
    chr_types=c("x","7","AllAut") #,"AllAut","7") #,
    #chr_types=c("x") #,
    # chr_types=c("x") #aut
    nphen=c(3)
    z=c(2,2.5,3)
    z=2.5
    CADD<-c(0,15)
    risk_cat<-c("relative_risk") #,"absolute_risk")
   
    # filter_version=c("typesSeenTwice","typesALL","typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
    #                  "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
    #                  "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
    filter_version=c("typesGQ5BlacklistRemovedALL")
    filter_version=c("typesALL")
    windows=c(500,5000,10000)
    collapsetypes="collapsed"

    # maxtomin_dic<-c("0","0","0","0","0","0")
    #windows=c(100,500,1000,5000,10000)
    max_outliers=3
    #outlier_types=c("outliersTOP","outliers")
    outlier_types=c("outliers")
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
                          if(this_chr=="x"){
                            file=paste0(mydir,"/relative_risk_x_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                        "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
                                        cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
                            this_filt="GQ5BlacklistRemovedALL"
                          }else if(grepl("sub",this_chr)){
                            file=paste0(mydir,"/relative_risk_",this_chr,"_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                        "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",
                                        cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
                          } else{
                            file=paste0(mydir,"/relative_risk_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                        "_",this_group,"_",this_chr,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
                                        "_linc_prot_maxoutliers",max_outliers,"_window",this_window,".txt")
                          }
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
                          all_risks<-rbind(all_risks,risks)
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
    }
    return(all_risks)
  }
get_all_risks_p<-function(all_risks){
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
  all_risks_p<-merge(all_risks,risks_to_compare_side,all=TRUE)
  all_risks_p$gene_window <- factor(all_risks_p$gene_window, levels = c("100","500","1000","5000","10000"))
  return(all_risks_p)
}

all_risks=get_all_risks(mydir,maxtomin_dic)
all_risks_p=get_all_risks_p(all_risks)



### fancy version
plot_sex_stratified<-function(all_risks,this_chr,my_ggtitle){
  to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>%dplyr::filter(var_location=="all") %>%
    dplyr::filter(chr==this_chr)%>%
    dplyr::filter(z==2.5)%>%filter(collapse_method=="collapsed") %>%dplyr::filter(num_outliers>1)%>% dplyr::filter(sex!="both")%>%
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
    ggtitle(my_ggtitle)+
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
    scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
    geom_hline(yintercept=1,color="red",linetype="dashed")
}
plot_sex_stratified(all_risks,this_chr="x",my_ggtitle="Sex-Stratified Enrichments X-Chromosome")
plot_sex_stratified(all_risks,this_chr="AllAut",my_ggtitle="Sex-Stratified Enrichments Autosomes")

  

