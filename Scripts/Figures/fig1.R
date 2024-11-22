library("ggplot2")
library(tidyverse)
library(data.table)
library("forcats")
library(ggpubr)
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR"

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
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
# filter_version=c("typesSeenTwice","typesALL","typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
#                  "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
#                  "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesGQ5BlacklistRemovedALL")
filter_version=c("typesALL")
windows=c(500,5000,10000)
collapsetypes="collapsed"
maxmaf<-c("0.1", "0.05","0.01", "0.001","0.0001") #,"0.001","0.0001")
maxtomin_dic<-c("0.05","0.01","0.001","0"     ,"0")
# maxtomin_dic<-c("0","0","0","0","0","0")
names(maxtomin_dic)<-maxmaf
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
# all_risks_p$outlierType <-str_replace(all_risks_p$outlierType,"outliers","")
chr_lens=c(748,972,20868)
names(chr_lens)<-c("x","7","AllAut")

inds_both_half=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_both_half.sex_regress.txt",header=F)%>%pull(V1)
inds_allboth=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_allboth.sex_regress.txt",header=F)%>%pull(V1)


###fig 1a
data_plot_fig1a=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5)%>% 
  dplyr::filter(outlierType=="outliers")%>% dplyr::filter(sex=="allboth") %>% dplyr::filter(veptype=="all"&var_location=="all")%>%
  dplyr::filter(gene_window==5000)%>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>%
  mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% 
  mutate(outliers_per_ind=outliers_tested/length(inds_both_half))%>%
  mutate(outliers_per_gene_per_ind=outliers_per_ind/chr_lens[chr])%>%
  mutate(outliers_per_gene=outliers_tested/chr_lens[chr])%>%
  dplyr::filter(exp_type=="all") # %>%
#dplyr::filter(chr=="autosomes")
#dplyr::filter(chr=="chrX" | chr=="chr7")
# ggplot(plot_fig1a,aes(x=as.factor(chr),y=outliers_tested, group=sex, fill=sex  ))+
plot_fig1a_pt1=ggplot(data_plot_fig1a%>%dplyr::filter(chr!="autosomes"),
                      aes(x=as.factor(chr),y=outliers_tested, group=sex, fill=sex ))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+  
  theme_classic(base_size=20)+
  # ggtitle("Number of outliers across chromosomes")+
  geom_bar(stat="identity", position=position_dodge())+
  xlab("") +ylab("Number Multi-Tissue Outliers")+ guides(fill="none")
plot_fig1a_pt2=ggplot(data_plot_fig1a%>%dplyr::filter(chr=="autosomes"),
                      aes(x=as.factor(chr),y=outliers_tested, group=sex, fill=sex ))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+  
  theme_classic(base_size=20)+
  # ggtitle("Number of outliers across chromosomes")+
  geom_bar(stat="identity", position=position_dodge())+
  xlab("") +ylab("")
ggarrange(plot_fig1a_pt1,plot_fig1a_pt2, widths = c(4, 4),labels="b")


###fig 1b

data_plot_fig1b=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>%
  filter(cadd==15) %>% dplyr::filter(sex=="allboth")  %>% dplyr::filter(z==2.5) %>%# dplyr::filter(chr!="AllAut")%>%
  dplyr::filter(veptype=="all"&gene_window==5000&var_location=="all") %>%filter(exp_type=="all")  %>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes")) %>% dplyr::filter(chr !="8")
plot_fig1b=ggplot(data_plot_fig1b,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=chr,color=chr,shape=chr,label=num_outliers)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_hline(yintercept=1,color="red",linetype="dashed")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_classic(base_size=20)+
  geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  # theme(legend.position="none")+
  scale_y_continuous(trans="log10")+
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  #scale_color_manual(values=c("#99176e","#926fa8","#47265c"))+
   scale_color_manual(values=c("#1B998B","#7E6B8F","#07004D")) #chrs
plot_fig1b 
  # scale_color_manual(values=c("#c7bccf","#926fa8","#47265c","#99176e"))+
  
  #guides(colour=FALSE)+
  #facet_wrap(~paste0(exp_type),ncol =3) + #,scales="free_y"
ggarrange(ggarrange(plot_fig1a_pt1,plot_fig1a_pt2, widths = c(4, 4),labels="b"),
         plot_fig1b,labels=c("b","c"), widths = c(3, 4))
##########supplemental
###fig s1a
data_plot_figs1a=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5)%>% 
  dplyr::filter(outlierType=="outliers")%>% dplyr::filter(sex=="both") %>%
  dplyr::filter(veptype=="all"&gene_window==5000&var_location=="all")%>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% 
  dplyr::filter(exp_type!="all") # %>%
#dplyr::filter(chr=="chrX" | chr=="chr7")
#dplyr::filter(chr=="autosomes") 
###s1a
plot_figs1a_pt1=ggplot(data_plot_figs1a%>%dplyr::filter(chr!="autosomes") ,aes(x=as.factor(exp_type),y=outliers_tested, group=sex, fill=exp_type ))+
  scale_fill_manual(values=c("#6fb35d","#abc9a3"))+  #"#296818",
  # scale_alpha_discrete(range=c(1,0.4,0.1))+
  theme_linedraw(base_size=12)+
  facet_wrap(~chr,ncol=3)+ guides(fill="none") +
  geom_bar(stat="identity", position=position_dodge())+
  xlab("") +ylab("Number of Outliers")
plot_figs1a_pt2=ggplot(data_plot_figs1a%>%dplyr::filter(chr=="autosomes") ,aes(x=as.factor(exp_type),y=outliers_tested, group=sex, fill=exp_type ))+
  scale_fill_manual(values=c("#6fb35d","#abc9a3"))+  #"#296818",
  # scale_alpha_discrete(range=c(1,0.4,0.1))+
  theme_linedraw(base_size=12)+
  facet_wrap(~chr,ncol=3) +
  geom_bar(stat="identity", position=position_dodge())+
  xlab("") +ylab("")
ggarrange(plot_figs1a_pt1,plot_figs1a_pt2, widths = c(4, 4))

data_plot_outlier_parameterizations=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 )%>% 
  dplyr::filter(outlierType=="outliers")%>% dplyr::filter(sex=="both") %>%
  dplyr::filter(veptype=="all"&var_location=="all")%>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% 
  dplyr::filter(exp_type!="all") # %>%
plot_outlier_parameterizations=ggplot(data_plot_outlier_parameterizations,aes(x=z,y=outliers_tested, group=sex, fill=exp_type ))+
  geom_point()
plot_outlier_parameterizations

data_plot_figs1b=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>%filter(outlierType=="outliers") %>%
  filter(cadd==0) %>% dplyr::filter(max_outliers==3) %>% dplyr::filter(sex=="both")  %>% dplyr::filter(z==2.5) %>%filter(exp_type!="all") %>%
  dplyr::filter(veptype=="all"&gene_window==5000&var_location=="all")%>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes")) %>% dplyr::filter(chr !="8")%>% mutate(is_sig=ifelse(Pval<0.01,T,F)) 
plot_figs1b=ggplot(data_plot_figs1b,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=exp_type,alpha=exp_type,color=chr,shape=chr,label=num_outliers)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  scale_alpha_discrete(range = c(1, 0.6), guide = guide_legend(override.aes = list(fill = "black"))) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_classic(base_size=12)+
  geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
  scale_size_discrete(range=c(2,5))+
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  scale_color_manual(values=c("#7E6B8F","#07004D","#1B998B"))+ #chrs
  scale_y_continuous(trans="log10")+
  facet_wrap(~paste0(chr),ncol =3)+
  geom_hline(yintercept=1,color="red",linetype="dashed")
plot_figs1b


