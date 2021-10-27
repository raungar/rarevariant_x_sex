library("ggplot2")
library(tidyverse)
library("forcats")
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8redo/RR"

#groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress")
groups=c("m","f","both") #, "both_half","both_half.sex","both_half.regress
my_cat=c("xci","par","strata","")
#my_cat=c(".xci",".par",".strata")
my_cat=""
# chr_types=c(as.character(c(1:5,7:11,15:21)),"x") #,
# chr_types=c(as.character(1:22),"x") #,
# chr_types=c(as.character(1:22),"x") #,
chr_types=c("x","7","AllAut") #,"AllAut","7") #,
#chr_types=c("x") #,
# chr_types=c("x") #aut
nphen=c(3) #,5)
z=c(2,2.5,3)
# z=2.5
CADD<-c(0,15)
risk_cat<-c("relative_risk") #,"absolute_risk")
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
filter_version=c("typesSeenTwice","typesALL","typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                 "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                 "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                 "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                 "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesGQ5BlacklistRemovedALL")
filter_version=c("typesALL")
outlier_types=c("outliersTOP","outliers")
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
maxtomin_dic<-c("0.05","0.01","0.001","0.0001","0")
maxtomin_dic<-c("0","0","0","0","0")
names(maxtomin_dic)<-as.factor(maxmaf)
all_risks<-data.frame(matrix(nrow=0,ncol=10))
colnames(all_risks)<-c("Risk","Lower","Upper" ,"Pval","Cat","Type" , "z","nphen","chr","sex")

for(this_cat in my_cat){
  for(this_outlier in outlier_types){
    print(this_outlier)
    for(this_maxmaf in maxmaf){
      for (this_group in groups){
        for(this_chr in chr_types){
          for(this_nphen in nphen){
            for(this_z in z){
              for(this_filt in filter_version){
                for(cadd_min in CADD){
                  # if(this_chr=="x"){
                  #   file=paste0(mydir,"/relative_risk_x_",this_outlier,"_z",this_z,"_nphen",this_nphen,"_",this_chr,"_",this_group,
                  #               "_min0max",this_maxmaf,"_CADD","typesGQ5BlacklistRemovedALL","_CADD",cadd_min,"_linc_prot",this_cat,".csv")
                  #   #relative_risk_x_outliers_z2.5_nphen5_x_both_min0max0.001_CADDtypesGQ5BlacklistRemovedALL_CADD15_linc_prot.csv
                  #   this_filt="GQ5BlacklistRemovedALL"
                  # }else{
                  # file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,"_",this_group,"_",this_chr,
                  #             "_min0max",this_maxmaf,"_CADD",this_filt,"_CADD",cadd_min,"_linc_prot",this_cat,".txt")
                  # } 
                  if(this_chr=="x"){
                    file=paste0(mydir,"/relative_risk_x_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
                                cadd_min,"_linc_prot.csv")
                    this_filt="GQ5BlacklistRemovedALL"
                    #_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                    #"-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                    #this_isnull,"-",use_chrgroup,"-binary.txt.gz")
                  }else{
                    file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                "_",this_group,"_",this_chr,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
                                "_linc_prot.txt")
                  }
                  
                  risks=read.csv(file[1])
                  #load(file)
                  risks$z<-this_z
                  risks$nphen<-this_nphen
                  risks$chr<-this_chr
                  risks$sex<-this_group
                  risks$maxmaf<-this_maxmaf
                  risks$filt<-this_filt
                  risks$cadd<-cadd_min
                  risks$nphen<-this_nphen
                  risks$thiscat<-this_cat
                  risks$outlierType<-this_outlier
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


all_risks$outliers_tested<-all_risks$exp_yn+all_risks$exp_yy
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/ following this
risks_to_compare<-all_risks %>% group_by(outlierType,CATEGORY,z,nphen,cadd,chr,sex,maxmaf,real_or_null,exp_type) %>% 
  dplyr::summarise(across(c(Risk,Lower,Upper),log)) %>%
  mutate(SE=abs(Lower-Upper)/(2*1.96))%>%dplyr::filter(sex=="m" | sex=="f")
risks_to_compare_side<-pivot_wider(risks_to_compare,names_from=c(sex),
                                   id_cols=c(outlierType,CATEGORY,z,nphen,cadd,chr,maxmaf,real_or_null,exp_type) ,
                                   values_from=c(Risk,SE)) #id_cols=, names_from=,
risks_to_compare_side$RiskDiff<-risks_to_compare_side$Risk_f-risks_to_compare_side$Risk_m
risks_to_compare_side$SEDiff<-sqrt(risks_to_compare_side$SE_f**2+risks_to_compare_side$SE_m**2)
risks_to_compare_side$zDiff<-(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
risks_to_compare_side$p<-pnorm(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
#n is 2 for cadd * 3 for z * 3 for nphen * 3 for MAF * 24 for chr  (23 + allAut)=1296 OR 54* 2 for outliers
risks_to_compare_side$padj<-p.adjust(risks_to_compare_side$p,method = "BH") #,n=108)
all_risks_p<-merge(all_risks,risks_to_compare_side[,c("outlierType","exp_type","CATEGORY","z","nphen","cadd","chr","maxmaf","p","padj")])

to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% filter(outlierType=="outliers") %>%
  filter(cadd==15)%>% dplyr::filter(z==2.5)# %>%# dplyr::filter(sex=="f")## %>% filter(chr=="x") #%>%filter(maxmaf==0.001) %>% dplyr::filter(z==2.5)   %>% dplyr::filter(sex=="both")  
ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=chr,color=chr,shape=chr,label=num_outliers)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=8)+
  geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),", CADD=",unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  scale_color_manual(values=c("#926fa8","#47265c","#99176e"))+  
  #guides(colour=FALSE)+
  facet_wrap(~paste0("z_type=",exp_type)*paste0("z=",z)*paste0("sex"=sex),ncol =3,scales="free_y") + #
  geom_hline(yintercept=1,color="red",linetype="dashed")# +
  geom_text(aes(label=paste0("RR=",round(Risk,2), "(n=",outliers_tested,")")), size=3, position = position_dodge(width = 1), vjust=-1.25) 
#geom_text(aes(label=paste0("n=",outliers_tested,"\nRR=",round(Risk,2))),  position = position_dodge(width = 1), vjust=-1.25) 

to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% #%>%filter(outlierType=="outliers") %>% dplyr::filter(z==2.5) 
   dplyr::filter(outlierType=="outliers")%>% dplyr::filter(exp_type=="all")   %>% filter(chr=="7") #%>%filter(maxmaf==0.001) filter(cadd==15)  %>%
ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=sex,color=sex,shape=sex,label=num_outliers)) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=5)+
  geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),", CADD=", unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen),",z=",unique(to_plot$z),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  geom_text(aes(label=paste0("RR=",round(Risk,2), "\n(n=",outliers_tested,")")), size=2, position = position_dodge(width = 1), vjust=-1.25) +

  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+ #"#296818",
  facet_wrap(~paste0("z_type=",exp_type)*paste0("z=",z)*paste0("chr=",chr)*paste0("cadd=",cadd),ncol =2,scales="free_y")+
  geom_hline(yintercept=1,color="red",linetype="dashed") 
  


#OUTLIERS
to_plot<-all_risks  %>% dplyr::filter(cadd==0) %>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
#%>% dplyr::filter(maxmaf=="0.01") %>% dplyr::filter(chr=="5")
all_risks_renamed<-all_risks
chr_rename<-c("chrX","chr7","autosomes")
names(chr_rename)<-c("x","7","AllAut")
all_risks_renamed$chr_renamed<-chr_rename[all_risks_renamed$chr]



# %>% dplyr::filter(exp_type=="under" | exp_type=="over") # %>%dplyr::filter(exp_type=="all") #
to_plot_spread<-to_plot %>% select(exp_type,chr,outlierType,maxmaf,outliers_tested,sex,z,nphen,cadd) %>% spread(exp_type,outliers_tested)
to_plot_spread$over_prop=to_plot_spread$over/to_plot_spread$all
to_plot_spread$under_prop=to_plot_spread$under/to_plot_spread$all
to_plot_gather<-to_plot_spread %>% gather(exp_type,outliers_tested,all:under_prop)
to_plot_gather_filt<-to_plot_gather %>% dplyr::filter(exp_type=="over_prop" | exp_type=="under_prop") %>% dplyr::filter(outlierType=="outliers")

to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5)%>% 
  dplyr::filter(outlierType=="outliers")%>% #dplyr::filter(sex!="both") %>% 
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% dplyr::filter(exp_type!="all") # %>%
  #  dplyr::filter(chr=="autosomes") #dplyr::filter(chr=="autosomes") 
 # dplyr::filter(chr=="chrX" | chr=="chr7") #dplyr::filter(chr=="autosomes") 

ggplot(to_plot,aes(x=interaction(sex),y=outliers_tested,group=sex,
                   alpha=as.factor(exp_type),
                   fill=sex,size=1))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
  #scale_fill_manual(values=c("#dbab3b","#5d8596"))+
  #scale_fill_manual(values=c("#dbab3b","#5d8596"))+
  theme_linedraw(base_size=10)+
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
  facet_wrap(~as.factor(chr),ncol=3,scales="free_y") +xlab("sex") +ylab("Number of Outliers")

