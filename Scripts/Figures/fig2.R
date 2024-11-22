library("ggplot2")
library(tidyverse)
library(data.table)
library("forcats")
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"
#####Outliers
sexes=c("m","f","both") #,"allboth")
chrtypes=c("aut","x")
zs=c(2.5)
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered"
all_outliers=data.frame()
for(this_sex in sexes){
  for(this_chrtype in chrtypes){
    for(this_zthresh in zs){
      this_outliers=fread(paste0(mydir,"/outliers_noglobal_medz_zthresh",this_zthresh,
             "_nphen3_",this_chrtype,"_",this_sex,"_maxoutliers3.txt.gz"))
      all_outliers=rbind(all_outliers,
                         cbind(this_outliers,sex=this_sex,chrtype=this_chrtype,zthresh=this_zthresh))
    }
  }
}

##plot distribution of outliers per individual
all_outliers=all_outliers%>%group_by(Ind,chrtype,sex,zthresh)%>%mutate(num_outliers_ind=sum(Y=="outlier"))
all_outliers_summ=all_outliers%>%group_by(Ind,chrtype,sex,zthresh)%>%summarise(num_outliers_ind=sum(Y=="outlier"))
ggplot(all_outliers_summ%>%dplyr::filter(zthresh==2.5),aes(x=num_outliers_ind,color=sex,alpha=0.9))+
  facet_wrap(~chrtype,scales="free")+
  geom_density()+
  xlab("Number of outliers for a given individual")+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+
  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+ 
  theme_classic(base_size=14)
###number of outliers by z thresholding
all_outliers_summ_pergroup=all_outliers%>%group_by(chrtype,sex,zthresh)%>%summarise(num_outliers=sum(Y=="outlier"), num_nonoutliers=sum(Y!="outlier"))
ggplot(all_outliers_summ_pergroup,aes(y=num_outliers,fill=as.factor(zthresh),x=sex))+
  facet_wrap(~chrtype,scales = "free_y")+
  geom_bar(stat="identity",position = 'dodge')+
  scale_fill_manual(values=c("#ed85d6","#c445a8","#850368"))+
  theme_classic(base_size=18)

###fishers test
ft_all=all_outliers_summ_pergroup%>%dplyr::filter((sex=='m'| sex =='f')&zthresh==2.5)
ft_all_contigencytable_x=ft_all%>%dplyr::filter(chrtype=='aut')%>%ungroup%>%select(sex,num_outliers,num_nonoutliers)%>%as.data.frame()
rownames(ft_all_contigencytable_x)=ft_all_contigencytable_x[,'sex']; ft_all_contigencytable_x=ft_all_contigencytable_x%>%select(-sex)
fisher.test(ft_all_contigencytable_x)

ft_all_contigencytable_aut=ft_all%>%dplyr::filter(chrtype=='aut')%>%ungroup%>%select(sex,num_outliers,num_nonoutliers)%>%as.data.frame()
rownames(ft_all_contigencytable_aut)=ft_all_contigencytable_aut[,'sex']; ft_all_contigencytable_aut=ft_all_contigencytable_aut%>%select(-sex)
fisher.test(ft_all_contigencytable_aut)

all_outliers_summ_pergroup_overunder=all_outliers%>%mutate(over_under=ifelse(MedZ>0,"over","under"))%>%
  group_by(chrtype,sex,zthresh,over_under)%>%summarise(num_outliers=sum(Y=="outlier"), num_nonoutliers=sum(Y!="outlier"))
ft_all_overunder=all_outliers_summ_pergroup_overunder%>%dplyr::filter((sex=='f'| sex =='m')&zthresh==2.5)
ft_all_contigencytable_x_under=ft_all_overunder%>%dplyr::filter(chrtype=='aut'&over_under=='over')%>%ungroup%>%select(sex,num_outliers,num_nonoutliers)%>%as.data.frame()
rownames(ft_all_contigencytable_x_under)=ft_all_contigencytable_x_under[,'sex']; ft_all_contigencytable_x_under=ft_all_contigencytable_x_under%>%select(-sex)
fisher.test(ft_all_contigencytable_x_under)


####
# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR"

#groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress")
groups=c("m","f","both","allboth") #, "both_half","both_half.sex","both_half.regress
my_cat=c("xci","par","strata","")
#my_cat=c(".xci",".par",".strata")
my_cat=""
# chr_types=c(as.character(c(1:5,7:11,15:21)),"x") #,
# chr_types=c(as.character(1:22),"x") #,
# chr_types=c(as.character(1:22),"x") #,
chr_types=c("x","AllAut","7") #,

# chr_types=c("x") #,
# chr_types=c("x") #aut
nphen=c(3) #,5)
z=c(2,2.5,3)
z=2.5
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
outlier_types=c("outliers")
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
maxtomin_dic<-c("0.05","0.01","0.001","0.0001","0")
 # maxtomin_dic<-c("0","0","0","0","0")
names(maxtomin_dic)<-as.factor(maxmaf)
maxoutliers=c("3") #,"4")
windows=c(500,5000,10000)
max_outliers=3
collapsetypes="collapsed"
all_risks<-data.frame(matrix(nrow=0,ncol=10))
colnames(all_risks)<-c("Risk","Lower","Upper" ,"Pval","Cat","Type" , "z","nphen","chr","sex")


inds_both_half=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_both_half.txt",header=F)%>%pull(V1)
inds_f=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_f.txt",header=F)%>%pull(V1)
inds_m=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_m.txt",header=F)%>%pull(V1)


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
                      risks$max_outliers<-max_outliers
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
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/ following this
risks_to_compare<-all_risks %>% group_by(outlierType,CATEGORY,z,nphen,cadd,chr,sex,maxmaf,real_or_null,exp_type,var_location,veptype,gene_window,max_outliers) %>% 
  dplyr::summarise(across(c(Risk,Lower,Upper),log)) %>%
  mutate(SE=abs(Lower-Upper)/(2*1.96))%>%dplyr::filter(sex=="m" | sex=="f")
risks_to_compare_side<-pivot_wider(risks_to_compare,names_from=c(sex),
                                   id_cols=c(outlierType,CATEGORY,z,nphen,cadd,chr,maxmaf,real_or_null,exp_type,var_location,veptype,gene_window,max_outliers) ,
                                   values_from=c(Risk,SE)) #id_cols=, names_from=,
risks_to_compare_side$RiskDiff<-risks_to_compare_side$Risk_f-risks_to_compare_side$Risk_m
risks_to_compare_side$SEDiff<-sqrt(risks_to_compare_side$SE_f**2+risks_to_compare_side$SE_m**2)
risks_to_compare_side$zDiff<-(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
risks_to_compare_side$p<-pnorm(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
#n is 2 for cadd * 3 for z * 3 for nphen * 3 for MAF * 24 for chr  (23 + allAut)=1296 OR 54* 2 for outliers
risks_to_compare_side$padj<-p.adjust(risks_to_compare_side$p,method = "BH") #,n=108)
all_risks_p<-merge(all_risks,risks_to_compare_side[,c("outlierType","exp_type","CATEGORY","z","nphen","cadd","chr","maxmaf","p","padj","max_outliers")])

to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F))%>%filter(chr=="7" | chr=="AllAut") %>% 
  filter(outlierType=="outliers")   %>% dplyr::filter(sex=="f")%>% dplyr::filter(var_location!="all") %>% #dplyr::filter(sex=="f") %>% 
  filter(cadd==0)%>% dplyr::filter(z==2.5)%>% dplyr::filter(exp_type=="all") %>% dplyr::filter(max_outliers==3) #%>% dplyr::filter(sex=="both" & exp_type=="all") ### #%>%filter(maxmaf==0.001) %>% dplyr::filter(z==2.5)   %>%   
ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=as.numeric(gene_window),
                   color=as.numeric(gene_window),label=num_outliers)) +  #,shape=chr
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_point(aes(size=as.factor(is_sig)),alpha=0.9,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=8)+
  # scale_size(range = c(2,12))+
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),", CADD=",unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen), ",z=",unique(to_plot$z),", gene window=",unique(to_plot$gene_window),
                 ",maxoutliers=",unique(to_plot$max_outliers),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
   # scale_color_brewer(palette = "RdYlBu")+
  viridis::scale_color_viridis("magma",begin=0,end=0.95,direction=-1)+
  # scale_color_gradient2(low="#abab24",high="#234e82",mid="#23822d",midpoint=5000)+
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
 # scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+ #"#296818",
  #scale_color_manual(values=c("#abab24","#23822d","#234e82","#668223","#23827a"))+
  # scale_color_manual(values=c("#c7bccf","#926fa8","#47265c","#99176e"))+
  # scale_color_manual(values=c("#47265c","#99176e"))+asa
  # scale_color_manual(values=c("#926fa8","#47265c","#99176e"))+  
  #guides(colour=FALSE)+
  facet_wrap(~paste0("var_location=",collapse_method)*paste0("z=",z),ncol =2,scales="free_y") + #
  geom_hline(yintercept=1,color="red",linetype="dashed") +
   #geom_text(aes(label=paste0("RR=",round(Risk,2), "(n=",outliers_tested,")")), size=3, position = position_dodge(width = 1), vjust=-1.25) 
  # geom_text(aes(label=paste0(round(Risk,2))),color=c("black"), size=2, position = position_dodge(width = 0.5), vjust=-1.25)
 geom_text(aes(label=paste0("n=",outliers_tested,"\nRR=",round(Risk,2))),color="black" ,size=2, position = position_dodge(width = 1), vjust=-1.25)


###facet by chr
to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% 
  filter(outlierType=="outliers")%>% # filter(chr=="x" | chr=="AllAut")   %>% #%>% dplyr::filter(sex=="m")
  filter(cadd==0)%>% dplyr::filter(z==2.5)%>% dplyr::filter(exp_type=="all")# %>% 
  dplyr::filter(max_outliers==3) #%>% dplyr::filter(sex=="both" & exp_type=="all") ### #%>%filter(maxmaf==0.001) %>% dplyr::filter(z==2.5)   %>%   
to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F))%>%filter(chr=="7" ) %>% # | chr=="AllAut") %>% 
    filter(outlierType=="outliers")   %>% dplyr::filter(gene_window==5000)%>% dplyr::filter(var_location!="all") %>% #dplyr::filter(sex=="f") %>% 
    filter(cadd==15)%>% dplyr::filter(z==2.5)%>% dplyr::filter(exp_type=="all") %>% dplyr::filter(max_outliers==3) #%>% dplyr::filter(sex=="both" & exp_type=="all") ### #%>%filter(maxmaf==0.001) %>% dplyr::filter(z==2.5)   %>%   
  
ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=var_location,shape=var_location,
                   color=sex,label=num_outliers)) +  #,shape=chr
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_point(aes(size=as.factor(is_sig)),alpha=0.7,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=8)+
  # scale_size(range = c(2,12))+
  
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),", CADD=",unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen), ",z=",unique(to_plot$z),
                 ",maxoutliers=",unique(to_plot$max_outliers),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  # scale_color_manual(values=c("#c7bccf","#926fa8","#47265c","#99176e"))+
  # scale_color_manual(values=c("#47265c","#99176e"))+
  scale_shape_manual(values = c(21:24))+ 
  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
  #guides(colour=FALSE)+
  facet_wrap(~paste0("chr=",chr)*paste0("z=",z)*paste0("sex"=sex),ncol =3,scales="free_y") + #
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  #geom_text(aes(label=paste0("RR=",round(Risk,2), "(n=",outliers_tested,")")), size=3, position = position_dodge(width = 1), vjust=-1.25) 
  # geom_text(aes(label=paste0(round(Risk,2))),color=c("black"), size=2, position = position_dodge(width = 0.5), vjust=-1.25)
  geom_text(aes(label=paste0("n=",outliers_tested,"\nRR=",round(Risk,2))), size=2,color="black", position = position_dodge(width = 1), vjust=-1.25)

####BY DISTANCe
to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F))%>%#filter(chr=="x" | chr=="AllAut") %>% 
  filter(outlierType=="outliers")   %>% dplyr::filter(var_location=="all")%>% dplyr::filter(gene_window==1000) %>% #dplyr::filter(sex=="f") %>% 
  filter(cadd==15)%>% dplyr::filter(z==2.5)%>% dplyr::filter(exp_type=="all") %>% dplyr::filter(max_outliers==3)
ggplot((to_plot),aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,shape=as.factor(is_sig),
                   group=collapse_method,color=(collapse_method), #group=paste0(maxtomin_dic[maxmaf],"-",maxmaf)
                   label=gene_window)) +  #,shape=chr
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  geom_point(aes(size=as.factor(is_sig)),alpha=0.7,position=position_dodge(width=1))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=1)) +
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=7)+  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  ggtitle(paste0("Relative Risk: (outlierType=",unique(to_plot$outlierType),", CADD=",unique(to_plot$cadd),
                 ",nphen=",unique(to_plot$nphen), ",z=",unique(to_plot$z),
                 ",maxoutliers=",unique(to_plot$max_outliers),"), num outlier range [",
                 min(unique(to_plot$num_outliers)),",",max(unique(to_plot$num_outliers)),"]"))+
  # scale_color_manual("MAFs",values=c("#320461","#512480","#6f439c","#a78fbf","#dcc4f5"))+
  # scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
  facet_wrap(~paste0("chr=",chr)*paste0("sex"=sex),ncol =3,scales="free_y") + #
  geom_hline(yintercept=1,color="red",linetype="dashed") + 
  geom_text(aes(label=paste0(round(Risk,1))), 
            size=3, position = position_dodge(width = 1), vjust=-10.25)


##by sex
to_plot=all_risks_p%>%dplyr::filter(CATEGORY=="all")%>% mutate(is_sig=ifelse(Pval<0.01,T,F)) %>% #%>%filter(outlierType=="outliers") %>% dplyr::filter(z==2.5) 
   dplyr::filter(outlierType=="outliers") %>%filter(cadd==0)#%>% # dplyr::filter(exp_type=="all")   
    filter(chr=="AllAut") #filter(maxmaf==0.001)   %>%
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
  facet_wrap(~paste0("z_type=",exp_type)*paste0("z=",z)*paste0("chr=",chr)*paste0("cadd=",cadd),ncol =3,scales="free_y")+
  geom_hline(yintercept=1,color="red",linetype="dashed") 
  
##by sex
to_plot=all_risks%>%dplyr::filter(CATEGORY=="all") %>% #%>%filter(outlierType=="outliers") %>% dplyr::filter(z==2.5) 
  dplyr::filter(outlierType=="outliers"&var_location=="all"&veptype=="all"&sex!='allboth'&exp_type!='all'&gene_window==5000) %>%filter(cadd==15) %>%
  mutate(is_sig=ifelse(Pval<0.01,T,F))#%>% # dplyr::filter(exp_type=="all")   
to_plot=to_plot%>%mutate(adj.pval=Pval*nrow(to_plot))%>%mutate(is_sig=ifelse(adj.pval<0.01,T,F)) 
#filter(chr=="AllAut") #filter(maxmaf==0.001)  
##by over/under
ggplot(to_plot,aes(x=paste0(maxtomin_dic[maxmaf],"-",maxmaf),y=Risk,group=sex,color=sex,shape=sex,label=num_outliers)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
  theme_bw()+
  geom_line(position=position_dodge(width=0.5))+
  theme_bw(base_size=5)+
  geom_point(aes(size=is_sig),position=position_dodge(width=0.5))+ #,
  xlab("MAF")+
  ylab("Relative Risk")+ #xlim(c(0,0.1))+
  labs(fill="Group")+
  scale
  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
  facet_wrap(~paste0(exp_type)*paste0("z=",z)*paste0("chr=",chr),ncol =3,scales="free_y")+
  geom_hline(yintercept=1,color="red",linetype="dashed") 



#OUTLIERS
to_plot<-all_risks  %>% dplyr::filter(cadd==0) %>%dplyr::filter(CATEGORY=="all") %>% dplyr::filter(max_outliers==3) # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
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

fig2b=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5)%>% 
  dplyr::filter(outlierType=="outliers")%>% dplyr::filter(chr !="8") %>%
  mutate(chr=fct_relevel(chr,"x","7","AllAut")) %>% 
  mutate(chr=recode(chr,x="chrX","7"="chr7",AllAut="autosomes"))  %>% 
  dplyr::filter(exp_type!="all")  #%>%  #dplyr::filter(chr=="chrX")# %>%
   # dplyr::filter(chr=="autosomes") #dplyr::filter(chr=="autosomes")
 dplyr::filter(chr=="chrX" | chr=="chr7") #dplyr::filter(chr=="autosomes")

 ggplot(fig2b,aes(x=interaction(sex),y=outliers_tested,group=sex,
  # ggplot(fig2b,aes(x=interaction(sex),y=outliers_per_ind,group=sex,
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
  geom_text(size=3, aes(group=exp_type,label=round(outliers_per_ind,2)),alpha=1, position = position_stack(vjust = .5))+
  facet_wrap(~as.factor(chr),ncol=3,scales="free_y") +xlab("sex") +ylab("Multi-Tissue Outlier Per Individual")

