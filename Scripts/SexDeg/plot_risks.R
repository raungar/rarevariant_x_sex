library("ggplot2")
library(dplyr)
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

mydir="/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8eqtl/RR"


chrgroups=c("aut","x") #,"7") #aut
#chrgroups="aut"
chr_types=chrgroups
tissues<-unique(c("ADPSBQ","ADPSBQ","ADPVSC","ADPVSC","ADRNLG","ADRNLG","ARTAORT","ARTAORT","ARTCRN","ARTCRN","ARTTBL","ARTTBL","BREAST","BREAST","BRNACC","BRNACC","BRNAMY","BRNAMY","BRNCDT","BRNCDT","BRNCHA","BRNCHA","BRNCHB","BRNCHB","BRNCTXA","BRNCTXA","BRNCTXB","BRNCTXB","BRNHPP","BRNHPP","BRNHPT","BRNHPT","BRNNCC","BRNNCC","BRNPTM","BRNPTM","BRNSNG","BRNSNG","BRNSPC","BRNSPC","CLNSGM","CLNSGM","CLNTRN","CLNTRN","ESPGEJ","ESPGEJ","ESPMCS","ESPMCS","ESPMSL","ESPMSL","FIBRBLS","FIBRBLS","HRTAA","HRTAA","HRTLV","HRTLV","KDNCTX","KDNCTX","LCL","LCL","LIVER","LIVER","LUNG","LUNG","MSCLSK","MSCLSK","NERVET","NERVET","PNCREAS","PNCREAS","PTTARY","PTTARY","SKINNS","SKINNS","SKINS","SKINS","SLVRYG","SLVRYG","SNTTRM","SNTTRM","SPLEEN","SPLEEN","STMACH","STMACH","THYROID","THYROID","WHLBLD","WHLBLD"))
tissues="BREAST"
mafmin=0
z=c(2,2.5,3)
z=2.5
nphen=c(3)
mafmin=c("0")
mafmax=c("0.01","0.001","0.0001")
#beta_vals=c("0.10842")
beta_vals=c("0.111") #,"0.04954") #"0.10842",
#mafmax=c(0.001)
cadd=c(0,15)
sex=c("m","f","both")
#sex=c("f")
outliergroups=c("outliers") #,"outliersTOP")
sex_degs=c("male","all","female")
cats=c("lincRNA","all","protein_coding")
# tissue_risks<-data.frame(matrix(ncol = 20))
tissue_risks<-data.frame()
# colnames(tissue_risks)<-c("Risk","Lower","Upper","Pval","Cat","sexdeg","exp_nn","exp_ny","exp_yn","exp_yy",
#                           "cadd","chr",
#                           "tissue","num_outliers","sex","z","nphen","Type","padj","mafmax")
num_tests=length(unique(tissues))*length(z)*length(nphen)*length(sex)*length(sex_degs)*length(cats)
this_nphen=3
for (this_tissue in tissues){
  for(beta_val in beta_vals){
  for(this_group in chrgroups){
  for(this_outlier in outliergroups){
  for(this_sex in sex){
    for(this_chr in chr_types){
      for(this_cadd in cadd){
        for(this_mafmax in mafmax){
          for(this_mafmin in mafmin){
            for(this_z in z){
    #rr_aut_beta0.111_ADPSBQ_z3_nphen2_f_0_0.01.txt #rr_aut_beta0.111_LCL_z4_m_0_0.01.txt
        # file=paste0(mydir,"/rr_",chr_types,"_beta0.111_",this_tissue,"_z",this_z,"_nphen",this_nphen,"_",this_sex,"_",mafmin,"_",mafmax,".txt")
    # rr_outliers_noglobal_varAnnot_medz_zthresh3_nphen3_x_f_beta0.111_cadd0_x_0.111_min0max0.001.txt.gz
         # file=paste0(mydir,"/rr_",this_chr,"_beta0.111_",this_tissue,"_z",this_z,"_",this_sex,"_",mafmin,"_",this_mafmax,"_cadd",this_cadd,".txt")
         file=paste0(mydir,"/rr_",this_outlier,"_noglobal_varAnnot_medz_zthresh",this_z,"_nphen",this_nphen,"_",this_chr,"_",this_sex,
                     "_beta",beta_val,"_cadd",this_cadd,"_",this_chr,"_",beta_val,"_min",this_mafmin,"max",this_mafmax,".txt.gz")
        risks=read.csv(file)
        risks$beta_val<-beta_val
       # risks$mafmax=this_mafmax
        #risks$padj<-p.adjust(risks$Pval,n=num_tests,method="BH")
        ###risks<-risks %>% select(-one_of("sex"))
        #full_risk<-cbind(risks[,-1],sex=this_sex,Chr=chr_types) #Tissue=this_tissue,Chr=chr_types,sex=this_sex,z=this_z,nphen=this_nphen
        # full_risk<-cbind(risks,Chr=chr_types) #Tissue=this_tissue,Chr=chr_types,sex=this_sex,z=this_z,nphen=this_nphen
        tissue_risks<-rbind(tissue_risks,risks)
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
tissue_risks$outliers_tested<-tissue_risks$exp_yn+tissue_risks$exp_yy

#tissue_risks$padj<-p.adjust(tissue_risks$Pval,n=num_tests,method="BH")
#tissue_risks_Rvfile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/all_risks_RVONLY.txt"
#tissue_risks<-fread(tissue_risks_Rvfile,header = T)
#tissue_risks$padj<-p.adjust(tissue_risks$Pval,method="BH")
to_plot<-tissue_risks%>%dplyr::filter(Cat=="all") %>% dplyr::filter(sexdeg_type=="any") %>% dplyr::filter(beta_val==0.111)%>%
  dplyr::filter(outlierType=="outliers") %>% dplyr::filter(chr=="aut")  %>% #dplyr::filter(sex=="m") %>%
  dplyr::filter(cadd==15)# %>% dplyr::filter(max_maf==0.001) # dplyr::filter(beta_val==0.10842 )
#& max_maf ==0.01) %>%

#to_plot$chr<-"NA"
#  dplyr::filter(z==2.5)  #%>% dplyr::filter(sexdeg=="all")
#dplyr::filter(sexdeg=="all")%>%
#  dplyr::filter(nphen==5)%>%dplyr::filter(z==2) %>%dplyr::filter(sex=="both")
#to_plot$padj<-pad.just(to_plot$Pval,method="BH")
ggplot(to_plot,aes(x=as.factor(max_maf),y=Risk,group=as.factor(outlier_direction),fill=sex, alpha=as.factor(outlier_direction))) +#, 
                            # label = ifelse(padj < 0.05, "*", ""))) + 
  ggtitle(paste0("SEXDEG RR for sex: ", unique(to_plot$sex),", z=", unique(to_plot$z),", chr=",unique(to_plot$chr),", method=",unique(to_plot$outlierType),
                ", cadd=",unique(to_plot$cadd),
                ", chr=",unique(to_plot$chr) , " beta=", unique(to_plot$beta_val)))+
                 #", MAF=0.0-",unique(to_plot$max_maf)))+
  scale_alpha_discrete(range=c(1,0.1,0.05))+
  
  #facet_wrap(~sexdeg_type*paste0("maf: ",max_maf),ncol=3,scales = "free_y")+
  facet_wrap(~beta_val*sex,ncol=3)+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("MAF")+
  ylab("Relative Risk")+
  #labs(fill="Group")+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+  
  # scale_fill_manual(values=c("#D1AFF0","#F4B0E0"))+ #
  # scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  #scale_fill_manual(values=c("#FCF5A9","#97D6F2"))+
  #scale_fill_manual(values=c("#e88c46","#9c4706"))+
  #scale_fill_manual(values=c("#d15a5a","#871b1b"))+
  #guides(colour=FALSE)+ #labs(fill = "sexdeg")+
  # ylim(c(0.2, 3))+ #6.1))+
  #geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  #geom_crossbar(aes(x=Cat,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  geom_crossbar(aes(x=as.factor(max_maf),y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  geom_text(size=4,alpha=1, aes(group=outlier_direction,label=ifelse((Pval)<0.05,round(Risk,1),"NS")), position = position_dodge(width = 1),vjust = -0.5) +
#  geom_text(size=10, aes(y = 0.5), position = position_dodge(width = 1),vjust = -0.5) +
  #geom_crossbar(aes(x=PAR_BINARY,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  geom_hline(yintercept=1,color="red") #+  #+#+

  #ylim(c(0.9,1.5))
#stat_compare_means(aes(group = chr,x=Cat,y=Risk), label = "p.format",size=3, method = "anova")
####OUTLIERS
to_plot<-tissue_risks  %>% dplyr::filter(cadd==0) %>%dplyr::filter(CATEGORY=="all") #%>% dplyr::filter(outlierType=="outliersTOP") %>%dplyr::filter(CATEGORY=="all")
##%>% dplyr::filter(maxmaf=="0.01") %>% dplyr::filter(chr=="5")
##all_risks_renamed<-all_risks
#chr_rename<-c("chrX","chr7","autosomes")
#names(chr_rename)<-c("x","7","AllAut")
#all_risks_renamed$chr_renamed<-chr_rename[all_risks_renamed$chr]
to_plot=tissue_risks%>%dplyr::filter(Cat=="all")%>%filter(cadd==15) %>% select(-sexdeg_type) %>% dplyr::filter(z==2.5)

prev_res<-(all_risks) %>% dplyr::filter(chr=="AllAut" & CATEGORY=="all")
colnames(prev_res)[which(colnames(prev_res)=="maxmaf")]<-"max_maf"
colnames(prev_res)[which(colnames(prev_res)=="CATEGORY")]<-"Cat"
prev_tobind<-prev_res %>%select(colnames(to_plot))# %>% select(-real_or_null) %>% select(-discrete_or_continuous)
to_plot_all<-rbind(prev_tobind,to_plot)
  #d#plyr::filter(outlierType=="outliers") #%>% mutate(chr=fct_relevel(chr_renamed,"chrX","chr7","autosomes"))

to_plot=tissue_risks%>%dplyr::filter(Cat=="all")%>%filter(cadd==15) %>%filter(max_maf==0.001)%>% filter(outlierType=="outliers") %>%  #dplyr::filter(chr=="x")
  dplyr::filter(outlier_direction=="all")%>% dplyr::filter(beta_val==0.111) #%>% filter(sexdeg_type=="any")

ggplot(to_plot,aes(x=as.factor(sexdeg_type),y=outliers_tested,group=interaction(sex,sexdeg_type),
                   # alpha=as.factor(sexdeg_type),
                   fill=sex,size=1))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+  
 #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
  theme_linedraw(base_size=10)+
  scale_alpha_discrete(range=c(1,0.4,0.1))+
  geom_bar(stat="identity", position=position_dodge()) +
  #guides(linetype = guide_legend(override.aes = list(fill = NA
   #                                                  , col = "black")))+
  #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
  ggtitle(paste0("Outliers at MAF ",unique(to_plot$maxmaf), " for Category ", unique(to_plot$CATEGORY), " with outlier type=",paste(unique(to_plot$outlierType),sep=","),
                 ", &+ z=",paste(unique(to_plot$z),collapse=",")," & nphen=",paste(unique(to_plot$nphen),collapse=",")))+
  geom_text(size=3, aes(group=sex,label=outliers_tested), 
            position = position_dodge(width = 1),vjust = 0.1) +
  
  facet_wrap(~beta_val*z*chr,ncol=2,scales = "free_y") +xlab("sex") +ylab("Number of Outliers")
#,scales="free_y"
to_plot=tissue_risks%>%dplyr::filter(Cat=="all")%>%filter(cadd==15) %>%filter(max_maf==0.001)%>% filter(outlierType=="outliers" & beta_val==0.111)
ggplot(to_plot,aes(x=as.factor(sexdeg_type),y=outliers_tested,group=interaction(sex,sexdeg_type),
                   alpha=as.factor(outlier_direction),
                   fill=sex,size=1))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+  
  theme_linedraw(base_size=10)+
  scale_alpha_discrete(range=c(1,0.4,0.1))+
  geom_bar(stat="identity", position=position_dodge()) +
 ggtitle(paste0("Outliers at MAF ",unique(to_plot$maxmaf), " for Category ", unique(to_plot$CATEGORY), " with outlier type=",paste(unique(to_plot$outlierType),sep=","),
                 ", &+ z=",paste(unique(to_plot$z),collapse=",")," & nphen=",paste(unique(to_plot$nphen),collapse=",")))+
  geom_text(size=3, aes(group=sex,label=outliers_tested), 
            position = position_dodge(width = 1),vjust = 0.1) +
  
  facet_wrap(~beta_val*z*chr*outlier_direction,ncol=3,scales = "free_y") +xlab("sex") +ylab("Number of Outliers")
#,scales="free_y"

sig<-tissue_risks%>%dplyr::filter(padj<0.05) %>%dplyr::filter(Cat=="all")%>%dplyr::filter(z==2)%>%dplyr::filter(nphen==5)
sig$padj<-p.adjust(sig$Pval)
sig_df<-data.frame(in_male=c(sig%>%dplyr::filter(sex=="m")%>%dplyr::filter(sexdeg=="male_biased")%>%nrow(),
                          sig%>%dplyr::filter(sex=="f")%>%dplyr::filter(sexdeg=="male_biased")%>%nrow(),
                          sig%>%dplyr::filter(sex=="both")%>%dplyr::filter(sexdeg=="male_biased")%>%nrow()),
                   in_female=c(sig%>%dplyr::filter(sex=="m")%>%dplyr::filter(sexdeg=="female_biased")%>%nrow(),
                          sig%>%dplyr::filter(sex=="f")%>%dplyr::filter(sexdeg=="female_biased")%>%nrow(),
                          sig%>%dplyr::filter(sex=="both")%>%dplyr::filter(sexdeg=="female_biased")%>%nrow()),
                   in_all=c(sig%>%dplyr::filter(sex=="m")%>%dplyr::filter(sexdeg=="all")%>%nrow(),
                          sig%>%dplyr::filter(sex=="f")%>%dplyr::filter(sexdeg=="all")%>%nrow(),
                          sig%>%dplyr::filter(sex=="both")%>%dplyr::filter(sexdeg=="all")%>%nrow())
                   )
#sig_df$in_any<-rowSums(sig_df)
rownames(sig_df)<-c("male","female","both")
plot_sig_df<-melt(sig_df); 
colnames(plot_sig_df)<-c("sexdegs","sig_tissues")
plot_sig_df$sex<-rep(c("male","female","both"),3)
ggplot(plot_sig_df,aes(x=sexdegs,y=sex,fill=sig_tissues))+geom_tile()+scale_fill_gradient(low="white", high="#9B7826") +
  geom_text(aes(label = sig_tissues))+ggtitle("Significant Tissues (z=2.5,nphen=2,all)") #+facet_wrap(~z)+
  



##heatmap
library(viridis)
protein_coding_numrvs_allSexDEGs<-tissue_risks%>%dplyr::filter(Cat=="protein_coding" & sexdeg=="all" & Type=="RR_RV_ALLOutliers" & z==4 )
protein_coding_numrvs_allSexDEGs$frac_outlier_rare<-protein_coding_numrvs_allSexDEGs$exp_nn/(protein_coding_numrvs_allSexDEGs$exp_nn+protein_coding_numrvs_allSexDEGs$exp_ny)
protein_coding_numrvs_allSexDEGs$frac_outlier_common<-protein_coding_numrvs_allSexDEGs$exp_ny/(protein_coding_numrvs_allSexDEGs$exp_nn+protein_coding_numrvs_allSexDEGs$exp_ny)
protein_coding_numrvs_allSexDEGs$frac_nonoutlier_rare<-protein_coding_numrvs_allSexDEGs$exp_yn/(protein_coding_numrvs_allSexDEGs$exp_yn+protein_coding_numrvs_allSexDEGs$exp_yy)
protein_coding_numrvs_allSexDEGs$frac_nonoutlier_common<-protein_coding_numrvs_allSexDEGs$exp_yy/(protein_coding_numrvs_allSexDEGs$exp_yn+protein_coding_numrvs_allSexDEGs$exp_yy)
colnames(protein_coding_numrvs_allSexDEGs)[7:10]<-c("outlier_rare_sexdeg","outlier_common_sexdeg","outlier_rare_nonsexdeg","nonoutlier_common_nonsexdeg")
to_plot<- melt(protein_coding_numrvs_allSexDEGs,id.vars =c("Risk","Lower","Upper","Pval","Cat","sexdeg", #"outlier_rare", "outlier_common","nonoutlier_rare","nonoutlier_common",
                                                           "tissue","num_outliers","sex","z","nphen","Type","padj","Chr"))
ggplot(to_plot,aes(x=tissue,y=variable,fill=value))+ #scale_fill_viridis(direction = -1,option = "rainbow") +
  geom_tile()+ggtitle("RR_RV_ALLOutliers: protein_coding (Z=4)")+scale_fill_gradient(low="white", high="orange") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),axis.text.y= element_text(angle = 45,  hjust=1))+
  geom_text(aes(label = round(value,3)),angle=90)+facet_wrap(~sex)
  
  #+#scale_fill_gradient(low="white", high="#9B7826") #+
  # geom_text(aes(label = sig_tissues))+ggtitle("Significant Tissues (z=2.5,nphen=2,all)") #+facet_wrap(~z)+



