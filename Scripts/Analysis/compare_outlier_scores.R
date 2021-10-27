library(data.table) #for fread
library(ggplot2)

dir="/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Outliers"
dir="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8redo/Current"
#dir="/Volumes/groups/smontgom/raungar/Sex/Output/en/Outliers"
#sexdeg_v8/Outliers/outliers_noglobal_varAnnot_medz_zthresh2.5_nphen3_x_f_beta0.111_cadd15_x.txt.gz 
groups=c("both","m","f","both.sex","both.regress", "both_half","both_half.sex","both_half.regress")
groups=c("m","f","both")

types=c("aut","x")
zs=c(2.5)
outliers<-data.frame()
for(this_z in zs){
for (this_group in groups){
  for(this_type in types){
    this_file=paste0(dir,"/outliers_noglobal_medz_zthresh",this_z,"_nphen3_",
                 this_type,"_",this_group,".txt.gz")
    # this_file=paste0(dir,"/outliers_noglobal_varAnnot_medz_zthresh",this_z,"_nphen3_",
    #              this_type,"_",this_group,"_beta0.111_cadd15_",this_type,".txt.gz")
    # varname=paste0(this_type,"_",this_group,"_z",this_z)
    print(this_file)
    df1<-fread(this_file)
    df1$z<-this_z
    df1$sex<-this_group
    df1$aut_or_x<-this_type
    outliers<-rbind(outliers,df1)
    # df2<-cbind(df1,"varname"=varname)
    # assign(varname,df2)
  }
  }
}
library(tidyverse)
outliers_all<-(outliers) ##%>% rename(Ind=ind,Gene=ensg) #%>% dplyr::filter(AbsBetaSex != -Inf)
# %>% mutate(aut_or_x=ifelse(chr=="x","x","aut"))
outliers_spread<-(outliers_all[,c("Ind","Gene","MedZ","z","sex","aut_or_x")]) %>% 
  group_by(Ind,Gene,aut_or_x) %>%
  spread(key=sex,value=MedZ)
outliers_spread$zdiff_m_both<-outliers_spread$m-outliers_spread$both
outliers_spread$zdiff_f_both<-outliers_spread$f-outliers_spread$both
outliers_spread_x<-outliers_spread%>%dplyr::filter(aut_or_x=="x")
outliers_spread_aut<-outliers_spread%>%dplyr::filter(aut_or_x=="aut")

outliers_m<-outliers %>% dplyr::filter(sex=="m") %>% select(ind,ensg,AbsBetaSex,aut_or_x,sex) 
outliers_f<-outliers %>% dplyr::filter(sex=="f") %>% select(ind,ensg,AbsBetaSex,aut_or_x,sex)
outliers_spread_x_m<-outliers_spread_x %>% dplyr::filter(!is.na(m))
outliers_spread_x_f<-outliers_spread_x%>% dplyr::filter(!is.na(f))
outliers_spread_aut_m<-outliers_spread_aut %>% dplyr::filter(!is.na(m))
outliers_spread_aut_f<-outliers_spread_aut%>% dplyr::filter(!is.na(f))

x_m_outliers_withDEGinfo=merge(outliers_spread_x_m,outliers_m, by.x=c("Ind","Gene"),by.y=c("ind","ensg"))
x_f_outliers_withDEGinfo=merge(outliers_spread_x_f,outliers_f, by.x=c("Ind","Gene"),by.y=c("ind","ensg"))
aut_m_outliers_withDEGinfo=merge(outliers_spread_aut_m,outliers_m, by.x=c("Ind","Gene"),by.y=c("ind","ensg"))
aut_f_outliers_withDEGinfo=merge(outliers_spread_aut_f,outliers_f, by.x=c("Ind","Gene"),by.y=c("ind","ensg"))
absmax <- function(x) { x[which.max( abs(x) )]}

x_m_outliers_withDEGinfo_beta=x_m_outliers_withDEGinfo%>% group_by(Gene) %>% summarise_at(c("AbsBetaSex"),absmax)
aut_m_outliers_withDEGinfo_beta=aut_m_outliers_withDEGinfo%>% group_by(Gene) %>% summarise_at(c("AbsBetaSex"),absmax)
x_f_outliers_withDEGinfo_beta=x_f_outliers_withDEGinfo%>% group_by(Gene) %>% summarise_at(c("AbsBetaSex"),absmax)
aut_f_outliers_withDEGinfo_beta=aut_m_outliers_withDEGinfo%>% group_by(Gene) %>% summarise_at(c("AbsBetaSex"),absmax)

x_m_outliers_withDEGinfo_max<-x_m_outliers_withDEGinfo %>% mutate(zdiff_m_both=m-both) %>% group_by(Gene) %>% summarise_at(c("zdiff_m_both"),absmax) 
x_m_outliers_withDEGinfo_max_beta<-merge(x_m_outliers_withDEGinfo_max,x_m_outliers_withDEGinfo_beta,by="Gene")
aut_m_outliers_withDEGinfo_max<-aut_m_outliers_withDEGinfo %>% mutate(zdiff_m_both=m-both) %>% group_by(Gene) %>% summarise_at(c("zdiff_m_both"),absmax) 
aut_m_outliers_withDEGinfo_max_beta<-merge(aut_m_outliers_withDEGinfo_max,aut_m_outliers_withDEGinfo_beta,by="Gene")

x_f_outliers_withDEGinfo_max<-x_f_outliers_withDEGinfo  %>% mutate(zdiff_f_both=f-both) %>% group_by(Gene) %>% summarise_at(c("zdiff_f_both"),absmax) 
x_f_outliers_withDEGinfo_max_beta<-merge(x_f_outliers_withDEGinfo_max,x_f_outliers_withDEGinfo_beta,by="Gene")
aut_f_outliers_withDEGinfo_max<-aut_f_outliers_withDEGinfo  %>% mutate(zdiff_f_both=f-both) %>% group_by(Gene) %>% summarise_at(c("zdiff_f_both"),absmax) 
aut_f_outliers_withDEGinfo_max_beta<-merge(aut_f_outliers_withDEGinfo_max,aut_f_outliers_withDEGinfo_beta,by="Gene")

combined_data<-rbind(cbind("sex"="F","chr"="x",x_f_outliers_withDEGinfo_max_beta %>% rename(zdiff=zdiff_f_both)),
                     cbind("sex"="M","chr"="x",x_m_outliers_withDEGinfo_max_beta%>% rename(zdiff=zdiff_m_both)),
                     cbind("sex"="F","chr"="aut",aut_f_outliers_withDEGinfo_max_beta %>% rename(zdiff=zdiff_f_both)),
                     cbind("sex"="M","chr"="aut",aut_m_outliers_withDEGinfo_max_beta%>% rename(zdiff=zdiff_m_both)) )
# to_plot<-x_f_outliers_withDEGinfo_max_beta %>% mutate(sexDEG=ifelse(AbsBetaSex==-Inf,0,AbsBetaSex))  %>%
#   mutate(zdiff_f_both=ifelse(is.infinite(zdiff_f_both),0,zdiff_f_both)) %>% 
# dplyr::filter( AbsBetaSex<6) # & abs(zdiff_m_both)>1)
to_plot<-combined_data %>% mutate(sexDEG=ifelse(AbsBetaSex==-Inf,0,AbsBetaSex))  %>%
  mutate(zdiff=ifelse(is.infinite(zdiff),0,zdiff)) %>%
dplyr::filter( sexDEG<6 ) # & abs(zdiff_m_both)>1)

get_form<-function(to_plot,sex,chr){
form=paste0("CHR= ",chr,", SEX=", sex, " : ",
            round((lm(abs(zdiff)~sexDEG,to_plot))$coefficients[2],4),"x+",
            round(lm(abs(zdiff)~sexDEG,to_plot)$coefficients[1],4),", R^2=",
            round(summary(lm(abs(zdiff)~sexDEG,to_plot))$r.squared,4))
  return(form)
}
form_x_f<-get_form(to_plot %>% dplyr::filter(chr=="x" & sex=="F"),"F","x")
form_x_m<-get_form(to_plot %>% dplyr::filter(chr=="x" & sex=="M"),"M","x")
form_aut_f<-get_form(to_plot %>% dplyr::filter(chr=="aut" & sex=="F"),"F","aut")
form_aut_m<-get_form(to_plot %>% dplyr::filter(chr=="aut" & sex=="M"),"M","aut")
ggplot(to_plot,aes(x=sexDEG,y=abs(zdiff),color=(sex),,shape=chr,fill=sex,linetype=chr,alpha=0.2))+
  #stat_summary(fun.data=mean_cl_normal) +
  scale_color_manual(values=c("#5d8596","#dbab3b"))+
  scale_fill_manual(values=c("#5d8596","#dbab3b"),aes(alpha=0.15))+
  annotate("text", x=2,y=4,label = form_x_f) +
  annotate("text", x=2,y=3.8,label = form_x_m) +
  annotate("text", x=2,y=3.6,label = form_aut_f) +
  annotate("text", x=2,y=3.4,label = form_aut_m) +
  geom_smooth(method='lm',alpha=0.1)+
  theme_bw()+
  geom_point()+ggtitle("SexDEG Only") #+xlim(c(0,3))

# form=paste0(round((lm(abs(zdiff_f_both)~sexDEG,to_plot))$coefficients[2],4),"x+",
#             round(lm(abs(zdiff_f_both)~sexDEG,to_plot)$coefficients[1],4),", R^2=",
#             round(summary(lm(abs(zdiff_f_both)~sexDEG,to_plot))$r.squared,4))

# ggplot(to_plot,aes(x=sexDEG,y=abs(zdiff_f_both),alpha=0.2))+
#   #stat_summary(fun.data=mean_cl_normal) + 
#   #annotate("text", x=2,y=4,label = form) + 
#   geom_smooth(method='lm')+
#   theme_bw()+
#   geom_point()+ggtitle("x_f_ (SexDEG only)") #+xlim(c(0,3))


xist="ENSG00000229807"
LINC01545="ENSG00000204904"
PUDP="ENSG00000130021"

outliers_spread_x_m<-outliers_spread_x %>% dplyr::filter(!is.na(m)) %>% mutate(isGENE=ifelse(sapply(strsplit(Gene,"\\."),"[[",1)==PUDP,"PUDP","NO"))
outliers_spread_x_f<-outliers_spread_x%>% dplyr::filter(!is.na(f)) %>% mutate(isGENE=ifelse(sapply(strsplit(Gene,"\\."),"[[",1)==PUDP,"PUDP","NO"))

outliers_spread_x_genes<-outliers_spread_x %>% 
  dplyr::filter(abs(zdiff_m_both)>0.5 | abs(zdiff_f_both>0.5 | abs(m)>2.5 | abs(f)>2.5)) %>% 
  # mutate(isGENE=ifelse(sapply(strsplit(Gene,"\\."),"[[",1)==PUDP,"PUDP","NO")) %>%
  mutate(zdiff=max(zdiff_m_both,zdiff_f_both,na.rm = T)) %>% mutate(within_sex=max(f,m,na.rm=T))
outliers_spread_aut_genes<-outliers_spread_aut %>% 
  dplyr::filter(abs(zdiff_m_both)>0.5 | abs(zdiff_f_both>0.5 | abs(m)>2.5 | abs(f)>2.5)) %>% 
  # mutate(isGENE=ifelse(sapply(strsplit(Gene,"\\."),"[[",1)==PUDP,"PUDP","NO")) %>%
  mutate(zdiff=max(zdiff_m_both,zdiff_f_both,na.rm = T)) %>% mutate(within_sex=max(f,m,na.rm=T))



to_plot<-outliers_spread_x_genes %>% dplyr::filter(abs(zdiff_m_both)>0.5 | abs(m)>2.5 | abs(f)>2.5 | abs(both)>2.5 ) #| is_XIST=="PUDP")
topguys=to_plot %>% dplyr::filter(zdiff != -Inf & (abs(within_sex)>2.5 | abs(both)>2.5)) %>% arrange(-abs(zdiff))
topguys=to_plot %>% dplyr::filter((abs(within_sex)>2.5 | abs(both)>2.5)) %>% arrange(-abs(zdiff))
#to_plot<-outliers_spread_x_genes #%>% mutate(isGENE=ifelse(sapply(strsplit(Gene,"\\."),"[[",1)==PUDP,"PUDP","NO")) %>% dplyr::filter( isGENE=="PUDP") #,size=is_XIST,shape=is_XIST)
#ggplot(to_plot,aes(x=within_sex,y=both,color=interaction(sign(zdiff),is.na(m)),alpha=abs(zdiff)))+ #
  ggplot(to_plot,aes(x=within_sex,y=both,color=interaction(is.na(m)),alpha=abs(zdiff)))+
  # geom_abline(linetype="dashed")+
  geom_hline(yintercept = c(-2.5,2.5),linetype="dashed",color="#6b0505")+
  geom_vline(xintercept =  c(-2.5,2.5),linetype="dashed",color="#6b0505")+
  # scale_alpha_continuous(range=c(0.2,1))+
  geom_hline(yintercept = c(0))+
  geom_vline(xintercept =  c(0))+
  geom_point() + # xlim(c(-4.8,4.8))+ ylim(c(-4.8,4.8))+
  theme_bw()+
  ggtitle("Z-scores Across Aut")+
 # scale_color_gradient(low="#296818",high="#5d8596")
  #scale_color_gradient(low="#296818",high="#dbab3b") #female
# scale_color_gradient(values=c("#926fa8","#47265c","#99176e"))  
    scale_color_manual(values=c("#5d8596","#dbab3b"))
 # scale_color_manual(values=c("#296818","#5d8596",
 #                             "#296818","#dbab3b"))
 #  
ggplot(outliers_spread_aut,aes(x=zdiff_m_both))+ggtitle("SexDEGs: Aut ZDiff M-Both") + 
  geom_density()

ggplot(outliers_spread_aut,aes(x=zdiff_m_both))+ggtitle("SexDEGs: Aut ZDiff M-Both") + 
  geom_density()
ggplot(outliers_spread_x,aes(x=abs(zdiff_f_both)))+ggtitle("SexDEGs: ABS X ZDiff F-Both") + 
  geom_histogram(bins = 100)
#
  

gene_SMS<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000102172.15")
gene_XK<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000047597.5")
gene_GAGE10<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000215274.5")
gene_STS<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000101846.6")
gene_SOX3<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000134595.8")
gene_LINC01545<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000204904.7")
gene_ADGRG7<-outliers_spread_aut %>% dplyr::filter(Gene=="ENSG00000144820.7")
gene_NPIPA8<-outliers_spread_aut %>% dplyr::filter(Gene=="ENSG00000214940.8") #GTEX-1GMRU 
gene_ZFX<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000005889.15") #GTEX-1GMRU 
gene_PUDP<-outliers_spread_x %>% dplyr::filter(Gene=="ENSG00000130021.13") #GTEX-1GMRU 
gene_XIST<-outliers_spread_x %>% dplyr::filter(str_detect(Gene,"ENSG00000229807")) #GTEX-1GMRU 

#largediff<-outliers_spread_x %>% dplyr::filter(abs(zdiff_m_both)>2 | abs(zdiff_f_both)> 2)
largediff<-outliers_spread_x %>% dplyr::filter(abs(zdiff_m_both)>2 | abs(zdiff_f_both)> 2)


to_plot<-gather(gene_XIST,sex,zval,both:zdiff_f_both) %>%
  dplyr::filter(sex=="both" | sex=="m" | sex =="f")
ggplot(to_plot,aes(color=sex,x=zval,fill=sex,alpha=0.1))+ geom_density()+
  ggtitle("XIST")+ 
  theme_bw()+
  scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
  scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
  xlim(c(-5,5)) #+
  #geom_jitter(data = to_plot[to_plot$Ind == "GTEX-X585" | to_plot$Ind == "GTEX-16XZZ"| to_plot$Ind == "GTEX-OHPN",],
                  geom_point(data = to_plot[to_plot$Ind == "GTEX-11EMC",],
                              aes(zval,y=0,shape=Ind),height=0.1,size=5,alpha=1)



#







get_melted_df<-function(group1,group2){
  group1_cast<-dcast(group1[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  group1_cast_med<-group1_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group2_cast<-dcast(group2[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  group2_cast_med<-group2_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group1vsgroup2<-merge(group1_cast_med,group2_cast_med)
  # group1vsgroup2<-rbind(group1,group2)[,c(1,2,5,7)]
  # group1vsgroup2_cast<-dcast(group1vsgroup2,  Gene + Ind ~varname, value.var="MedZ")
  # #group1vsgroup2_cast_no_na<-na.omit(group1vsgroup2_cast)
  # group1vsgroup2_cast_no_na_med<-group1vsgroup2_cast_no_na[ , lapply(.SD, median), by = Gene,.SDcols=3:4]
  return(group1vsgroup2)
}

get_sex_vs_both<-function(group1,group2){
  group1_cast<-dcast(group1[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  #group1_cast_med<-group1_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group2_cast<-dcast(group2[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  #group2_cast_med<-group2_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group1vsgroup2<-merge(group1_cast,group2_cast)
  zdiff<-abs(group1vsgroup2[,3]-group1vsgroup2[,4])
  group1vsgroup2_full<-group1vsgroup2[,"zdiff":=zdiff]
  # group1vsgroup2<-rbind(group1,group2)[,c(1,2,5,7)]
  # group1vsgroup2_cast<-dcast(group1vsgroup2,  Gene + Ind ~varname, value.var="MedZ")
  # #group1vsgroup2_cast_no_na<-na.omit(group1vsgroup2_cast)
  # group1vsgroup2_cast_no_na_med<-group1vsgroup2_cast_no_na[ , lapply(.SD, median), by = Gene,.SDcols=3:4]
  return(group1vsgroup2)
}


#just histogram this shit
# z_hist_df<-as.data.frame(rbind(cbind(chr="x",sex="m",zscore=x_m$MedZ),cbind(chr="x",sex="f",zscore=x_f$MedZ),cbind(chr="x",sex="both",zscore=x_both_half.regress$MedZ),
#                  cbind(chr="x",sex="m",zscore=x_m$MedZ),cbind(chr="x",sex="f",zscore=x_f$MedZ),cbind(chr="x",sex="both",zscore=x_both_half.regress$MedZ),
#                  cbind(chr="aut",sex="m",zscore=aut_m$MedZ),cbind(chr="aut",sex="f",zscore=aut_f$MedZ),cbind(chr="aut",sex="both",zscore=aut_both_half.regress$MedZ),
#                   cbind(chr="aut",sex="m",zscore=aut_m$MedZ),cbind(chr="aut",sex="f",zscore=aut_f$MedZ),cbind(chr="aut",sex="both",zscore=aut_both_half.regress$MedZ)))
# #just histogram this shit baby time
# z_hist_df_head<-rbind(data.frame(chr="x",sex="m",zscore=head(x_m$MedZ)),data.frame(chr="x",sex="f",zscore=head(x_f$MedZ)),data.frame(chr="x",sex="both",zscore=head(x_both_half.regress$MedZ)),
#                                     data.frame(chr="x",sex="m",zscore=head(x_m$MedZ)),data.frame(chr="x",sex="f",zscore=head(x_f$MedZ)),data.frame(chr="x",sex="both",zscore=head(x_both_half.regress$MedZ)),
#                                     data.frame(chr="aut",sex="m",zscore=head(aut_m$MedZ)),data.frame(chr="aut",sex="f",zscore=head(aut_f$MedZ)),data.frame(chr="aut",sex="both",zscore=head(aut_both_half.regress$MedZ)),
#                                     data.frame(chr="aut",sex="m",zscore=head(aut_m$MedZ)),data.frame(chr="aut",sex="f",zscore=head(aut_f$MedZ)),data.frame(chr="aut",sex="both",zscore=head(aut_both_half.regress$MedZ)))
# gg_myhist<-ggplot(z_hist_df, aes(x=zscore, color=interaction(chr,sex))) +
#   geom_bar(alpha=0.2, position="identity",fill="white") + 
#   scale_color_manual(values=c("#D5BAE7","#AB63EC","#7126A9",
#                              "#FCC6C0","#EC7063","#A93226"))
#plot(gg_myhist)
#ggsave("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/hist_medz_all.png",
#        width=10,height=10)



png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_x.png")
hist(sample(x_f_z3$MedZ,length(x_both_z3$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(252/255, 245/255, 169/255,.7), main="Histogram of Z-Scores on the X",xlab = "Z-Score")
hist(sample(x_m$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(151/255, 214/255, 242/255,.7),add=T)
hist(x_both$MedZ,breaks=250,xlim=c(-2,2),
     col=rgb(177/255, 234/255, 162/255,.7),add=T)
legend("topright", c("Male", "Female","Both"), col=c("#97D6F2","#FCF5A9","#B1EAA2"), lwd=10)
dev.off()

png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_aut.png")
hist(sample(aut_f$MedZ,length(aut_both$MedZ)),breaks=250,xlim=c(-5,5),
     col=rgb(252/255, 245/255, 169/255,.7),main="Histogram of Z-Scores on the Autosomes",xlab = "Z-Score")
hist(sample(aut_m$MedZ,length(aut_both$MedZ)),breaks=250,xlim=c(-5,5),
     col=rgb(151/255, 214/255, 242/255,.7),add=T)
hist(aut_both$MedZ,breaks=250,xlim=c(-5,5),
     col=rgb(177/255, 234/255, 162/255,.7),add=T)
legend("topright", c("Female", "Male","Both"),
       col=c("#FCF5A9","#97D6F2","#B1EAA2", "#FCC6C0","#EC7063","#A93226"), lwd=10)
dev.off()


png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_aut_x_all.png")
hist(sample(aut_f$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(213/255, 186/255, 231/255,.4),main="Histogram of Z-Scores on the Autosomes and X",xlab = "Z-Score")
hist(sample(aut_m$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(171/255, 99/255, 236/255,.4),add=T)
hist(sample(aut_both$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(113/255, 38/255, 169/255,.4),add=T)
hist(sample(x_m$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(252/255, 198/255, 192/255,.4),add=T)
hist(sample(x_f$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(236/255, 112/255, 99/255,.4),add=T)
hist(x_both$MedZ,breaks=250,xlim=c(-2,2),
     col=rgb(169/255, 50/255, 38/255,.4),add=T)
legend("topright", c("X Female", "X Male","X Both","Autosomes Female", "Autosomes Male","Autosomes Both"),
       col=c("#D5BAE7","#AB63EC","#7126A9", "#FCC6C0","#EC7063","#A93226"), lwd=10)
dev.off()

png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_aut_x_male.png")
hist(sample(aut_m$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(171/255, 99/255, 236/255,.4),main="Histogram of Z-Scores on the Autosomes and X: Male",xlab = "Z-Score")
hist(sample(x_m$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(252/255, 198/255, 192/255,.4),add=T)
legend("topright", c("X Male", "Autosomes Male"),
       col=c("#AB63EC","#EC7063"), lwd=10)
dev.off()

png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_aut_x_both.png")
hist(sample(aut_both$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(113/255, 38/255, 169/255,.4),main="Histogram of Z-Scores on the Autosomes and X: Both",xlab = "Z-Score")
hist(x_both$MedZ,breaks=250,xlim=c(-2,2),
     col=rgb(169/255, 50/255, 38/255,.4),add=F)
legend("topright", c("X Both","Autosomes Both"),
       col=c("#7126A9", "#A93226"), lwd=10)
dev.off()

png("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/HistogramZscores/hist_medz_aut_x_female.png")
hist(sample(aut_f$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(213/255, 186/255, 231/255,.4),main="Histogram of Z-Scores on the Autosomes and X: Females",xlab = "Z-Score")
hist(sample(x_f$MedZ,length(x_both$MedZ)),breaks=250,xlim=c(-2,2),
     col=rgb(236/255, 112/255, 99/255,.4),add=T)
legend("topright", c("X Female", "Autosomes Female"),
       col=c("#D5BAE7", "#FCC6C0"), lwd=10)
dev.off()

if(TRUE){stop("THE ENDS")}
print("DONE BUT DIDNT STOP")

### MEDIAN OF MEDIANS
#X 
x_m_bothregress<-get_melted_df(x_m_z3,x_both_z3)
x_f_bothregress<-get_melted_df(x_f,x_both)
x_m_f<-get_melted_df(x_m,x_f)
x_bothregress_both<-get_melted_df(x_both,x_both)
# x_bothsex_both<-get_melted_df(x_both.sex,x_both)
# x_bothregress_bothsex<-get_melted_df(x_both_half.regress,x_both.sex)



#AUTOSOME
aut_m_both<-get_melted_df(aut_m,aut_both)
aut_f_both<-get_melted_df(aut_f,aut_both)
aut_m_f<-get_melted_df(aut_m,aut_f)
# aut_bothregress_both<-get_melted_df(aut_both_half.regress,aut_both)
# aut_bothsex_both<-get_melted_df(aut_both.sex,aut_both)
# aut_bothregress_bothsex<-get_melted_df(aut_both_half.regress,aut_both.sex)
#   
####JUST MEDIANS, KEEP INDIVIDUALS
###X
#SAME NUMBER (both=F+M)
x_m_bothregress_med<-get_sex_vs_both(x_m,x_both)
x_f_bothregress_med<-get_sex_vs_both(x_f,x_both)
# x_bothregress_both_med<-get_sex_vs_both(x_both_half.regress,x_both_half.regress)
# x_bothsex_both_med<-get_sex_vs_both(x_both.sex,x_both_half.regress)
# x_bothregress_bothsex_med<-get_sex_vs_both(x_both_half.regress,x_both_half.regress.sex)
# #half (both=0.5F+0.5M)
# x_m_halfboth.regress_med<-get_sex_vs_both(x_m,x_both_half.regress)
# x_f_halfboth.regress_med<-get_sex_vs_both(x_f,x_both_half.regress)
# # x_halfboth.regress_halfboth_med<-get_sex_vs_both(x_both_half.regress,x_both_half)
# x_halfboth.sex_halfboth_med<-get_sex_vs_both(x_both_half.sex,x_both_half)
# # x_halfboth.regress_halfboth.sex_med<-get_sex_vs_both(x_both_half.regress,x_both_half.sex)
# #half compare to full
# x_halfboth.regress_both.regress_med<-get_sex_vs_both(x_both_half.regress,x_both_half.regress)
# x_halfboth.sex_both.sex_med<-get_sex_vs_both(x_both_half.sex,x_both.sex)
# x_halfboth_both.sex_med<-get_sex_vs_both(x_both_half,x_both)

####AUT
aut_m_bothregress_med<-get_sex_vs_both(aut_m,aut_both)
aut_f_bothregress_med<-get_sex_vs_both(aut_f,aut_both)
# aut_bothregress_both_med<-get_sex_vs_both(aut_both_half.regress,aut_both)
# aut_bothsex_both_med<-get_sex_vs_both(aut_both.sex,aut_both)
# aut_bothregress_bothsex_med<-get_sex_vs_both(aut_both_half.regress,aut_both.sex)
# #half (both=0.5F+0.5M)
# aut_m_halfboth.regress_med<-get_sex_vs_both(aut_m,aut_both_half.regress)
# aut_f_halfboth.regress_med<-get_sex_vs_both(aut_f,aut_both_half.regress)
# aut_halfboth.regress_halfboth_med<-get_sex_vs_both(aut_both_half.regress,aut_both_half)
# aut_halfboth.sex_halfboth_med<-get_sex_vs_both(aut_both_half.sex,aut_both)
# aut_halfboth.regress_halfboth.sex_med<-get_sex_vs_both(aut_both_half.regress,aut_both_half.sex)
# #half compare to full
# aut_halfboth.regress_both.regress_med<-get_sex_vs_both(aut_both_half.regress,aut_both_half.regress)
# aut_halfboth.sex_both.sex_med<-get_sex_vs_both(aut_both_half.sex,aut_both.sex)
# aut_halfboth_both.sex_med<-get_sex_vs_both(aut_both_half,aut_both)

highz<-(apply(x_m_bothregress_med,1,function(x){any(abs(as.numeric(x[5]))>2)}))
highz_genes<-x_m_bothregress_med[highz,]
highz_genes

plot_list<-list("x_m_bothregress"=x_m_bothregress,"x_f_bothregress"=x_f_bothregress,"x_m_f"=x_m_f) #,
                #"aut_m_bothregress"=aut_m_bothregress,"aut_f_bothregress"=aut_f_bothregress,"aut_m_f"=aut_m_f,
                #"aut_bothregress_both"=aut_bothregress_both,"aut_bothsex_both"=aut_bothsex_both,"aut_bothregress_bothsex"=aut_bothregress_bothsex)

plot_list_med<-list("x_m_bothregress"=x_m_bothregress_med,"x_f_bothregress"=x_f_bothregress_med,
                "x_bothregress_both"=x_bothregress_both_med,"x_bothsex_both"=x_bothsex_both_med,"x_bothregress_bothsex"=x_bothregress_bothsex_med,
                "aut_m_bothregress"=aut_m_bothregress_med,"aut_f_bothregress"=aut_f_bothregress_med,
                "aut_bothregress_both"=aut_bothregress_both_med,"aut_bothsex_both"=aut_bothsex_both_med,"aut_bothregress_bothsex"=aut_bothregress_bothsex_med)

plot_list_med_half<-list(
  # "aut_m_halfboth.regress_med"=aut_m_halfboth.regress_med,"aut_f_halfboth.regress_med"=aut_f_halfboth.regress_med,
  #    "aut_halfboth.regress_halfboth_med"=aut_halfboth.regress_halfboth_med,"aut_halfboth.sex_halfboth_med"=aut_halfboth.sex_halfboth_med,
  #    "aut_halfboth.regress_halfboth.sex_med"=aut_halfboth.regress_halfboth.sex_med,"aut_halfboth.regress_both.regress_med"=aut_halfboth.regress_both.regress_med,
  #    "aut_halfboth.sex_both.sex_med"=aut_halfboth.sex_both.sex_med,"aut_halfboth_both.sex_med"=aut_halfboth_both.sex_med,
  #    
     "x_m_halfboth.regress_med"=x_m_halfboth.regress_med,"x_f_halfboth.regress_med"=x_f_halfboth.regress_med,
     "x_halfboth.regress_halfboth_med"=x_halfboth.regress_halfboth_med,"x_halfboth.sex_halfboth_med"=x_halfboth.sex_halfboth_med,
     "x_halfboth.regress_halfboth.sex_med"=x_halfboth.regress_halfboth.sex_med,"x_halfboth.regress_both.regress_med"=x_halfboth.regress_both.regress_med,
     "x_halfboth.sex_both.sex_med"=x_halfboth.sex_both.sex_med,"x_halfboth_both.sex_med"=x_halfboth_both.sex_med
     )

tmp_plot<-list("x_m_bothregress_med"=x_m_bothregress_med)

my_r2<-cor(aut_halfboth.regress_halfboth.sex_med[,3],aut_halfboth.regress_halfboth.sex_med[,4])^2
my_r2
this_plot=aut_halfboth.regress_halfboth.sex_med
i<-0
filter_me=F
for (this_plot in plot_list){
  i<-i+1
  this_name<-names(tmp_plot)[i]
  all_nrow<-nrow(this_plot)
  if(filter_me==T){
    this_plot<-this_plot %>% dplyr::filter(zdiff>0.02)
    if(nrow(this_plot)==0){next}
  }
  percent_passed=paste0(round(nrow(this_plot)/all_nrow*100,2),"%")
  g<-ggplot(this_plot,
         aes_string(x=colnames(this_plot)[3], y=colnames(this_plot)[2],
                    alpha=0.2))+
     theme(legend.position = "none")+
    ggtitle(paste0("Median Z-score Across Tissues: ",this_name, " (",percent_passed,")"))+
    geom_abline(intercept = 0, slope = 1,color="blue",alpha=0.2)+
    geom_abline(intercept = 0, slope = 0,color="blue",alpha=0.2)+
    geom_vline(xintercept = 0,color="blue",alpha=0.2)+
    # xlim(c(-6,10))+
    # ylim(-6,9)+
    geom_point()
#  ggsave(paste0("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/MedzCompare/",this_name,"_genemedz_only_zdiff1.5.png"),
 #        width=10,height=10)
#  ggsave(paste0("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/MedzCompare/",this_name,"_genemedz_only_zdiff1.5.png"),
 #        width=10,height=10)
}
