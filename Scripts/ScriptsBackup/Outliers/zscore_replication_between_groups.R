library("readr")
library("tidyverse")
library("data.table")
#outliers_zthresh2_nphen5_globalOutliersRemoved_x_both_half_regress.txt

z=2.5
nphen=3
m_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2.5_nphen3_globalOutliersRemoved_x_m.txt"
f_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2.5_nphen3_globalOutliersRemoved_x_f.txt"
both_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2.5_nphen3_globalOutliersRemoved_x_both_half_regress.txt"

m<-fread(m_file)
f<-fread(f_file)
both<-fread(both_file)


get_sex_vs_both<-function(group1,group2){
  group1_cast<-dcast(cbind(group1[,c(1,2,4,5)],zscore="zscore"),  Gene + Ind ~zscore, value.var="MedZ")
  group2_cast<-dcast(cbind(group2[,c(1,2,4,5)],zscore="zscore"),  Gene + Ind ~zscore, value.var="MedZ")
  group1vsgroup2<-merge(group1_cast,group2_cast,all.x = T)
  zdiff<-abs(group1vsgroup2[,3]-group1vsgroup2[,4])
  group1vsgroup2_full<-group1vsgroup2[,"zdiff":=zdiff]
  return(group1vsgroup2_full)
}

#get the individuals who were actually used in both
both_inds<-unique(both$Ind)
m_inboth<-m %>% dplyr::filter(Ind %in% both_inds)
f_inboth<-f %>% dplyr::filter(Ind %in% both_inds)
#for these groups, get outliers
m_inboth_outlier<-m_inboth%>% dplyr::filter(Y=="outlier")
f_inboth_outlier<-f_inboth%>% dplyr::filter(Y=="outlier")
both_outlier<-both%>% dplyr::filter(Y=="outlier")



#get number of outliers
n_outliers<-data.frame(m=m%>% dplyr::filter(Y=="outlier") %>%nrow(),
                       f=f%>% dplyr::filter(Y=="outlier") %>%nrow(),
                       both=both%>% dplyr::filter(Y=="outlier") %>%nrow())
#get number of outliers from the individuals that are in the both group
inboth_n_outliers<-data.frame(m=m_inboth%>% dplyr::filter(Y=="outlier") %>%nrow(),
           f=f_inboth%>% dplyr::filter(Y=="outlier") %>%nrow(),
           both=both%>% dplyr::filter(Y=="outlier") %>%nrow())

#get outlier comparison
zscores_m_both<-get_sex_vs_both(m_inboth_outlier,both)
zscores_f_both<-get_sex_vs_both(f_inboth_outlier,both)
replicated_n_outliers<-data.frame(m=zscores_m_both %>% dplyr::filter(abs(zscore.x)>3) %>% dplyr::filter(abs(zscore.y)>3) %>%nrow(),
                                  f=zscores_f_both %>% dplyr::filter(abs(zscore.x)>3) %>% dplyr::filter(abs(zscore.y)>3) %>%nrow(),
                                  both=both%>% dplyr::filter(Y=="outlier") %>%nrow())

#get inds that were outliers in given sex and both group

n_outliers_toplot<-data.frame(cbind("n_outliers"=c(as.numeric(inboth_n_outliers),as.numeric(n_outliers),as.numeric(replicated_n_outliers)),
                                "sex"=colnames(inboth_n_outliers),
                                "my_group"=c(rep("in_both",3),rep("original",3),rep("replicated",3))
                                ))

ggplot(n_outliers_toplot,aes(x=sex,y=as.numeric(as.character(n_outliers)),fill=my_group))+
  geom_bar(stat = "identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  ylab("number of outliers")+
  geom_text(aes(label=n_outliers), position=position_dodge(width=0.9), vjust=-0.25)+
  ggtitle("Num Outliers: Z2.5,NPHEN3,X")




ggplot(zscores_f_both,aes(x=zscore.x,y=zscore.y,color=zdiff))+
  theme(legend.position = "none")+
  ggtitle(paste0("Median Z-score Across Tissues: "))+
  geom_abline(intercept = 0, slope = 1,color="blue",alpha=0.2)+
  geom_abline(intercept = 0, slope = 0,color="blue",alpha=0.2)+
  geom_vline(xintercept = 0,color="blue",alpha=0.2)+
  xlab("F")+ylab("Both")+
  # xlim(c(-6,10))+
  # ylim(-6,9)+
  scale_color_gradient2(midpoint = 1, low="gray", high="red",mid="pink")+
  geom_point()
