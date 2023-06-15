library("ggplot2")
library(tidyverse)
library(data.table)
library(reshape2)
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

all_outliers_x<-dplyr::filter(all_outliers,chrtype=="x")
all_outliers_aut<-dplyr::filter(all_outliers,chrtype=="aut")
all_outliers_7<-dplyr::filter(all_outliers,chr=="chr7")

to_plot_all_outliers_aut<-dcast(all_outliers_aut,value.var = "MedZ",Ind+Gene~sex)%>%
  dplyr::mutate(sex_stratified=ifelse(is.na(f),"male","female"))%>%mutate(sex_stratified_score=ifelse(sex_stratified=="male",m,f))%>%
  mutate(zdiff=sex_stratified_score-both) %>%
  dplyr::filter(abs(both)>2.5 | abs(sex_stratified_score)>2.5 |
                  (abs(zdiff)>0.5 & abs(sex_stratified_score)<2.5 & abs(both)<2.5))
to_plot_all_outliers_x<-dcast(all_outliers_x,value.var = "MedZ",Ind+Gene~sex)%>%
  dplyr::mutate(sex_stratified=ifelse(is.na(f),"male","female"))%>%mutate(sex_stratified_score=ifelse(sex_stratified=="male",m,f))%>%
  mutate(zdiff=sex_stratified_score-both) %>%
  dplyr::filter(abs(both)>2.5 | abs(sex_stratified_score)>2.5 |
                  (abs(zdiff)>0.5 & abs(sex_stratified_score)<2.5 & abs(both)<2.5))

to_plot_all_outliers_x%>%dplyr::filter(abs(both)>2.5 & abs(sex_stratified_score)<=2.5)%>%nrow() #gt in both
to_plot_all_outliers_x%>%dplyr::filter(abs(both)<=2.5 & abs(sex_stratified_score)>2.5)%>%nrow() #gt in sex-stratified
to_plot_all_outliers_aut%>%dplyr::filter(abs(both)>2.5 & abs(sex_stratified_score)<=2.5)%>%nrow() #gt in both
to_plot_all_outliers_aut%>%dplyr::filter(abs(both)<=2.5 & abs(sex_stratified_score)>2.5)%>%nrow() #gt in sex-stratified

#sex-stratified direction 
(to_plot_all_outliers_x%>%dplyr::filter(abs(both)>2.5 & abs(sex_stratified_score)<=2.5)%>%nrow()+to_plot_all_outliers_x%>%dplyr::filter(abs(both)<=2.5 & abs(sex_stratified_score)>2.5)%>%nrow())/nrow(to_plot_all_outliers_x)
(to_plot_all_outliers_aut%>%dplyr::filter(abs(both)>2.5 & abs(sex_stratified_score)<=2.5)%>%nrow()+to_plot_all_outliers_aut%>%dplyr::filter(abs(both)<=2.5 & abs(sex_stratified_score)>2.5)%>%nrow())/nrow(to_plot_all_outliers_x)

ggplot((to_plot_all_outliers_x),aes(x=sex_stratified_score,y=both,color=sex_stratified))+
  theme_bw(base_size = 20)+
  geom_hline(yintercept=-2.5,color="#b02010")+geom_hline(yintercept=2.5,color="#b02010")+geom_vline(xintercept=-2.5,color="#b02010")+geom_vline(xintercept=2.5,color="#b02010")+
  xlab("z-score sex-stratified")+ylab("z-score both")+
  geom_point(aes(alpha=0.7))+
  ggtitle("X-chromosome")+
  scale_color_manual(values=c("#dbab3b","#5d8596")) #"#296818",
plot_gene_distributions<-function(outliers,gene,ind,gene_name,inds_f,inds_m){
  # this_gene=outliers%>%dplyr::filter(Gene==gene)
  # myplot=ggplot(this_gene,aes(x=MedZ,fill=sex,color=sex))+geom_density(aes(alpha=0.5))+
  #   scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
  #   scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
  #   geom_point(data = outliers%>%dplyr::filter(Ind==ind&Gene==gene),size=7,aes(color=sex,y=0))+
  #   theme_bw(base_size = 20)+
  #   ggtitle(gene_name)
  this_gene=outliers%>%dplyr::filter(Gene==gene)
  myplot=ggplot(this_gene,aes(x=MedZ,fill=sex,color=sex))+geom_density(aes(alpha=0.5))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
    geom_point(data = outliers%>%dplyr::filter(Ind==ind&Gene==gene),size=7,aes(color=sex,y=0))+
    theme_bw(base_size = 20)+
    ggtitle(gene_name)
  return(myplot)
}
inds_both_half=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_both_half.txt",header=F)%>%pull(V1)
inds_f=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_f.txt",header=F)%>%pull(V1)
inds_m=fread("/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/gtex_2017-06-05_v8_individuals_passed_m.txt",header=F)%>%pull(V1)

plot_gene_distributions(all_outliers_x,"ENSG00000171388.11","GTEX-U8XE","APLN",inds_f,inds_m)
plot_gene_distributions(all_outliers_x,"ENSG00000189108.12","GTEX-13O3Q","IL1RAPL2")
plot_gene_distributions(all_outliers_x,"ENSG00000187268.11","GTEX-1GF9X","FAM9C")
plot_gene_distributions(all_outliers_x,"ENSG00000214827.9","GTEX-1GF9X","MTCP1")
plot_gene_distributions(all_outliers_aut,"ENSG00000152254.10","GTEX-1GMR8","G6PC2")
plot_gene_distributions(all_outliers_aut,"ENSG00000171435.13","GTEX-QCQG","KSR2")
plot_gene_distributions(all_outliers_aut,"ENSG00000166118.7","GTEX-131YS","SPATA19")
plot_gene_distributions(all_outliers_aut,"ENSG00000198610.10","GTEX-YFC4","AKR1C4")
