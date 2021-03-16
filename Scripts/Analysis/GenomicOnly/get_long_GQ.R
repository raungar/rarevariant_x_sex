library(ggplot2)
library(GenomicRanges)
library("RColorBrewer")
library("data.table")
library("crunch")
library(dplyr)
#read everything in

#get par region
make_subregion_df<-function(){
  #PAR REGIONS
  par1_s=10001; par1_e=2781479; par2_s=155701383; par2_e=156030895; nonpar_s=2781480; nonpar_e=155701382;
  xar_s=2731479;xar_e=58555579
  xcr1_s=62462543; xcr1_e=89140830;
  xtr_s=89140830;xtr_e=93428068;
  xcr2_s=93428068;xcr2_e=155701383
  len_par1=par1_e-par1_s;len_nonpar=par2_s-par1_e; len_par2=par2_e-par2_s
  len_xcr1=xcr1_e-xcr1_s; len_xcr2=xcr2_e-xcr2_s
  len_xtr=xtr_e-xtr_s; len_xar=xar_e-xar_s; 
  
  region_len<-c(len_par1, len_nonpar,len_par2,len_xcr1,len_xcr2,len_xtr,len_xar)
  names(region_len)<-c("PAR1","NONPAR","PAR2", "XCR1","XCR2","XTR","XAR")
  start=c(par1_s,nonpar_s,par2_s,xcr1_s,xcr2_s,xtr_s,xar_s)
  end=c(par1_e,nonpar_e,par2_e,xcr1_e,xcr2_e,xtr_e,xar_e)
  
  subregion_df<-cbind(data.frame(region_len),start,end)
  subregion_df$region<-rownames(subregion_df)
  
  par_hash=region_len
  return(list(subregion_df, par_hash))
  
}
subregion_df<-cbind(make_subregion_df()[[1]],chr="chrX")
subregion_df<-subregion_df[c("PAR1","PAR2","NONPAR","XTR"),]
par_hash<-make_subregion_df()[[2]]


output_dir<-"/oak/stanford/groups/smontgom/raungar/Sex/Output"
if(FALSE){
genotypeq<-fread(paste0(output_dir,"/analysis_v8/genomic_only/qualitychecks/chrX_GQ_MAF.txt"),fill=T)
colnames(genotypeq)[1:3]<-c("chrom","pos","MAF")
sex<-fread("/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt")
sex_dic<-sex$SEX
sex_dic[sex_dic==2]<-"FEMALE"
sex_dic[sex_dic==1]<-"MALE"
names(sex_dic)<-sex$SUBJID


## Genotype quality
head(genotypeq[1:5,1:5])
genoq_to_melt<-genotypeq[,-c(3)]
head(genoq_to_melt[1:5,1:5])
print("abt to melt")
genoq_melt<-data.table::melt.data.table(genoq_to_melt, id.vars =c("chrom","pos"))
head(genoq_melt)
colnames(genoq_melt)[3:4]<-c("ind","GQ")
print("colnamed")
data.table::fwrite(genoq_melt,file="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_melted.txt.gz", bom = TRUE)
##genoq_melt$GQ<-as.numeric(genoq_melt$GQ)
#genoq_melt$sex<-sex_dic[genoq_melt$ind]
genoq_sexdiff<-genoq_melt[,sex:=sex_dic[ind]]

print("melted...")
head(genoq_melt)
fwrite(genoq_sexdiff,file="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_long_wsex.txt.gz",bom = TRUE)

print("HIII")
head(genoq_melt[1:5,1:5])
##genoq_sexdiff<-dcast(genoq_melt,pos~sex)
print("Dcaseted")
####genoq_melt$GQ<-as.numeric(genoq_melt$GQ)
genoq_melt_numeric <- genoq_melt[, GQ:=as.numeric(GQ)]

genoq_melt_med<-(genoq_melt_numeric)[, as.list(summary(GQ)[1:6]), by = c("pos","sex")]
print("donezo")
head(genoq_melt_med)
data.table::fwrite(genoq_melt_med,file="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_melted_median.txt.gz", bom = TRUE)
genoq_sexdiff<-dcast((genoq_melt_med[,c("pos","sex","Median"),]),pos~sex)
genoq_sexdiff$sexdiff<-genoq_sexdiff$MALE-genoq_sexdiff$FEMALE
print("sex dfif")
head(genoq_sexdiff)
fwrite(genoq_sexdiff,file="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_long.txt.gz",bom = TRUE)
##system("gzip /oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_long.csv")
#######genoq_sexdiff<-dcast((genoq_melt[,c("pos","sex","Median"),]),pos~sex)

}
print("reading in")
genoq_long<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_melted_median.txt.gz")
print("read")
print(head(genoq_long))
ggplot(genoq_long%>%filter(Median<=30),aes(x=pos,y=Median,color=sex,alpha=0.8))+
  scale_color_manual(values=c("#FCD12A","#0492c2"))+ 
  geom_point()+
  ggtitle("GQ<=30: PAR1")+xlim(c(10000,2781479))
ggsave("/oak/stanford/groups/smontgom/raungar/Sex/Plots/GQ_xchr_lessthan30.png")
print("SAVED")



genoq_wsexdiff<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/qualitychecks/GQ_chrX_long.txt.gz")
ggplot(genoq_wsexdiff%>%filter(abs(sexdiff)>10),aes(x=pos,y=sexdiff,color=sexdiff))+
  scale_color_gradient2(low="#FCD12A",mid="gray",high="#0492c2")+
  geom_point()+ggtitle("sex diff (male-female) >10") #+xlim(c(10000,2781479))
ggsave("/oak/stanford/groups/smontgom/raungar/Sex/Plots/GQ_xchr_sexdiff.png")


