library(ggplot2)
library(GenomicRanges)
library("RColorBrewer")

#read everything in
output_dir<-"/oak/stanford/groups/smontgom/raungar/Sex/Output"
mapq<-fread(paste0(output_dir,"/analysis_v8/genomic_only/qualitychecks/chrX_MQ.txt"))
colnames(mapq)<-c("pos","MQ")
genotypeq<-fread(paste0(output_dir,"/analysis_v8/genomic_only/qualitychecks/chrX_GQ_MAF.txt"),fill=T)
colnames(genotypeq)[1:3]<-c("chrom","pos","MAF")
sex<-fread("/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt")
sex_dic<-sex$SEX
sex_dic[sex_dic==2]<-"FEMALE"
sex_dic[sex_dic==1]<-"MALE"
names(sex_dic)<-sex$SUBJID

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

#make into genomic ranges
subregion_df_gr<-makeGRangesFromDataFrame(subregion_df,ignore.strand = T,keep.extra.columns = T)
mapq_gr<-makeGRangesFromDataFrame(cbind(mapq,chr="chrX",start=mapq$pos,end=mapq$pos),ignore.strand = T,keep.extra.columns = T)
mapq_subregion_overlap_gr<-findOverlaps(mapq_gr,subregion_df_gr,select = "first")
# mapq_region<-subregion_df_gr[(subjectHits(mapq_subregion_overlap_gr)),]$region
mapq$subregion<-subregion_df_gr[mapq_subregion_overlap_gr,]$region

ggplot(mapq %>% dplyr::filter(MQ<30),(aes(x=pos/1000000,y=MQ,color=subregion)))+geom_point()

low_qual<-mapq %>% dplyr::filter(MQ<30)
low_qual_par<-low_qual %>%dplyr::filter(subregion=="PAR1" | subregion=="PAR2" )%>% nrow()
low_qual_nonpar<-low_qual %>%dplyr::filter(subregion=="NONPAR" | subregion=="XTR" )%>% nrow()
total_par<-mapq %>%dplyr::filter(subregion=="PAR1" | subregion=="PAR2" )%>% nrow()
total_nonpar<-mapq  %>%dplyr::filter(subregion=="NONPAR"| subregion=="XTR" ) %>% nrow()

low_qual_par1_percent<-(low_qual %>%dplyr::filter(subregion=="PAR1" )%>% nrow())/(mapq  %>%dplyr::filter(subregion=="PAR1" )%>% nrow())*100
low_qual_par2_percent<-(low_qual %>%dplyr::filter(subregion=="PAR2"  )%>% nrow())/(mapq  %>%dplyr::filter(subregion=="PAR2" )%>% nrow())*100
#low_qual_xtr_percent<-(low_qual  %>%dplyr::filter( subregion=="XTR" ) %>% nrow())/(mapq  %>%dplyr::filter(subregion=="XTR" )%>% nrow())*100

print(paste0("So ",round(low_qual_par/total_par,3)*100, "% of PAR regions have low MapQ (MQ<30) reads, and ", round(low_qual_nonpar/total_nonpar,3)*100, "% of nonPAR regions have low MAPQ reads"))
print(paste0("For low MapQ (MQ<30) reads, this includes ",round(low_qual_par1_percent,3), "% of PAR1 regions, and ", 
             round(low_qual_par2_percent,3), "% of PAR2 regions." ))



## Genotype quality
genoq_melt<-melt(genotypeq, id.vars =c("chrom","pos"))
colnames(genoq_melt)[3:4]<-c("ind","GQ")
genoq_melt$GQ<-as.numeric(genoq_melt$GQ)
genoq_melt$sex<-sex_dic[genoq_melt$ind]
tmp<-genoq_melt[runif(10000, 1, nrow(genoq_melt)),]
#genoq_melt_med<-setDT(tmp)[, .(meanage = summary(GQ)), by = c("pos","sex")]
# genoq_melt_med<-setDT(tmp)[, c("min","IQR1","med","ave","IQR3","max") := (summary(GQ)), by = c("pos","sex")]
 genoq_melt_med<-(genoq_melt)[, as.list(summary(GQ)[1:6]), by = c("pos","sex")]
genoq_sexdiff<-dcast((genoq_melt_med[,c("pos","sex","Median"),]),pos~sex)
genoq_sexdiff$sexdiff<-genoq_sexdiff$MALE-genoq_sexdiff$FEMALE
ggplot(genoq_melt_med%>%filter(Median<=30),aes(x=pos,y=Median,color=sex,alpha=0.8))+
  scale_color_manual(values=c("#FCD12A","#0492c2"))+ 
  geom_point()+
  ggtitle("GQ<=30: PAR1")+xlim(c(10000,2781479))
ggplot(genoq_sexdiff%>%filter(abs(sexdiff)>10),aes(x=pos,y=sexdiff,color=sexdiff))+
  scale_color_gradient2(low="#FCD12A",mid="gray",high="#0492c2")+
  geom_point()+ggtitle("sex diff (male-female) >10") #+xlim(c(10000,2781479))


