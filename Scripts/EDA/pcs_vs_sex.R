library("ggplot2")
library("data.table")

pc_file<-"/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt"
pcs<-fread(pc_file)

m_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_individuals_passed_m.txt"
f_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_individuals_passed_f.txt"
males<-read_csv(m_file,col_names = "COL")$COL
females<-read_csv(f_file,col_names = "COL")$COL
sexdic=rep("female",length(females))
names(sexdic)<-females
sexdic_m=rep("male",length(males))
names(sexdic_m)<-males
sexdic<-c(sexdic,sexdic_m)

pcs$sample=paste0("GTEX-",sapply(str_split((pcs$FID),"-"),"[[",c(2)))
pcs$sex=sexdic[pcs$sample]
ggplot(pcs,aes(x=PC5,y=PC6,color=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  geom_point()
