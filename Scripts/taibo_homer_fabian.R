library(data.table)
library(tidyverse)
library(ggrepel)


file_tissue_stats="/oak/stanford/groups/smontgom/raungar/Sex/Files/SexBiased_GTEx_tissuestatistics.csv"
file_fabian="/oak/stanford/groups/smontgom/raungar/Sex/Files/TFTargetScores_SummarizedTFbySex_results.csv"

fabian=fread(file_fabian)
ggplot(fabian,aes(x=FemaleMean,y=MaleMean,alpha=0.8,label=TF))+
  geom_abline(slope=1,intercept=0,color="gray")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept = 0,color="gray")+
  geom_point()+
  theme_bw()+
  geom_text_repel()+
  xlab("Female mean FABIAN score x Watershed score")+
  ylab("Female mean FABIAN score x Watershed score")

ggplot(fabian,aes(x=FemaleMean,y=MaleMean,alpha=0.8,label=TF))+
  geom_abline(slope=1,intercept=0,color="gray")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept = 0,color="gray")+
  geom_point()+
  theme_bw()+
  geom_text_repel()+
  xlab("Female mean FABIAN score x Watershed score")+
  ylab("Female mean FABIAN score x Watershed score")

ggplot(fabian,aes(x=FemaleMax,y=MaleMax,alpha=0.8,label=TF))+
  geom_abline(slope=1,intercept=0,color="gray")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept = 0,color="gray")+
  geom_point()+
  theme_bw()+
  geom_text_repel()+
  xlab("Female max FABIAN score x Watershed score")+
  ylab("Male max FABIAN score x Watershed score")

