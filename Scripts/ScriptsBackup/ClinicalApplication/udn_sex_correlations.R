library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(tidyverse)

dir<-"/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/eOutliers"

###TPM values
m_file<-paste0(dir,"/f_Blood_zscores_hg38_UDN.txt")
f_file<-paste0(dir,"/f_Blood_zscores_hg38_UDN.txt")
both_file<-paste0(dir,"/both_Blood_zscores_hg38_UDN.txt")

z_f<-fread(f_file)
z_m<-fread(m_file)
z_both<-fread(both_file)

metadata_file_udn<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/2019_12_05_Rare_Disease_Metadata.tsv")
affected=metadata_file_udn$affected_status
names(affected)<-metadata_file_udn$sample_id
sex_dic=metadata_file_udn$sex
names(sex_dic)<-metadata_file_udn$sample_id
z_f_info<-z_f%>%mutate(case_or_control=if_else(affected[sample_id]=="Case","Case","Control")) %>% mutate(sex="F")%>% mutate(subgroup="female")
z_m_info<-z_m%>%mutate(case_or_control=if_else(affected[sample_id]=="Case","Case","Control")) %>% mutate(sex="M")%>% mutate(subgroup="male")
z_both_info<-z_both%>%mutate(case_or_control=if_else(affected[sample_id]=="Case","Case","Control"))%>% mutate(sex=if_else(sex_dic[sample_id]=="F","F","M"))%>% mutate(subgroup="both")
z_all<-rbind(z_f_info,z_m_info,z_both_info)

z_all_spread<-z_all %>% 
  group_by(sample_id,gene,case_or_control,sex) %>%
  spread(key=subgroup,value=zscore)

z_all_spread$zdiff_m_both<-z_all_spread$male-z_all_spread$both
z_all_spread$zdiff_f_both<-z_all_spread$female-z_all_spread$both


ggplot(z_all_spread,aes(x=zdiff_m_both))+ggtitle("SexDEGs: Aut ZDiff M-Both") + 
  geom_density()
ggplot(z_all_spread,aes(x=(zdiff_m_both)))+ggtitle(" ZDiff M-Both") + 
  geom_histogram(bins = 100)
#
 z_f_highdiff=z_all_spread%>% dplyr::filter(abs(zdiff_f_both)>1)
 
 solved<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases.txt")#RD058
 solved_genes=z_all_spread %>% dplyr::filter(gene %in% solved$ENSG)
 
  