library(data.table)
library(tidyverse)
library(ggplot2)

infile_f_cadd0="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/OutliersAndRVs/x_collapsed_outliers_noglobal_medz_varAnnot_zthresh2.5_nphen3_f_CADDtypesGQ5BlacklistRemovedALL_CADD0_linc_prot_maxoutliers3_window10000.txt.gz"
infile_m_cadd0="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/OutliersAndRVs/x_collapsed_outliers_noglobal_medz_varAnnot_zthresh2.5_nphen3_m_CADDtypesGQ5BlacklistRemovedALL_CADD0_linc_prot_maxoutliers3_window10000.txt.gz"

infile_f_cadd15="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/OutliersAndRVs/x_collapsed_outliers_noglobal_medz_varAnnot_zthresh2.5_nphen3_f_CADDtypesGQ5BlacklistRemovedALL_CADD15_linc_prot_maxoutliers3_window10000.txt.gz"
infile_m_cadd15="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/OutliersAndRVs/x_collapsed_outliers_noglobal_medz_varAnnot_zthresh2.5_nphen3_m_CADDtypesGQ5BlacklistRemovedALL_CADD15_linc_prot_maxoutliers3_window10000.txt.gz"

subset_to_vartype<-function(myfile_f,myfile_m,vepfilter,cadd_filter){
  mydata_f<-fread(myfile_f)
  colnames(mydata_f)[21:22]<-c("veptype","vepconsq")
  filt_data_f<-mydata_f%>%dplyr::filter(grepl(vepfilter,veptype) & cadd_phred>=cadd_filter)
  filt_data_notvartype_f<-mydata_f%>%dplyr::filter(!grepl(vepfilter,veptype) & !is.na(veptype)& !is.na(sex))
  random_rows_f<-floor(runif(nrow(filt_data_f),1,nrow(filt_data_notvartype_f)))
  
  mydata_m<-fread(myfile_m)
  colnames(mydata_m)[21:22]<-c("veptype","vepconsq")
  filt_data_m<-mydata_m%>%dplyr::filter(grepl(vepfilter,veptype)& cadd_phred>=cadd_filter)
  filt_data_notvartype_m<-mydata_m%>%dplyr::filter(!grepl(vepfilter,veptype) & !is.na(veptype) & !is.na(sex))
  random_rows_m<-floor(runif(nrow(filt_data_m),1,nrow(filt_data_notvartype_m)))
  
  all_data<-rbind(cbind(filt_data_f,is_veptype=T),
                  cbind(filt_data_m,is_veptype=T),
                  cbind(filt_data_notvartype_m[random_rows_m,],is_veptype=F),
                  cbind(filt_data_notvartype_f[random_rows_f,],is_veptype=F))
  return(all_data)
}
veptype_filter="stop_gained"
cadd0_stopgained<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,veptype_filter,0)
cadd15_stopgained<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,veptype_filter,15)
cadd0_stop<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,veptype_filter,0)
cadd15_stop<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"stop",15)
cadd0_NMD<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,"NMD_transcript_variant",0)
cadd15_NMD<-subset_to_vartype(infile_f_cadd15,infile_m_cadd0,"NMD_transcript_variant",15)
cadd0_splice<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,"splice",0)
cadd15_splice<-subset_to_vartype(infile_f_cadd15,infile_m_cadd0,"splice",15)
cadd0_splicedonor<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,"splice_donor",0)
cadd15_splicedonor<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"splice_donor",15)
cadd0_spliceacceptor<-subset_to_vartype(infile_f_cadd0,infile_m_cadd0,"splice_acceptor",0)
cadd15_spliceacceptor<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"splice_acceptor",15)
cadd15_5UTR<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"5_prime_UTR_variant",15)
cadd15_3UTR<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"3_prime_UTR_variant",15)
cadd15_regregion<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"regulatory_region_variant",15)
cadd15_TFBS<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"TF_binding_site_variant",15)
cadd15_intron<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"intron_variant",15)
cadd15_upstream<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"upstream_gene_variant",15)
cadd15_downstream<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"downstream_gene_variant",15)
cadd15_syn<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"synonymous_variant",15)
cadd15_mis<-subset_to_vartype(infile_f_cadd15,infile_m_cadd15,"missense_variant",15)

veptype_filter="intr"
cadd_filter=15
ggplot(cadd15_mis,aes(x=MedZ,fill=is_veptype,alpha=0.7))+
  geom_density(bins=50)+theme_bw()+facet_wrap(~sex)+
  ggtitle(paste0(veptype_filter), paste0(" for cadd=",cadd_filter))




