library(dplyr)
library(ggplot2)
library(data.table)
library('epitools')

betas=fread('/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8eqtl/alltissues_genes_sexDEGs_beta0.04954.txt.gz')
colnames(betas)=c('tissue','chr','Gene','beta')
betas_summ=betas%>%group_by(Gene)%>%summarise(max_beta=max(abs(beta)))

get_sig<-function(filepath,myfilt='sexDEG0111'){
  outliers=fread(filepath)
  
  outliers_w_beta=merge(outliers,betas_summ,all.x = T,by='Gene')
  outliers_w_beta=outliers_w_beta%>%mutate(sexDEG0111=ifelse(!is.na(max_beta)&max_beta>0.111,"sexDEG","non-sexDEG"))%>%
    mutate(sexDEG04954=ifelse(!is.na(max_beta),"sexDEG","non-sexDEG"))
  
  sexDEG_outlier=outliers_w_beta%>%dplyr::filter(get(myfilt)=='sexDEG'&Y=='outlier')
  sexDEG_nonoutlier=outliers_w_beta%>%dplyr::filter(get(myfilt)=='sexDEG'&Y!='outlier')
  nonsexDEG_outlier=outliers_w_beta%>%dplyr::filter(get(myfilt)!='sexDEG'&Y=='outlier')
  nonsexDEG_nonoutlier=outliers_w_beta%>%dplyr::filter(get(myfilt)!='sexDEG'&Y!='outlier')
  contigency_table=rbind(sexDEG=c(sexDEG_outlier%>%nrow(),sexDEG_nonoutlier%>%nrow()),nonsexDEG=c(nonsexDEG_outlier%>%nrow(),nonsexDEG_nonoutlier%>%nrow()))
  colnames(contigency_table)<-c('outlier','non-outlier')
  
  epi=epitab(contigency_table,method="riskratio")
  risk_ratio=(epi$tab[2,5:ncol(epi$tab)]) #lower prop of outliers are sexDEGs in both gorup
  fisher_test=fisher.test(contigency_table)
  return(list(risk_ratio,fisher_test))
}

#allboth
outliers_aut_allboth0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_allboth_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_aut_allboth0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_allboth_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_allboth0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_allboth_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_allboth0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_allboth_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*

#female
outliers_aut_f0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_f_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_aut_f0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_f_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_f0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_f_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_f0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_f_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*

#male
outliers_aut_m0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_m_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_aut_m0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_m_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_m0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_m_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_m0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_m_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*

#both
outliers_aut_both0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_both_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_aut_both0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_both_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_both0.111=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_both_maxoutliers3.txt.gz",'sexDEG0111') #outliers_*zthresh2.5*nphen3*maxoutliers3*
outliers_x_both0.04954=get_sig("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_both_maxoutliers3.txt.gz",'sexDEG04954') #outliers_*zthresh2.5*nphen3*maxoutliers3*

bind_info_rr<-function(my_list,my_chr,my_sex,my_sexdeg){
  cbind(data.frame(as.list(my_list[[1]])),chr=my_chr,sex=my_sex,sexdeg=my_sexdeg)
}

rr_df=rbind(bind_info_rr(outliers_aut_allboth0.04954,"aut","allboth",0.04954),
            bind_info_rr(outliers_aut_allboth0.111,"aut","allboth",0.111),
            bind_info_rr(outliers_x_allboth0.04954,"x","allboth",0.04954),
            bind_info_rr(outliers_x_allboth0.111,"x","allboth",0.111),
            
            bind_info_rr(outliers_aut_both0.04954,"aut","both",0.04954),
            bind_info_rr(outliers_aut_both0.111,"aut","both",0.111),
            bind_info_rr(outliers_x_both0.04954,"x","both",0.04954),
            bind_info_rr(outliers_x_both0.111,"x","both",0.111),
            
            bind_info_rr(outliers_aut_f0.04954,"aut","f",0.04954),
            bind_info_rr(outliers_aut_f0.111,"aut","f",0.111),
            bind_info_rr(outliers_x_f0.04954,"x","f",0.04954),
            bind_info_rr(outliers_x_f0.111,"x","f",0.111),
            
            bind_info_rr(outliers_aut_m0.04954,"aut","m",0.04954),
            bind_info_rr(outliers_aut_m0.111,"aut","m",0.111),
            bind_info_rr(outliers_x_m0.04954,"x","m",0.04954),
            bind_info_rr(outliers_x_m0.111,"x","m",0.111))%>%mutate(significance=ifelse(p.value<0.01,'sig','not-sig'))
            
ggplot(rr_df%>%dplyr::filter(sexdeg==0.04954),aes(x=sex,y=riskratio,group-chr,color=chr,shape=significance))+
  geom_point(position=position_dodge(width=0.7),size=6)+
  geom_hline(yintercept=1,color='red')+
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.7))+theme_bw()+
  ggtitle("RR of proportion of (non)sexDEG vs (non)outlier 0.04954")
