library("ggplot2")
library("dplyr")
library("mltools") #ecdf
library("data.table")



infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/all_rvs_inds_types.txt.gz"
euro_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
sex_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"

#Read everything in
exp_data = fread(infile,data.table=F)
colnames(exp_data)<-c("chr","start","end","maf","gtex_sample","vartype")
euro<-fread(euro_file)
euro_vec_unchecked<-as.character(data.frame(euro)[,1])
wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE", 
"GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")
euro_vec<-euro_vec_unchecked[!(euro_vec_unchecked %in% wrong_self_reported_ancestry )]
exp_data_euro<-exp_data[exp_data$gtex_sample %in% euro_vec,]

sex_df<-fread(sex_file,data.table=F)
sex_hash<-sex_df$SEX
names(sex_hash)<-sex_df$SUBJID
#temp<-rbind(head(exp_data,5000),head(tail(head(exp_data,100000),5000)),tail(head(exp_data,35000),5000))
#this uses data table (efficient for large dataset)
#uses mltools empirical_cdf function to get a CDF that also contains num cumulative
# such that it outputs the cumulative number of RV at a certain MAF
#such that that min_y=0 and max_y=(total #RV)

#takes in cumulative MAF matrix and then across individuals for the same vartype/chr
#calculatates the mean, median, and sd at a given MAF and returns this df
#which has cols vartype,maf,chr,this_mean,this_median,this_sd
get_summary<-function(this_cum_maf){

  cum_maf_sd<-as.data.table(this_cum_maf)[,
                                     Reduce(c,lapply(.SD,stats::sd)), 
                                     by=.(vartype,UpperBound,chr),.SDcols="N.cum"]
  cum_maf_summary<-as.data.table(this_cum_maf)[,
                                          Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                          by=.(vartype,UpperBound,chr),.SDcols="N.cum"]
  
  cum_maf_df<-merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr"), all=T)
  colnames(cum_maf_df)<-c("vartype","maf","chr","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
  upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
  lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
  cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
  return(as.data.frame(cum_maf_df_CI))
}

#same as get_summary but includes PAR as consideration
get_summary_par<-function(this_cum_maf, this_sdcol){
  
  cum_maf_sd<-as.data.table(this_cum_maf)[,
                                          Reduce(c,lapply(.SD,stats::sd)), 
                                          by=.(vartype,UpperBound,chr,par_region),.SDcols=this_sdcol]
  cum_maf_summary<-as.data.table(this_cum_maf)[,
                                               Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                               by=.(vartype,UpperBound,chr,par_region),.SDcols=this_sdcol]
  
  cum_maf_df<-merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr","par_region"), all=T)
  colnames(cum_maf_df)<-c("vartype","maf","chr","par","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
  upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
  lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
  cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
  return(as.data.frame(cum_maf_df_CI))
}
#mann-whitney-wilcoxon test across MAF on m vs. F for a given vartype
get_mww_test_m_f<-function(maf_m,maf_f,vartype,mult_test_type){
  mw_test_p<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                           function(x)
                             wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$UpperBound==x & maf_m$vartype==vartype,"N.cum"]),
                                         as.numeric(as.data.frame(maf_f)[maf_f$UpperBound==x & maf_f$vartype==vartype ,"N.cum"]))$p.value)
                    )
  mw.adj.p<-p.adjust(mw_test_p,method =mult_test_type)
  mw_df<-data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p,"padj"=mw.adj.p, "vartype"=vartype)
  return(mw_df)
  
}
#again, par considered
get_mww_test_m_f_par<-function(maf_m,maf_f,this_var,mult_test_type="BH"){
  print(this_var)
  mw_test_p_par1<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                           function(x)
                             wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR1" & UpperBound==x) %>% pull(N.cum)),
                                         as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="PAR1" & UpperBound==x) %>% pull(N.cum)))$p.value
  ))
  mw_test_p_nonpar<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                  function(x)
                                    wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="NONPAR" & UpperBound==x) %>% pull(N.cum)),
                                                         as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="NONPAR" & UpperBound==x) %>% pull(N.cum)))$p.value
  ))
  if(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR2") %>% nrow() ==0){
    mw_df<-rbind(
      data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par1,"padj"=p.adjust(mw_test_p_par1,method =mult_test_type), "vartype"=this_var, "par_region"="PAR1"),
      data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_nonpar,"padj"=p.adjust(mw_test_p_nonpar,method =mult_test_type), "vartype"=this_var, "par_region"="NONPAR")
    )
  }else{
    mw_test_p_par2<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                  function(x)
                                    wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR2" & UpperBound==x) %>% pull(N.cum)),
                                                as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="PAR2" & UpperBound==x) %>% pull(N.cum)))$p.value
    ))
    mw_df<-rbind(
      data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par1,"padj"=p.adjust(mw_test_p_par1,method =mult_test_type), "vartype"=this_var, "par_region"="PAR1"),
      data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_nonpar,"padj"=p.adjust(mw_test_p_nonpar,method =mult_test_type), "vartype"=this_var, "par_region"="NONPAR"),
      data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par2,"padj"=p.adjust(mw_test_p_par2,method =mult_test_type), "vartype"=this_var, "par_region"="PAR2")
      )
  }
               
  return(mw_df)
  
}

#uses mltools package to calculate the cumulative number of RVs at each MAF for
#each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
cum_maf<-as.data.table(exp_data_euro)[,
                    Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
                           by=.(gtex_sample,vartype,chr),.SDcols="maf"]

#subset to m/f/both
this_cum_maf_both<-cum_maf
this_cum_maf_m<-cum_maf[sex_hash[cum_maf$gtex_sample]==1,]
this_cum_maf_f<-cum_maf[sex_hash[cum_maf$gtex_sample]==2,]
this_cum_maf_f_half<-this_cum_maf_f
this_cum_maf_f_half$N.cum<-this_cum_maf_f_half$N.cum/2


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
cum_maf_summ_all<-get_summary_par(this_cum_maf_both)
cum_maf_summ_m<-get_summary_par(this_cum_maf_m)
cum_maf_summ_f<-get_summary_par(this_cum_maf_f)
cum_maf_summ_f_half<-get_summary(this_cum_maf_f_half)
#combine m+f for plotting
cum_maf_summ_plot_m_and_f<-rbind(cbind(cum_maf_summ_m, "sex"="male"),
                                 cbind(cum_maf_summ_f_half, "sex"="female"))

#split by variant types
cum_maf_summ_plot_m_and_f_snps<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="SNPs")
cum_maf_summ_plot_m_and_f_indels<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="indels")
cum_maf_summ_plot_m_and_f_sv<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="SV")

#and plot!!!
ggplot(cum_maf_summ_plot_m_and_f_snps,aes(x=maf,y=this_mean,color=sex))+
  geom_ribbon(aes(ymin=lower,ymax=upper, fill=sex), alpha=0.2, colour = NA)+
  ggtitle("Average Number of RVs: SNPs [95% CI]")+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean))

ggplot(cum_maf_summ_plot_m_and_f_snps,aes(x=maf,y=this_mean,color=sex))+
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  ggtitle("Average Number of RVs: SNPs [IQR]")+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean))



### get p val and adj pval for m vs. f comparison
mmw_snps<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f_half,"SNPs","BH")
mmw_indels<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f_half,"indels","BH")
mmw_sv<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f_half,"SV","BH")
mmw_all<-rbind(mmw_snps[-1,],
               mmw_indels[-1,],
               mmw_sv[-1,])

#plot p-valuex x vartypes
ggplot(mmw_all,aes(x=maf,y=-log(padj),group=vartype,color=vartype))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
  xlim(0,0.01)+
  ggtitle("Significance of Difference in Males vs. Females of Number of Rare Variants")+
  geom_line()


#PAR REGIONS
par1_s=10001; par1_e=2781479; par2_s=155701383; par2_e=156030895
len_par1=par1_e-par1_s;len_nonpar=par2_s-par1_e; len_par2=par2_e-par2_s
par_hash<-c(len_par1, len_nonpar,len_par2)
names(par_hash)<-c("PAR1","NONPAR","PAR2")

par_region=mapply(function(s,e){
  if(s>=par1_s & e<=par1_e){
    "PAR1"
    }else if(s>=par2_s & e<=par2_e){
    "PAR2"
    }else{
    "NONPAR"
  }
},exp_data_euro$start,exp_data_euro$end)
par_binary=mapply(function(s,e){
  if(s>=par1_s & e<=par1_e){
    "PAR"
  }else if(s>=par2_s & e<=par2_e){
    "PAR"
  }else{
    "NONPAR"
  }
},exp_data_euro$start,exp_data_euro$end)
exp_data_euro_wpars<-cbind(exp_data_euro,par_region,par_binary)
cum_maf_parreg<-as.data.table(exp_data_euro)[,
                                      Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
                                      by=.(gtex_sample,vartype,chr,par_region),.SDcols="maf"]

#adjust number per 1000 bp
cum_maf_parreg$N.cum.adj<-cum_maf_parreg$N.cum/par_hash[cum_maf_parreg$par_region]*10000
#subset to m/f/both
this_cum_maf_parreg_m<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==1,]
this_cum_maf_parreg_f<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==2,]
this_cum_maf_parreg_f_half<-this_cum_maf_parreg_f
this_cum_maf_parreg_f_half$N.cum<-this_cum_maf_parreg_f$N.cum/2
this_cum_maf_parreg_f_half$N.cum.adj<-this_cum_maf_parreg_f$N.cum.adj/2


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
cum_maf_summ_parreg_m<-get_summary_par(this_cum_maf_parreg_m, "N.cum")
cum_maf_summ_parreg_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum")
cum_maf_summ_parreg_f_half_par<-get_summary_par(this_cum_maf_parreg_f_half,"N.cum")


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
cum_maf_summ_parreg_adj_m<-get_summary_par(this_cum_maf_parreg_m,"N.cum.adj")
cum_maf_summ_parreg_adj_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum.adj")
cum_maf_summ_parreg_f_half_adj_par<-get_summary_par(this_cum_maf_parreg_f_half,"N.cum.adj")


#combine m+f for plotting
cum_maf_summ_plot_m_and_f_parreg<-rbind(cbind(cum_maf_summ_parreg_m, "sex"="male"),
                                 cbind(cum_maf_summ_parreg_f_half_par, "sex"="female"))
cum_maf_summ_plot_m_and_f_parreg_adj<-rbind(cbind(cum_maf_summ_parreg_adj_m, "sex"="male"),
                                        cbind(cum_maf_summ_parreg_f_half_adj_par, "sex"="female"))

#split by variant types
cum_maf_summ_plot_m_and_f_parreg_snps<-cum_maf_summ_plot_m_and_f_parreg %>% dplyr::filter(vartype=="SNPs")
cum_maf_summ_plot_m_and_f_parreg_indels<-cum_maf_summ_plot_m_and_f_parreg %>% dplyr::filter(vartype=="indels")
cum_maf_summ_plot_m_and_f_parreg_sv<-cum_maf_summ_plot_m_and_f_parreg %>% dplyr::filter(vartype=="SV")
#adj by region len
cum_maf_summ_plot_m_and_f_parreg_adj_snps<-cum_maf_summ_plot_m_and_f_parreg_adj %>% dplyr::filter(vartype=="SNPs")
cum_maf_summ_plot_m_and_f_parreg_adj_indels<-cum_maf_summ_plot_m_and_f_parreg_adj %>% dplyr::filter(vartype=="indels")
cum_maf_summ_plot_m_and_f_parreg_adj_sv<-cum_maf_summ_plot_m_and_f_parreg_adj %>% dplyr::filter(vartype=="SV")

##and plot with pars!
ggplot(cum_maf_summ_plot_m_and_f_parreg_adj_snps,aes(x=maf,y=this_mean,color=sex,group=interaction(sex,par)))+
  geom_ribbon(aes(ymin=lower,ymax=upper, fill=sex), alpha=0.2, colour = NA)+
  ggtitle("Average Number of RVs/10kb: SNPs [95% CI]")+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean,linetype=par))

ggplot(cum_maf_summ_plot_m_and_f_parreg_adj_snps,aes(x=maf,y=this_mean,color=sex,group=interaction(sex,par)))+
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  #xlim(0,0.09)+ylim(0,1000)+
  ggtitle("Average Number of RVs/10kb: SNPs [IQR]")+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean,linetype=par))



### get p val and adj pval for m vs. f comparison
mmw_snps_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m,this_cum_maf_parreg_f_half,"SNPs","BH")
mmw_indels_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m,this_cum_maf_parreg_f_half,this_var = "indels","BH")
mmw_sv_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m,this_cum_maf_parreg_f_half,"SV","BH")
mmw_all_par<-rbind(mmw_snps_par[-1,],
               mmw_indels_par[-1,],
               mmw_sv_par[-1,]) %>% dplyr::filter(maf!=0)

#plot p-valuex x vartypes
ggplot(mmw_all_par,aes(x=maf,y=-log(padj),color=vartype,group=interaction(vartype,par_region)))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
  xlim(0,0.01)+
  ggtitle("Significance of Difference in Males vs. Females of Number of Rare Variants")+
geom_line(aes(linetype=par_region))



