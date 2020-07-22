library("dplyr")
library("mltools") #ecdf
library("data.table")



option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--euro_file"), type = 'character', default = NULL, help = "path of european file (gtex_2017-06-05_v8_euro_VCFids.txt)"),
  make_option(c("--sex_file"), type = 'character', default = NULL, help = "path of sample annotations file that includes sex"),
  make_option(c("--out_numrvs"), type = 'character', default = NULL, help = "outfile: number of rare variants"),
  make_option(c("--out_sigdif"), type = 'character', default = NULL, help = "outfile: sig difference between m/f")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
euro_file <- as.character(opt$euro_file)
sex_file <- as.character(opt$sex_file)
out_numrvs <- as.character(opt$out_numrvs)
out_sigdif <- as.character(opt$out_sigdif)

get_mww_test_m_f_par<-function(maf_m,maf_f,this_var,mult_test_type="BH",par_hash){
  print(this_var)
  mw_df<-data.frame(row.names = c("maf","pval","padj","vartype","par_region"))
  for (this_par_region in names(par_hash)){
    mw_test_p_reg<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                  function(x)
                                    wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region==this_par_region & UpperBound==x) %>% pull(N.cum)),
                                                as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region==this_par_region & UpperBound==x) %>% pull(N.cum)))$p.value
    ))
      mw_df_reg<- data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_reg,
                         "padj"=p.adjust(mw_test_p_reg,method =mult_test_type), "vartype"=this_var, 
                         "par_region"=this_par_region)
      mw_df<-rbind(mw_df,mw_df_reg)
  }
  
  return(mw_df)
  
}


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


subregion_df<-make_subregion_df()[[1]]
par_hash<-make_subregion_df()[[2]]

cum_maf_parreg<-data.frame(row.names=c("gtex_sample","vartype","chr","par_region","UpperBound","N.cum","CDF"))


#for this_subregion in rownames(subregion_df){
for (this_subregion in "XTR"){
  #first, get RVs in this subregion
  this_s<-subregion_df[this_subregion,"start"]
  this_e<-subregion_df[this_subregion,"end"]
  par_region=mapply(function(s,e){
    if(s>=this_s & e<=this_e){this_subregion}
    else{"NA"}
  },exp_data_euro$start,exp_data_euro$end)
  exp_data_euro_wpars<-cbind(exp_data_euro,par_region)
  
  #look only at this subregion
  exp_data_this_subregion<-exp_data_euro_wpars %>% dplyr::filter(par_region==this_subregion)
  #uses mltools package to calculate the cumulative number of RVs at each MAF for
  #each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
  cum_maf_parreg_this_subregion<-as.data.table(exp_data_this_subregion)[,
                                                   Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
                                                   by=.(gtex_sample,vartype,chr,par_region),.SDcols="maf"]
  cum_maf_parreg<-rbind(cum_maf_parreg,cum_maf_parreg_this_subregion)

}

#adjust number per 1000 bp

cum_maf_parreg$N.cum.adj<-cum_maf_parreg$N.cum/par_hash[cum_maf_parreg$par_region]*10000
#subset to m/f/both
this_cum_maf_parreg_m<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==1,]
this_cum_maf_parreg_f<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==2,]
this_cum_maf_parreg_m_double<-this_cum_maf_parreg_f
this_cum_maf_parreg_m_double$N.cum<-this_cum_maf_parreg_m_double$N.cum*2
this_cum_maf_parreg_m_double$N.cum.adj<-this_cum_maf_parreg_m_double$N.cum.adj*2


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
cum_maf_summ_parreg_m<-get_summary_par(this_cum_maf_parreg_m, "N.cum")
cum_maf_summ_parreg_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum")
cum_maf_summ_parreg_m_double<-get_summary_par(this_cum_maf_parreg_m_double,"N.cum")


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
cum_maf_summ_parreg_adj_m<-get_summary_par(this_cum_maf_parreg_m,"N.cum.adj")
cum_maf_summ_parreg_adj_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum.adj")
cum_maf_summ_parreg_m_double_adj_par<-get_summary_par(this_cum_maf_parreg_m_double,"N.cum.adj")

#combine m+f for plotting
cum_maf_summ_plot_m_and_f_parreg<-rbind(cbind(cum_maf_summ_parreg_m_double, "sex"="male"),
                                        cbind(cum_maf_summ_parreg_f, "sex"="female"))
cum_maf_summ_plot_m_and_f_parreg_adj<-rbind(cbind(cum_maf_summ_parreg_m_double_adj_par, "sex"="male"),
                                            cbind(cum_maf_summ_parreg_f_half_adj_par, "sex"="female"))
#and write!!!
write.table(cum_maf_summ_plot_m_and_f_parreg_adj,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 



###and significance:
mmw_snps_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m_double,this_cum_maf_parreg_f,"SNPs","BH",par_hash)
mmw_indels_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m_double,this_cum_maf_parreg_f,this_var = "indels","BH",par_hash)
mmw_sv_par<-get_mww_test_m_f_par(this_cum_maf_parreg_m_double,this_cum_maf_parreg_f,"SV","BH",par_hash)
mmw_all_par<-rbind(mmw_snps_par[-1,],
                   mmw_indels_par[-1,],
                   mmw_sv_par[-1,]) 

write.table(mmw_all,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 



