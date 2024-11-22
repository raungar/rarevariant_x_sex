library("dplyr")
library("mltools") #ecdf
library("data.table")
library(optparse)



option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--euro_file"), type = 'character', default = NULL, help = "path of european file (gtex_2017-06-05_v8_euro_VCFids.txt)"),
  make_option(c("--sex_file"), type = 'character', default = NULL, help = "path of sample annotations file that includes sex"),
  make_option(c("--effective_chrlen"), type = 'character', default = NULL, help = "path of file that has all chr lengths"),
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
chrlens <- as.character(opt$effective_chrlen)

###chrs!
chr_df<-fread(chrlens,data.table=F)
colnames(chr_df)<-c("chr","effective_length")
chr_hash<-chr_df$effective_length
names(chr_hash)<-chr_df$chr
par_hash=chr_hash

# infile<-"/Volumes/groups/smontgom/raungar/Sex/Output/features_v8/Combined/x_all_rvs_inds_CADDtypesGQ5BlacklistRemovedALL_linc_prot.txt.gz"
# euro_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids_notambiguous.txt"
# chrlens<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/effective_chr_len_prot_linc_xsubtype.txt"
#  sex_file<-"/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
#  out_numrvs<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_subtypes_all_numrv.txt.gz"
# out_sigdif<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_subtypes_all_sigdif.txt"
# chrtype="x"


wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE", 
                                "GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")


exp_data = fread(infile,data.table=F)
#colnames(exp_data)<-c("chr","start","end","maf","gtex_sample","vartype","ensg","genetype")
colnames(exp_data)<-c("chr","start","end","maf_gtex","maf_gnomad","maf","gtex_sample","vartype","ensg","genetype","sex")
euro<-fread(euro_file)
euro_vec_unchecked<-as.character(data.frame(euro)[,1])
euro_vec<-euro_vec_unchecked[!(euro_vec_unchecked %in% wrong_self_reported_ancestry )]
exp_data_euro<-exp_data[exp_data$gtex_sample %in% euro_vec,]
sex_df<-fread(sex_file,data.table=F)
sex_hash<-sex_df$SEX
names(sex_hash)<-sex_df$SUBJID


get_mww_test_m_f_par<-function(maf_m,maf_f,this_var,mult_test_type="BH",par_hash){
     mw_df<-data.frame(row.names = c("maf","pval","padj","vartype","par_region"))
     for (this_par_region in names(par_hash)){
       nrow_m<-maf_m %>% dplyr::filter(vartype==this_var & par_region==this_par_region) %>% nrow()
       nrow_f<-maf_f %>% dplyr::filter(vartype==this_var & par_region==this_par_region) %>% nrow()
       if(nrow_m == 0 | nrow_f ==0){
         print(paste0("not enough obs ", this_par_region, " for ", this_var))
         #this_mw_df<-data.frame(maf=seq(0,.25, 0.001),"pval"="NA","padj"="NA", "vartype"=this_var, par_region=this_par_region)
         #mw_df<-rbind(mw_df,this_mw_df)
         next
       }
       mw_test_p_reg<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                     function(x)
                                       wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region==this_par_region & UpperBound==x) %>% pull(N.cum.adj)),
                                                   as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region==this_par_region & UpperBound==x) %>% pull(N.cum.adj)))$p.value
       ))
      # mw_test_p_reg<- wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region==this_par_region) %>% pull(N.cum)),
      #                                          as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region==this_par_region ) %>% pull(N.cum)))$p.value
      # mw_df_reg<- data.frame("pval"=mw_test_p_reg,
      #                      "padj"=p.adjust(mw_test_p_reg,n=nrow_f,method =mult_test_type), "vartype"=this_var,
      #                      "par_region"=this_par_region)

    mw_df_reg<- data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_reg,
                       "padj"=p.adjust(mw_test_p_reg,method =mult_test_type), "vartype"=this_var,
                       "par_region"=this_par_region)
    mw_df<-rbind(mw_df,mw_df_reg)
  }
  
  return(mw_df)
  
}


get_summary_par<-function(this_cum_maf, this_sdcol){
  print(paste0("THIS SD COL: ",this_sdcol))
  print(head(this_cum_maf))
  print(tail(this_cum_maf))
  if(nrow(this_cum_maf)==0){
	print("no observations")
	return(data.frame("vartype"=character(),"maf"=numeric(),"chr"=character(),
	"par"=character(),"this_min"=numeric(),"IQR1"=numeric(),
	"this_median"=numeric(),"this_mean"=numeric(),"IQR3"=numeric(),
	"this_max"=numeric(),"this_sd"=numeric()))
  }
  cum_maf_sd<-as.data.table(this_cum_maf)[,
                                          Reduce(c,lapply(.SD,stats::sd)), 
                                          by=.(vartype,UpperBound,chr,par_region),.SDcols=c(this_sdcol)]
  print("then summary")
  cum_maf_summary<-as.data.table(this_cum_maf)[,
                                               Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                               by=.(vartype,UpperBound,chr,par_region),.SDcols=c(this_sdcol)]
  
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

# par_hash<-make_subregion_df()[[2]]

cum_maf_parreg<-data.frame(row.names=c("gtex_sample","vartype","chr","par_region","UpperBound","N.cum","CDF", "N.cum.adj"))


for (this_subregion in rownames(subregion_df)){
  print(paste0("running: ",this_subregion))
#for (this_subregion in "XTR"){
  #first, get RVs in this subregion
  this_s<-subregion_df[this_subregion,"start"]
  this_e<-subregion_df[this_subregion,"end"]
  par_region=mapply(function(s,e){
    if(s>=this_s & e<=this_e){this_subregion}
    else{"NA"}
  },as.numeric(exp_data_euro$start),as.numeric(exp_data_euro$end))
  colnames(exp_data_euro)[c(12,13,14)]<-c("cadd_raw","cadd_phred","num_rv")
  exp_data_euro_wpars<-cbind(exp_data_euro,"par_region"=par_region)
  
  #look only at this subregion
  exp_data_this_subregion<-exp_data_euro_wpars %>% dplyr::filter(par_region==this_subregion)
  #uses mltools package to calculate the cumulative number of RVs at each MAF for
  #each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
  print("further in bin list")
  bin_list<-list(c(-1,0),c(0,0.001),c(0.001,0.005),c(0.005,0.01), c(0.01,0.05),c(0.05,0.1),c(0.1,0.25))
  print(bin_list)
  for (this_bin in bin_list){
    this_bin_min=this_bin[1]
    this_bin_max=this_bin[2]
    
  	print(paste0("this bin min: ",this_bin_min," and this bin max: ", this_bin_max))

    cum_maf_parreg_this_subregion<-as.data.table(exp_data_this_subregion %>% dplyr::filter(maf > this_bin_min & maf <= this_bin_max))[,
                                                     Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=this_bin_max)))),
                                                     by=.(gtex_sample,vartype,chr,par_region),.SDcols="maf"] %>% mutate(range_min=this_bin_min,range_max=this_bin_max)
    #adjust number per 10000 bp
    cum_maf_parreg_this_subregion$N.cum.adj<-cum_maf_parreg_this_subregion$N.cum/par_hash[this_subregion]*10000
    
    cum_maf_parreg<-rbind(cum_maf_parreg,cum_maf_parreg_this_subregion)
  }

}


###NA FOR SOME REASON FOR XAR CHECK IT OUT!
#subset to m/f/both
this_cum_maf_parreg_m<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==1,]
this_cum_maf_parreg_f<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==2,]
this_cum_maf_parreg_m_double<-this_cum_maf_parreg_m
#double for males
this_cum_maf_parreg_m_double$N.cum<-this_cum_maf_parreg_m_double$N.cum*2
this_cum_maf_parreg_m_double$N.cum.adj<-this_cum_maf_parreg_m_double$N.cum.adj*2

# print("cum_maf_parreg head and tail")
# print(head(cum_maf_parreg))
# print(tail(cum_maf_parreg))
# #get mean, median, and sd at a given MAF summary for indiviudals
# #at same vartype/maf/chr
# print("N.cum")
# cum_maf_summ_parreg_m<-get_summary_par(this_cum_maf_parreg_m, "N.cum")
# print("f")
# cum_maf_summ_parreg_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum")
# print("m.double")
# cum_maf_summ_parreg_m_double<-get_summary_par(this_cum_maf_parreg_m_double,"N.cum")


#get mean, median, and sd at a given MAF summary for indiviudals
#at same vartype/maf/chr
print("adj")
# cum_maf_summ_parreg_adj_m<-get_summary_par(this_cum_maf_parreg_m,"N.cum.adj")
# print("f")
# cum_maf_summ_parreg_adj_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum.adj")
# print("m.double")
# cum_maf_summ_parreg_m_double_adj_par<-get_summary_par(this_cum_maf_parreg_m_double,"N.cum.adj")
cum_maf_summ_plot_m_and_f_parreg_adj_all<-data.frame(matrix(nrow=0,ncol=16))
colnames(cum_maf_summ_plot_m_and_f_parreg_adj_all)<-c("vartype","maf","chr","par","this_min", "IQR1","this_median",
                                                      "this_mean","IQR3",  "this_max","this_sd","upper","lower",
                                                      "sex","bin_min","bin_max")
mmw_all_par<-data.frame(matrix(nrow=0,ncol=7))
colnames(mmw_all_par)<-c("maf","pval","padj","vartype","par_region","bin_min","bin_max")

# bin_list<-list(c(-1,0),c(0,0.00698),c(0.00698,0.00699),c(0.00699,0.001), c(0.001,0.01),c(0.01,0.05),c(0.05,0.1),c(0.1,0.25))
print("bin list")
print(bin_list)
for (this_bin in bin_list){
  this_bin_min=this_bin[1]
  this_bin_max=this_bin[2]
  print(paste0("doing: ",this_bin_min," - ",this_bin_max))
  binned_f<-this_cum_maf_parreg_f  %>% dplyr::filter(UpperBound > this_bin_min & UpperBound <= this_bin_max)
  binned_m<-this_cum_maf_parreg_m_double  %>% dplyr::filter(UpperBound > this_bin_min & UpperBound <= this_bin_max)
  print("getting summary")
  summ_parreg_adj_binned_f<-get_summary_par(binned_f,"N.cum.adj")
   print("summ summ_parreg_adj_binned_f: ")
   print(summ_parreg_adj_binned_f)
  summ_parreg_adj_binned_m<-get_summary_par(binned_m,"N.cum.adj")
  print("summ_parreg_adj_binned_m")
  print(summ_parreg_adj_binned_m)


  nrow_summ_parreg_adj_binned_m<-nrow(summ_parreg_adj_binned_m)
  nrow_summ_parreg_adj_binned_f<-nrow(summ_parreg_adj_binned_f)

  print(paste0("nrow m: ",nrow_summ_parreg_adj_binned_m,"  --- nrow f" ,nrow_summ_parreg_adj_binned_f))
  if (nrow_summ_parreg_adj_binned_m!=0 & nrow_summ_parreg_adj_binned_f!=0){
   print("neither zero")
   this_cum_maf_summ_plot_m_and_f_parreg_adj<-cbind(rbind(cbind(summ_parreg_adj_binned_m, "sex"="male"),
                                              cbind(summ_parreg_adj_binned_f, "sex"="female")),
                                              "bin_min"=this_bin_min,"bin_max"=this_bin_max)
	  print("combined.....")
	  print(this_cum_maf_summ_plot_m_and_f_parreg_adj)
	  cum_maf_summ_plot_m_and_f_parreg_adj_all<-rbind(cum_maf_summ_plot_m_and_f_parreg_adj_all,this_cum_maf_summ_plot_m_and_f_parreg_adj)

  } else if (nrow_summ_parreg_adj_binned_m==0 & nrow_summ_parreg_adj_binned_f==0){
	print("both were zero")
  } else if(nrow_summ_parreg_adj_binned_m==0){
     print("male zero")
	this_cum_maf_summ_plot_m_and_f_parreg_adj<-cbind(rbind(cbind(summ_parreg_adj_binned_f, "sex"="female"),
                                              "bin_min"=this_bin_min,"bin_max"=this_bin_max))
       print("combined.....")
       print(this_cum_maf_summ_plot_m_and_f_parreg_adj)
      cum_maf_summ_plot_m_and_f_parreg_adj_all<-rbind(cum_maf_summ_plot_m_and_f_parreg_adj_all,this_cum_maf_summ_plot_m_and_f_parreg_adj)
  } else if (nrow_summ_parreg_adj_binned_f==0){
      print("females 0")
	this_cum_maf_summ_plot_m_and_f_parreg_adj<-cbind(rbind(cbind(summ_parreg_adj_binned_m, "sex"="male"),
                                              "bin_min"=this_bin_min,"bin_max"=this_bin_max))	
       print("combined.....")
       print(this_cum_maf_summ_plot_m_and_f_parreg_adj)
      cum_maf_summ_plot_m_and_f_parreg_adj_all<-rbind(cum_maf_summ_plot_m_and_f_parreg_adj_all,this_cum_maf_summ_plot_m_and_f_parreg_adj)
	
  } else{
	print("shouldnt get here.. error")
  }

  
  mmw_snps_par<-get_mww_test_m_f_par(binned_m,binned_f,"SNPs","BH",par_hash)
  #mmw_indels_par<-get_mww_test_m_f_par(binned_m,binned_f,this_var = "indels","BH",par_hash)
  #mmw_sv_par<-get_mww_test_m_f_par(binned_m,binned_f,"SV","BH",par_hash)
  this_mmw_all_par<-rbind(mmw_snps_par ) #,
                     #mmw_indels_par,
                     #mmw_sv_par) 
  this_mmw_all_par_bins<-cbind(this_mmw_all_par,"bin_min"=rep(this_bin_min,nrow(this_mmw_all_par)),
                               "bin_max"=rep(this_bin_max,nrow(this_mmw_all_par)))
  
  mmw_all_par<-rbind(mmw_all_par,this_mmw_all_par_bins)
}


#combine m+f for plotting
# cum_maf_summ_plot_m_and_f_parreg<-rbind(cbind(cum_maf_summ_parreg_m_double, "sex"="male"),
#                                         cbind(cum_maf_summ_parreg_f, "sex"="female"))
# cum_maf_summ_plot_m_and_f_parreg_adj<-rbind(cbind(cum_maf_summ_parreg_m_double_adj_par, "sex"="male"),
#                                             cbind(cum_maf_summ_parreg_adj_f, "sex"="female"))
#and write!!!
write.table(cum_maf_summ_plot_m_and_f_parreg_adj_all,gzfile(out_numrvs), quote = F,sep="\t",row.names = F) 



###and significance:

write.table(mmw_all_par,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 



