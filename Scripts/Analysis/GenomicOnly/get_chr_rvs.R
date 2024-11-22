library("dplyr")
library("mltools") #ecdf
library("data.table")
library(optparse)


option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--euro_file"), type = 'character', default = NULL, help = "path of european file (gtex_2017-06-05_v8_euro_VCFids.txt)"),
  make_option(c("--sex_file"), type = 'character', default = NULL, help = "path of sample annotations file that includes sex"),
  make_option(c("--chrlen_infile"), type = 'character', default = NULL, help = "path of file that has all chr lengths"),
  make_option(c("--out_numrvs"), type = 'character', default = NULL, help = "outfile: number of rare variants"),
  make_option(c("--out_sigdif"), type = 'character', default = NULL, help = "outfile: sig difference between m/f"),
  make_option(c("--chrtype"), type = 'character', default = NULL, help = "either x or autosomes")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
euro_file <- as.character(opt$euro_file)
sex_file <- as.character(opt$sex_file)
chrlens <- as.character(opt$chrlen_infile)
out_numrvs <- as.character(opt$out_numrvs)
out_sigdif <- as.character(opt$out_sigdif)
chrtype <- as.character(opt$chrtype)
# summary_stat <- as.character(opt$summary_stat)

print(paste0("my chr is : ",chrtype))
# if(chrtype != "aut" & chrtype != "x"){stop("ERROR: chrtype must be aut or x")}
# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/aut_all_rvs_inds_types.txt.gz"
# euro_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
# sex_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# out_numrvs<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_all_numrv.txt.gz"
# out_sigdif<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_all_sigdif.txt.gz"
# chrlens<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/GRCh38.p13_assembly_report.txt"
# chrtype="aut"


wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE", 
                                "GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")

#Read everything in and make proper hashes 
exp_data = fread(infile,data.table=F)
#exp_data=fread(cmd = paste("head -n 501 | tail -500", infile),data.table=F)
colnames(exp_data)<-c("chr","start","end","maf","gtex_sample","vartype", "ensg","genetype")
euro<-fread(euro_file)
euro_vec_unchecked<-as.character(data.frame(euro)[,1])
euro_vec<-euro_vec_unchecked[!(euro_vec_unchecked %in% wrong_self_reported_ancestry )]
exp_data_euro<-exp_data[exp_data$gtex_sample %in% euro_vec,]
sex_df<-fread(sex_file,data.table=F)
sex_hash<-sex_df$SEX
names(sex_hash)<-sex_df$SUBJID
chr_df<-fread(chrlens,data.table=F)
chr_hash<-chr_df$`Sequence-Length`
names(chr_hash)<-chr_df$`UCSC-style-name`
chr_hash["aut"]<-sum(chr_hash[1:22])
chr_hash["x"]<-chr_hash["chrX"]

#takes in cumulative MAF matrix and then across individuals for the same vartype/chr
#calculatates the mean, median, and sd at a given MAF and returns this df
#which has cols vartype,maf,chr,this_mean,this_median,this_sd
get_summary<-function(this_cum_maf, this_sdcol){
  
  cum_maf_sd<-as.data.table(this_cum_maf)[,
                                          Reduce(c,lapply(.SD,stats::sd)), 
                                          by=.(vartype,UpperBound,chr),.SDcols=c(this_sdcol)]
  cum_maf_summary<-as.data.table(this_cum_maf)[,
                                               Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                               by=.(vartype,UpperBound,chr),.SDcols=c(this_sdcol)]
  
  cum_maf_df<-merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr"), all=T)
  colnames(cum_maf_df)<-c("vartype","maf","chr","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
  upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
  lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
  cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
  return(as.data.frame(cum_maf_df_CI))
}

#mann-whitney-wilcoxon test across MAF on m vs. F for a given vartype
get_mww_test_m_f<-function(maf_m,maf_f,this_var,mult_test_type){
  mw_df<-data.frame(row.names = c("maf","pval","padj","vartype","chr"))
  for (this_chr in unique(maf_m$chr)){
  
    nrow_m<-maf_m %>% dplyr::filter(vartype==this_var & chr==this_chr) %>% nrow()
    nrow_f<-maf_f %>% dplyr::filter(vartype==this_var & chr==this_chr) %>% nrow()
    if(nrow_m == 0 | nrow_f ==0){
      print(paste0("not enough obs ", this_chr, " for ", this_var))
      this_mw_df<-do.call(rbind.data.frame,lapply(sort(unique(maf_m$UpperBound)),function(x){c(x,"NA","NA",this_var,this_chr)}))
      colnames(this_mw_df)<-c("maf","pval","padj","vartype","chr")
      mw_df<-rbind(mw_df,this_mw_df)
      #next
    }
    mw_test_p<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                   function(x)
                                     wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$UpperBound ==x  & maf_m$vartype==this_var &  maf_m$chr==this_chr,"N.cum.adj"]),
                                                 as.numeric(as.data.frame(maf_f)[maf_f$UpperBound ==x  &  maf_f$vartype==this_var &  maf_f$chr==this_chr,"N.cum.adj"]))$p.value)
    )
   this_mw_df<- data.frame("maf"=sort(unique(maf_m$UpperBound)), "pval"=mw_test_p,
                           "padj"=p.adjust(mw_test_p,n=nrow_f,method =mult_test_type), "vartype"=this_var,
                           "chr"=this_chr)
   mw_df<-rbind(mw_df,this_mw_df)
  }
  return(mw_df)
  
}

#uses mltools package to calculate the cumulative number of RVs at each MAF for
#each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
cum_maf<-as.data.table(exp_data_euro)[,
                                      Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
                                      by=.(gtex_sample,vartype,chr),.SDcols="maf"]

#subset to m/f/both
print("Head chr hash")
print(head(chr_hash[cum_maf$chr]))
cum_maf$N.cum.adj<-cum_maf$N.cum/chr_hash[cum_maf$chr]*10000
print("adjusted by length")

this_cum_maf_both<-cum_maf
this_cum_maf_m<-cum_maf[sex_hash[cum_maf$gtex_sample]==1,]
this_cum_maf_f<-cum_maf[sex_hash[cum_maf$gtex_sample]==2,]
print("Split into m/f")

#calculatates the mean, median, and sd at a given MAF and returns this df
#which has cols vartype,maf,chr,this_mean,this_median,this_sd
#cum_maf_summ_all<-get_summary(this_cum_maf_both,"N.cum.adj")
cum_maf_summ_m<-get_summary(this_cum_maf_m,"N.cum.adj")
cum_maf_summ_f<-get_summary(this_cum_maf_f,"N.cum.adj")
cum_maf_summ_combined<-rbind(cbind(cum_maf_summ_m, "sex"="male"),
                                 cbind(cum_maf_summ_f, "sex"="female"))

#only double males for the X!
if(chrtype == "x"){
  this_cum_maf_m_double<-this_cum_maf_m
  this_cum_maf_m_double$N.cum<-this_cum_maf_m_double$N.cum*2
  this_cum_maf_m_double$N.cum.adj<-this_cum_maf_m_double$N.cum.adj*2
  
  cum_maf_summ_m_double<-get_summary(this_cum_maf_m_double, "N.cum.adj")
  cum_maf_summ_combined<-rbind(cbind(cum_maf_summ_m_double, "sex"="male"),
                               cbind(cum_maf_summ_f, "sex"="female"))
  # cum_maf_summ_combined<-rbind(cbind(cum_maf_summ_m, "sex"="male"),
  #                              cbind(cum_maf_summ_f, "sex"="female"),
  #                              cbind(cum_maf_summ_m_double, "sex"="male_double"))
  print("chr x - males have been doubled")
  
}

##SAVE THIS
write.table(cum_maf_summ_combined,gzfile(out_numrvs), quote = F,sep="\t",row.names = F) 

### get p val and adj pval for m vs. f comparison
#only double males for the X!
if(chrtype == "x"){
  mmw_snps<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"SNPs","BH")
  #mmw_indels<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"indels","BH")
 # mmw_sv<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"SV","BH")
  mmw_all<-rbind(mmw_snps[-1,] ) #,
                 #mmw_indels[-1,],
                # mmw_sv[-1,])
  
}else{
  print("snps...")
  mmw_snps<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SNPs","BH")
  print("indels...")
  #mmw_indels<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"indels","BH")
  print("sv....")
 # mmw_sv<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SV","BH")
  print("combining....")
  mmw_all<-rbind(mmw_snps[-1,]) #,
               #  mmw_indels[-1,],
               #  mmw_sv[-1,])
  
}
print("writing....")
write.table(mmw_all,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 
print("DONE!")



