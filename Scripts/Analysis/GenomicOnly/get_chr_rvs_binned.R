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
print(paste0("chrtype is ",chrtype))
# summary_stat <- as.character(opt$summary_stat)
# 
# print(paste0("my chr is : ",chrtype))
# # if(chrtype != "aut" & chrtype != "x"){stop("ERROR: chrtype must be aut or x")}
# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/Combined/x_all_rvs_inds_typesALL_linc_prot.txt.gz"
# infile="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/Combined/chr7_all_rvs_inds_CADDtypesSeenTwice_linc_prot.txt.gz"
# euro_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
# sex_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# #out_numrvs<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/17_all_numrv.txt.gz"
# #out_sigdif<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_all_sigdif.txt.gz"
# chrlens<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/effective_chr_len_prot_linc.txt"
# chrtype="21"


wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE",
                                "GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")

#Read everything in and make proper hashes 
exp_data = fread(infile,data.table=F)
#exp_data=fread(cmd = paste("head -n 501 | tail -500", infile),data.table=F)
#if(chrtype=="x"){
  colnames(exp_data)<-c("chr","start","end","gtex_maf","gnomad_maf","maf","gtex_sample","vartype","ensg", "genetype","sex","cadd_raw","cadd_phred")
  head(exp_data)
  
#}else{
 # colnames(exp_data)<-c("chr","start","end","gtex_maf","gnomad_maf","maf","gtex_sample","vartype","ensg","genetype","ensg","genetype","sex","cadd_raw","cadd_phred")
#}
# # print("filter by: ","chr",toupper(chrtype))
# # exp_data<-exp_data %>% dplyr::filter(chr == paste0("chr",toupper(chrtype)))
# print("filtered chrtype")
head(exp_data)
euro<-fread(euro_file)
euro_vec_unchecked<-as.character(data.frame(euro)[,1])
euro_vec<-euro_vec_unchecked[!(euro_vec_unchecked %in% wrong_self_reported_ancestry )]
exp_data_euro<-exp_data[exp_data$gtex_sample %in% euro_vec,]
sex_df<-fread(sex_file,data.table=F)
sex_hash<-sex_df$SEX
names(sex_hash)<-sex_df$SUBJID
chr_df<-fread(chrlens,data.table=F)
colnames(chr_df)<-c("chr","effective_length")
chr_hash<-chr_df$effective_length
names(chr_hash)<-chr_df$chr
chr_hash["aut"]<-sum(chr_hash[1:22])
chr_hash["x"]<-chr_hash["chrX"]

exp_data_euro$maf<-as.numeric(exp_data_euro$maf)
print(head(exp_data_euro))

#takes in cumulative MAF matrix and then across individuals for the same vartype/chr
#calculatates the mean, median, and sd at a given MAF and returns this df
#which has cols vartype,maf,chr,this_mean,this_median,this_sd
get_summary<-function(this_cum_maf, this_sdcol){
  
  cum_maf_sd<-as.data.table(this_cum_maf)[,
                                          Reduce(c,lapply(.SD,stats::sd)), 
                                          by=.(vartype,UpperBound,chr,range_min, range_max),.SDcols=c(this_sdcol)]
  cum_maf_summary<-as.data.table(this_cum_maf)[,
                                               Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                               by=.(vartype,UpperBound,chr,range_min, range_max),.SDcols=c(this_sdcol)]
  
  cum_maf_df<-(merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr", "range_min","range_max"), all=T)) #[,-c("NA\'s")]
  colnames(cum_maf_df)<-c("vartype","maf","chr", "range_min","range_max","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
  upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
  lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
  cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
  return(as.data.frame(cum_maf_df_CI))
}

#mann-whitney-wilcoxon test across MAF on m vs. F for a given vartype
get_mww_test_m_f<-function(maf_m,maf_f,this_var,mult_test_type){
  mw_df<-data.frame(row.names = c("maf","pval","padj","vartype","chr", "range_min","range_max"))
  for (this_chr in unique(maf_m$chr)){
  
    nrow_m<-maf_m %>% dplyr::filter(vartype==this_var & chr==this_chr) %>% nrow()
    nrow_f<-maf_f %>% dplyr::filter(vartype==this_var & chr==this_chr) %>% nrow()
    if(nrow_m == 0 | nrow_f ==0){
      print(paste0("not enough obs ", this_chr, " for ", this_var))
      this_mw_df<-do.call(rbind.data.frame,lapply(sort(unique(maf_m$UpperBound)),function(x){c(x,"NA","NA",this_var,this_chr,maf_m$range_min,maf_m$range_max )}))
      colnames(this_mw_df)<-c("maf","pval","padj","vartype","chr", "range_min","range_max")
      mw_df<-rbind(mw_df,this_mw_df)
      next
    }
    mw_test_p<-unlist(lapply(sort(unique(maf_m$UpperBound)),
                                   function(x)
                                     wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$UpperBound ==x  & maf_m$vartype==this_var &  maf_m$chr==this_chr,"N.cum.adj"]),
                                                 as.numeric(as.data.frame(maf_f)[maf_f$UpperBound ==x  &  maf_f$vartype==this_var &  maf_f$chr==this_chr,"N.cum.adj"]))$p.value)
    )
   this_mw_df<- data.frame("maf"=sort(unique(maf_m$UpperBound)), "pval"=mw_test_p,
                           "padj"=p.adjust(mw_test_p,n=nrow_f,method =mult_test_type), "vartype"=this_var,
                           "chr"=this_chr, "range_min"=unique(maf_m$range_min), "range_max"=unique(maf_m$range_max))
   mw_df<-rbind(mw_df,this_mw_df)
  }
  return(mw_df)
  
}

#uses mltools package to calculate the cumulative number of RVs at each MAF for
#each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
# cum_maf_0<-as.data.table(exp_data_euro %>% dplyr::filter(maf==0))[,
#                                             Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=0)))),
#                                             by=.(gtex_sample,vartype,chr,ensg),.SDcols="maf"] %>% mutate(range_min=0,range_max=0)
cum_maf_0_0.001<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0 & maf <= 0.005))[,
                                                                                                Reduce(c,lapply(.SD, function(x) 
                                                                                                  as.list(empirical_cdf(x,ubounds=0.005)))),
                                                                                                by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0,range_max=0.001)
cum_maf_0.001_0.005<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0.001 & maf <= 0.005))[,
                                                                           Reduce(c,lapply(.SD, function(x) 
                                                                             as.list(empirical_cdf(x,ubounds=0.005)))),
                                                                           by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0.001,range_max=0.005)
cum_maf_0.005_0.01<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0.005 & maf <= 0.01))[,
                                                                                               Reduce(c,lapply(.SD, function(x) 
                                                                                                 as.list(empirical_cdf(x,ubounds=0.01)))),
                                                                                               by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0.005,range_max=0.01)
cum_maf_0.01_0.05<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0.01 & maf <= 0.05))[,
                                                                                               Reduce(c,lapply(.SD, function(x) 
                                                                                                 as.list(empirical_cdf(x,ubounds=0.05)))),
                                                                                               by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0.01,range_max=0.05)
cum_maf_0.05_0.1<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0.05 & maf <= 0.1))[,
                                                                                            Reduce(c,lapply(.SD, function(x) 
                                                                                              as.list(empirical_cdf(x,ubounds=0.1)))),
                                                                                            by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0.05,range_max=0.1)
cum_maf_0.1_0.25<-as.data.table(exp_data_euro %>% dplyr::filter(maf > 0.1 & maf <= 0.25))[,
                                                                                         Reduce(c,lapply(.SD, function(x) 
                                                                                           as.list(empirical_cdf(x,ubounds=0.25)))),
                                                                                         by=.(gtex_sample,vartype,chr),.SDcols="maf"] %>% mutate(range_min=0.1,range_max=0.25)

#subset to m/f/both
print("Head chr hash")
# print(head(chr_hash[cum_maf$chr]))

cum_maf_list<-list("cum_maf_0_0.001"=cum_maf_0_0.001,"cum_maf_0.001_0.005"=cum_maf_0.001_0.005,"cum_maf_0.005_0.01"=cum_maf_0.005_0.01,
                   "cum_maf_0.01_0.05"=cum_maf_0.01_0.05,"cum_maf_0.05_0.1"=cum_maf_0.05_0.1,"cum_maf_0.1_0.25"=cum_maf_0.1_0.25)
print("adjusted by length")
cum_maf_summ_combined_all<-data.frame(matrix(ncol = 15, nrow = 0))
colnames(cum_maf_summ_combined_all)<-c("vartype","maf","chr","range_min","range_max","this_min","IQR1",
            "this_median","this_mean","IQR3","this_max","this_sd","upper","lower","sex")
mmw_combined_all<-data.frame(matrix(ncol = 7, nrow = 0))
#mmw_all<-data.frame(matrix())
colnames(mmw_combined_all)<-c("maf","pval","padj","vartype","chr","range_min","range_max")
for(cum_maf in cum_maf_list){
  cum_maf$N.cum.adj<-cum_maf$N.cum/chr_hash[cum_maf$chr]*10000
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
  cum_maf_summ_combined_all<-rbind(cum_maf_summ_combined_all,cum_maf_summ_combined)

  ### get p val and adj pval for m vs. f comparison
  #only double males for the X!
  if(chrtype == "x"){
    print("X sig...")
    mmw_snps<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"SNPs","BH")
   # mmw_indels<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"indels","BH")
   # mmw_sv<-get_mww_test_m_f(this_cum_maf_m_double,this_cum_maf_f,"SV","BH")
    mmw_all<-rbind(mmw_snps )#,
                  # mmw_indels,
                  # mmw_sv)
    #
  }else{
    print("sig snps...")
    mmw_snps<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SNPs","BH")
    print(paste0("nrows snps: ",nrow(mmw_snps)))
    print("sig indels...")
   # mmw_indels<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"indels","BH")
    #print(paste0("nrows indels: ",nrow(mmw_indels)))
    print("sig sv....")
    #mmw_sv<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SV","BH")
    #print(paste0("nrows sv: ",nrow(mmw_sv)))
    print("combining....")
    #if(nrow(mmw_sv) != 0){
     #    mmw_all<-rbind(mmw_snps,
      #             mmw_indels,
       #            mmw_sv)
    #} else if(nrow(mmw_indels) == 0 ) {
   #   mmw_all<-rbind(mmw_all,mmw_snps)
    #}
    #else {
    #    mmw_all<-rbind(mmw_snps,
    #               mmw_indels)
    #}
  }

  #mmw_combined_all<-rbind(mmw_combined_all,mmw_all)
  mmw_combined_all<-rbind(mmw_combined_all,mmw_snps)

}

####FIX PVAL

#mmw_combined_all$padj<-p.adjust(mmw_combined_all$pval,n=length(cum_maf_list),method ="BH")

##SAVE THIS
write.table(cum_maf_summ_combined_all,gzfile(out_numrvs), quote = F,sep="\t",row.names = F) 
print("writing....")
write.table(mmw_combined_all,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 
print("DONE!")

