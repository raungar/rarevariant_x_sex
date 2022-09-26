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

print(paste0("my chr is : ",chrtype))
if(chrtype != "aut" & chrtype != "x"){stop("ERROR: chrtype must be aut or x")}

# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/x_all_rvs_inds_types.txt.gz"
# euro_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
# sex_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# out_numrvs<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_all_numrv.txt.gz"
# out_sigdif<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only/x_all_sigdif.txt.gz"
# chrlens<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/GRCh38.p13_assembly_report.txt"
# chrtype="x"


wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE", 
                                "GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")

#Read everything in and make proper hashes 
exp_data = fread(infile,data.table=F)
#exp_data=fread(cmd = paste("head -n 501 | tail -500", infile),data.table=F)
colnames(exp_data)<-c("chr","start","end","maf","gtex_sample","vartype")
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

#mann-whitney-wilcoxon test across MAF on m vs. F for a given vartype
get_mww_test_m_f<-function(maf_m,maf_f,this_var,mult_test_type){
  # print(paste0("vartype: ",this_var))
  # maf_df<-rbind(c(0,0.001),c(0.001,0.01),c(0.01,0.05),c(0.05,0.1),c(0.1,0.25))
  # rownames(maf_df)<-c("0-0.001","0.001-0.01","0.01-0.05","0.05-0.10","0.10-0.25")
  # colnames(maf_df)<-c("min","max")
  # mw_df<-data.frame(row.names = c("maf","pval","padj","vartype","chr"))
  # 
  # for (this_chr in unique(maf_m$chr)){
  #   nrow_m<-nrow(maf_m[ maf_m$vartype==this_var &  maf_m$chr==this_chr,"N.cum"])
  #   nrow_f<-nrow(maf_f[maf_f$vartype==this_var &  maf_f$chr==this_chr,"N.cum"])
  #   if(nrow_m == 0 | nrow_f ==0){
  #     print(paste0("not enough observation in ", this_chr, " for ", this_var))
  #     #next
  #   }
  #   for(maf_max in rownames(maf_df)){
  #     this_maf_min<-maf_df[maf_max,"min"]
  #     this_maf_max<-maf_df[maf_max,"max"]
  #     
  #     maf_m_filt<-maf_m %>% dplyr::filter(vartype ==this_var & chr==this_chr & UpperBound <= this_maf_max & UpperBound > this_maf_min) %>% as.data.frame()
  #     maf_f_filt<-maf_f %>% dplyr::filter(vartype ==this_var & chr==this_chr & UpperBound <= this_maf_max & UpperBound > this_maf_min) %>% as.data.frame()
  #     if(nrow(maf_m_filt) == 0 | nrow(maf_f_filt) ==0){
  #       print(paste0("skipping ", this_chr, " for ", vartype, " for MAF: ", maf_max))
  #       this_mw_df<-data.frame(maf=maf_max,"pval"="NA","padj"="NA", "vartype"=this_var, "chr"=this_chr)
  #       mw_df<-rbind(mw_df,this_mw_df)
  #       next
  #     }
  #     mw_test_p<-wilcox.test(as.numeric(maf_m_filt[,"N.cum"]), as.numeric(maf_f_filt[,"N.cum"]))$p.value
  #     
  #     mw.adj.p<-p.adjust(mw_test_p,n = length(maf_max), method =mult_test_type)
  #     this_mw_df<-data.frame(maf=maf_max,"pval"=mw_test_p,"padj"=mw.adj.p, "vartype"=this_var, "chr"=this_chr)
  #     mw_df<-rbind(mw_df,this_mw_df)
  #       # mw_test_p<-unlist(lapply(sort(unique(maf_m$UpperBound)),
  #       #                          function(x)
  #       #                            wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$vartype==vartype &  maf_m$chr==this_chr,"N.cum"]),
  #       #                                        as.numeric(as.data.frame(maf_f)[ maf_f$vartype==vartype &  maf_f$chr==this_chr,"N.cum"]))$p.value)
  #       # )
  #   }
  # }
  
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
                                     wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$UpperBound ==x  & maf_m$vartype==this_var &  maf_m$chr==this_chr,"N.cum"]),
                                                 as.numeric(as.data.frame(maf_f)[maf_f$UpperBound ==x  &  maf_f$vartype==this_var &  maf_f$chr==this_chr,"N.cum"]))$p.value)
    )
    mw_df<- data.frame("maf"=sort(unique(maf_m$UpperBound)), "pval"=mw_test_p,
                           "padj"=p.adjust(mw_test_p,n=nrow_f,method =mult_test_type), "vartype"=this_var,
                           "chr"=this_chr)
  }
  return(mw_df)
  
}

#uses mltools package to calculate the cumulative number of RVs at each MAF for
#each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
cum_maf<-as.data.table(exp_data_euro)[,
                                      Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
                                      by=.(gtex_sample,vartype,chr),.SDcols="maf"]

#subset to m/f/both
cum_maf$N.cum.adj<-cum_maf$N.cum/chr_hash[chrtype]*10000
print("n.cum.adj built")


this_cum_maf_both<-cum_maf
this_cum_maf_m<-cum_maf[sex_hash[cum_maf$gtex_sample]==1,]
this_cum_maf_f<-cum_maf[sex_hash[cum_maf$gtex_sample]==2,]
  mmw_snps<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SNPs","BH")
  mmw_indels<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"indels","BH")
  mmw_sv<-get_mww_test_m_f(this_cum_maf_m,this_cum_maf_f,"SV","BH")
  mmw_all<-rbind(mmw_snps[-1,],
                 mmw_indels[-1,],
                 mmw_sv[-1,])
  

write.table(mmw_all,gzfile(out_sigdif), quote = F,sep="\t",row.names = F) 




