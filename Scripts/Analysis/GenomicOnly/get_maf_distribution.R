library("dplyr")
library("mltools") #ecdf
library("data.table")
library(optparse)


option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--euro_file"), type = 'character', default = NULL, help = "path of european file (gtex_2017-06-05_v8_euro_VCFids.txt)"),
  make_option(c("--sex_file"), type = 'character', default = NULL, help = "path of sample annotations file that includes sex"),
  make_option(c("--chrlen_infile"), type = 'character', default = NULL, help = "path of file that has all chr lengths"),
  make_option(c("--outfile"), type = 'character', default = NULL, help = "outfile: maf"),
  make_option(c("--chrtype"), type = 'character', default = NULL, help = "chrtype (aut or x)"),
  make_option(c("--maf"), type = 'numeric', default = NULL, help = "minor allele frequency to pull out")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
euro_file <- as.character(opt$euro_file)
sex_file <- as.character(opt$sex_file)
chrlens <- as.character(opt$chrlen_infile)
outfile <- as.character(opt$outfile)
chrtype <- as.character(opt$chrtype)
maf <- as.numeric(opt$maf)
# summary_stat <- as.character(opt$summary_stat)
print(paste0("my chr is : ",chrtype))
#if(chrtype != "aut" & chrtype != "x"){stop("ERROR: chrtype must be aut or x")}

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

final_maf_f<-this_cum_maf_f %>% dplyr::filter(UpperBound==maf) %>% dplyr::mutate("sex"="female")

if(chrtype == "x"){
  this_cum_maf_m_double<-this_cum_maf_m
  this_cum_maf_m_double$N.cum<-this_cum_maf_m_double$N.cum*2
  this_cum_maf_m_double$N.cum.adj<-this_cum_maf_m_double$N.cum.adj*2
  final_maf_m <- this_cum_maf_m_double %>% dplyr::filter(UpperBound==maf) %>% dplyr::mutate("sex"="male")
}else{
  final_maf_m<-this_cum_maf_f %>% dplyr::filter(UpperBound==maf) %>% dplyr::mutate("sex"="male")
  
}

final_maf_all<-rbind(final_maf_f,final_maf_m)


write.table(final_maf_all,gzfile(outfile), quote = F,sep="\t",row.names = F) 


