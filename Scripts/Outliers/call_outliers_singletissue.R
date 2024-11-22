library(data.table)
library(tidyverse)
require(optparse)
##-- Read command line arguments
option_list = list(
  make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
  make_option(c('--min_num_tiss'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
  make_option(c('--tissue'), type = 'character', help = 'this tiss'),
  make_option(c('--N.PHEN'), type = 'character', help = 'tnphen'),
  make_option(c('--chrtype'), type = 'character', help = 'either aut or x'),
  make_option(c("--gene_chr_table"),type="character", default=NULL, help= "gtf with just ensg and chr"),
  make_option(c("--outfile"),type="character", default=NULL, help= "output file path"),
  make_option(c("--outfile_top"),type="character", default=NULL, help= "output file path"),
  make_option(c('--ZTHRESH'), type = 'numeric', default = 2, help = 'threshold for abs(MedZ) outliers')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

min_num_tiss = opt$min_num_tiss
zthresh = opt$ZTHRESH
outfile_any=opt$outfile
nphen=opt$N.PHEN
mytissue=opt$tissue
tissue=str_replace_all(mytissue,"-","_")
outfile_top=opt$outfile_top
chrtype=opt$chrtype
chr_file=opt$gene_chr_table
zscore=opt$Z.SCORES

# 
# 
# zscore="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_x_f.txt.gz"
# globalFile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2_nphen5_globalOutliersRemoved_x_f.txt"
# # chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_gtf_protclinc_wchr.txt"
#  chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/x_proteincoding_lncrna.gtf"
# chrtype="x"
# nphen=3
#tissue="Adipose_Subcutaneous"
# zthresh=2.5

if(chrtype!="x"){
  chrs_table=fread(chr_file,header = F)
  chrs_dic<-as.character(chrs_table$V2)
  names(chrs_dic)<-as.character(chrs_table$V1)
}
print(paste0("ZSCORES: ", zscore))
##-- Analysis
call.outliers<-function(data, nphen, zthresh,metric) {
  data.melted<-melt.data.table(data,id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z') %>%dplyr::filter(!is.na(Z))
  colnames(data.melted)[3:4] = c('Ind', 'Z')
  this.test.stats=data.melted%>%mutate(Df=1) %>%mutate(Y=ifelse(abs(Z)>zthresh,"outlier","control"))
  return(this.test.stats)
}
remove_global_outliers<-function(medz_data){
  medz_outliers = filter(medz_data, Y == 'outlier')
  total_num_genes=medz_data %>%select(Gene) %>%unique() %>% nrow()
  onepercent_genes<-floor(total_num_genes*0.01)
  # Look at count of outliers per individual
  medzCount = as.data.frame(table(medz_outliers$Ind)) %>% rename("Ind"="Var1")
  
  # Get count of genes per individual
  outliers_per_ind = medz_data %>% select(Gene,Ind) %>% group_by(Ind)%>%mutate(Freq=n())
  medzCountAll = merge(medzCount,outliers_per_ind, by='Ind')%>% mutate(PropOut = Freq.x/Freq.y) %>% dplyr::filter(PropOut !=0)
  medzCountAll$PropOut<-replace(medzCountAll$PropOut,is.na(medzCountAll$PropOut),0)
  #medzCount=medzCount
  # Remove individuals with more outliers than Q3 + 1.5*IQ
  q3 = quantile(medzCountAll$PropOut)[4]
  q1 = quantile(medzCountAll$PropOut)[2]
  qthres = q3 + 1.5*(q3-q1)
  indsRemove = unique(filter(medzCountAll, PropOut > qthres & Freq.x>onepercent_genes)$Var1)
  
  medz_data = filter(medz_data, !(Ind %in% indsRemove)) %>% arrange(desc(abs(Z))) %>%
    relocate(Tissue, .after = last_col()) %>% relocate(Gene,.after=Ind) %>% 
    rename(ind=Ind)%>%rename(ensg=Gene)%>%rename(MedZ=Z)
  
}

## Read in the normalized data
data = (fread( zscore)) %>% dplyr::filter(Tissue==tissue)
outliers_medz = call.outliers(data, nphen, zthresh,metric = 'medz')
#outliers_medz_top = pick.outliers(outliers_medz, nphen, zthresh)
outliers<-outliers_medz
outliers_top<-outliers_medz### outlier is only top gene + zthresh
outliers_noglobal<-remove_global_outliers(outliers)
outliers_top_noglobal<-remove_global_outliers(outliers_top)
write.outliers <- function(outliers, filename) {
  write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
}

if(chrtype!="x"){
  chrs_table=fread(chr_file,header = F)
  chrs_dic<-chrs_table$V2
  names(chrs_dic)<-chrs_table$V1
  outliers_noglobal$chr<-chrs_dic[outliers_noglobal$Gene]
  outliers_top_noglobal$chr<-chrs_dic[outliers_top_noglobal$Gene]
}
if(chrtype=="x"){
  outliers_noglobal$chr<-"x"
  outliers_top_noglobal$chr<-"x"
}


write.outliers(outliers_noglobal, outfile_any)
write.outliers(outliers_top_noglobal, outfile_top)

