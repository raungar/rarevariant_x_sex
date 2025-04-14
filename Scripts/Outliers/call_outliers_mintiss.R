#--------------- MAIN
## Load required packages
require(data.table)
require(optparse)
require(reshape2)
require(dplyr)
require(foreach)
require(doMC)
require(robustbase)

##-- Read command line arguments
option_list = list(
  make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
  make_option(c('--min_num_tiss'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
  make_option(c('--nphen'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
  make_option(c('--chrtype'), type = 'character', help = 'either aut or x'),
  make_option(c("--gene_chr_table"),type="character", default=NULL, help= "gtf with just ensg and chr"),
  make_option(c("--outfile"),type="character", default=NULL, help= "output file path"),
  make_option(c('--ZTHRESH'), type = 'numeric', default = 2, help = 'threshold for abs(MedZ) outliers')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

min_num_tiss = opt$min_num_tiss
nphen = opt$nphen
zthresh = opt$ZTHRESH
outfile_any=opt$outfile
chrtype=opt$chrtype
chr_file=opt$gene_chr_table
zscore=opt$Z.SCORES

# zscore="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_x_f.txt.gz"
# globalFile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2_min_num_tiss5_globalOutliersRemoved_x_f.txt"
# # chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_gtf_protclinc_wchr.txt"
#  chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/x_proteincoding_lncrna.gtf"
# chrtype="x"
# min_num_tiss=3
# zthresh=2.5
pick.outliers <- function(test.stats, nphen, zthresh,medz){
    this.test.stats = test.stats %>% filter(Df >= nphen) %>%
        group_by(Gene) %>%
        mutate(N = n())
        out = this.test.stats %>% arrange(desc(abs(MedZ))) %>%
            mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
            ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y)
        print("outliers all done")
        out_topgene <- out %>% 
          group_by(Gene) %>% mutate(maxMedZ = max(abs(MedZ))) %>% 
          mutate(topOutlier=ifelse(abs(MedZ)==maxMedZ & abs(MedZ)>zthresh,"outlier","control")) %>%ungroup()%>% 
        select(Ind, Gene, N, Df, MedZ, topOutlier)
        colnames(out_topgene)[6]<-"Y"
        print("outliers top done")
    return(list(out,out_topgene))
}

remove_global_outliers<-function(medz_data){
  medz_outliers = filter(medz_data, Y == 'outlier')
  total_num_genes=medz_data %>%select(Gene) %>%unique() %>% nrow() 
  onepercent_genes<-floor(total_num_genes*0.01)
  # Look at count of outliers per individual
  medzCount = as.data.frame(table(medz_outliers$Ind))
  
  # Get count of genes per individual
  outliers_per_ind = as.data.frame(table(medz_data$Ind))
  medzCountAll = merge(medzCount,outliers_per_ind, by='Var1')%>% mutate(PropOut = Freq.x/Freq.y)
  medzCountAll$PropOut<-replace(medzCountAll$PropOut,is.na(medzCountAll$PropOut),0)
  #medzCount=medzCount
  # Remove individuals with more outliers than Q3 + 1.5*IQ 
  q3 = quantile(medzCountAll$PropOut)[4] 
  q1 = quantile(medzCountAll$PropOut)[2]
  qthres = q3 + 1.5*(q3-q1)
  indsRemove = unique(filter(medzCountAll, PropOut > qthres & Freq.x>onepercent_genes)$Var1)
  
  medz_data = filter(medz_data, !(Ind %in% indsRemove))
  
}
#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, min_num_tiss, zthresh,nphen,metric) {
  # data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
  data.melted<-melt.data.table(data,id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
  colnames(data.melted)[3:4] = c('Ind', 'Z')
  print("melted data")
  #print(head(data.melted))

  my_med=data.melted[, MedZ := median(Z, na.rm = TRUE), by = list(Ind,Gene)]
  print("MedZ calculated")
  #test.stats=my_med%>%group_by(Ind,Gene)%>%mutate(Df=sum(!is.na(Z)))
  test.stats=my_med[,Df:=sum(!is.na(Z)),by=list(Ind,Gene)]
  print(head(test.stats))
  print('Df calculated')
  #test.stats=test.stats%>%select(-Z,-Tissue)
  # test.stats=test.stats[,c("Z","Tissue"):=NULL]
  colsToDelete=c("Z","Tissue")
  set(test.stats, , colsToDelete, NULL)
  print("remove columns")
  #test.stats=unique(test.stats)
  test.stats=(test.stats)%>%distinct()
  print('uniqued')
  test.stats=test.stats%>% mutate(Y=ifelse(MedZ>zthresh&Df>nphen,"outlier","control"))
  print('get outlier/cnotrol')
  test.stats=test.stats%>%arrange(desc(abs(MedZ)))
  print('arranged')
  test.stats=test.stats%>% group_by(Gene) %>%
    mutate(N = n())
    print('n=n()')

  print(head(test.stats%>%dplyr::filter(!is.na(MedZ))))


  # my.summs=function(x) list(MedZ=median(x,na.rm=T))
  # test.stats = head(data.melted,500)%>% group_by(Ind, Gene) %>%
  #   summarise(MedZ = median(Z, na.rm = T), 
  #              Df = sum(!is.na(Z)),
  #              MedZthreshpass=median(data.melted$Z[which(data.melted$Z>zthresh)]))%>%
  #   mutate(Y=ifelse(table(data.melted$Z>zthresh)["TRUE"]>min_num_tiss & Df>nphen,
  #                     "outlier","control"))%>%
  #   arrange(desc(abs(MedZ)))%>%
  #   group_by(Gene) %>%
  #   mutate(N = n()) %>% 
  #   select(Ind, Gene, N, Df, MedZ,MedZthreshpass, Y)
  # print("MedZ calculated")
  # print(head(test.stats))
  # stop("FOR NOW")
  # head(test.stats)
  outliers = pick.outliers(test.stats, min_num_tiss, zthresh, medz = T)
  return(outliers)
}

if(chrtype!="x"){
  chrs_table=fread(chr_file,header = F)
  chrs_dic<-as.character(chrs_table$V2)
  names(chrs_dic)<-as.character(chrs_table$V1)
}
print(paste0("ZSCORES: ", zscore))
##-- Analysis


## Read in the normalized data
data = (fread( zscore))
print("call outliers")
outliers_medz = call.outliers(data, min_num_tiss, zthresh,nphen,metric = 'mintis')
#outliers_medz_top = pick.outliers(outliers_medz, min_num_tiss, zthresh)
outliers<-outliers_medz[[1]]
outliers_top<-outliers_medz[[2]] ### outlier is only top gene + zthresh
print("remove global outliers")

outliers_noglobal<-remove_global_outliers(outliers)
print("write outliers")
write.outliers <- function(outliers, filename) {
    write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
}

write.outliers(outliers_noglobal, outfile_any)

