#!/usr/bin/env Rscript

## Script to call outliers given a file with Z-scores.
## The two first columns have the gene and tissue or phenotype, and the gene columne is named "Gene".
## Subsequent columns have Z-score data for each individual and the individual ID is the column name.

print("IN R SCRIPT")
## Get path to the GOATs data directories
#dir = Sys.getenv('RAREDIR')

## Load required packages
require(data.table)
require(optparse)
require(reshape2)
require(dplyr)
require(foreach)
require(doMC)
require(robustbase)

## Register parallel backend
registerDoMC(cores = 3)

#--------------- FUNCTIONS

#### Function to pick top outliers per gene
## Input: Data frame of outlier test statistics and the number of phenotypes required per test
## Output: Assignments of outliers and controls for genes with outliers
## Performs BF correction on the gene levels and Bh across genes
pick.outliers.old <- function(test.stats, nphen, zthresh,medz){
    test.stats = test.stats %>% filter(Df >= nphen) %>%
        group_by(Gene) %>%
        mutate(N = n())
        out = test.stats %>% arrange(desc(abs(MedZ))) %>%
            mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
            ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y,MaxBetaSex,MinBetaSex,AbsBetaSex)
    return(out)
}

#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers.old <- function(data, nphen=3, sexdeg,metric,zthresh=2) {
    data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
    colnames(data.melted)[3:4] = c('Ind', 'Z')
    data.melted_NoNA<-data.melted %>% dplyr::filter(!is.na(Z))
    data.melted.wsexdegs<-as.data.table(data.melted_NoNA)[sexdeg, on = c('Tissue','Gene'), BetaSex := beta_sex]
    data.melted.wsexdegs.inout<-data.melted.wsexdegs[,BetaSexOutlier:=ifelse(abs(Z)>zthresh,BetaSex,NA)]
    #filtering here will speed things up!
    data.melted.wsexdegs.inout.withN<-data.melted.wsexdegs[,Df:=.N,by=.(Ind, Gene)] %>% dplyr::filter(Df>=nphen)
    # head(data.melted)
    # test.stats = data.melted.wsexdegs %>% group_by(Ind, Gene) %>%
    #     summarise(MedZ = median(Z, na.rm = T),
    #               MaxBeta=max(BetaSex,na.rm=T),MinBeta=max(BetaSex,na.rm=T),AbsBeta=max(abs(BetaSex),na.rm=T),
    #               Df = sum(!is.na(Z)))
    # no_na_data<-na.omit(data.melted.wsexdegs,cols="Z")
    # test.stats=data.melted.wsexdegs.inout[, c("MedZ","Df","MaxBetaSex","MinBetaSex","AbsBetaSex"):=
    #                                    list(median(Z, na.rm = T),
    #                                         .N,
    #                                         max(BetaSexOutlier,na.rm=T),
    #                                         min(BetaSexOutlier,na.rm=T),
    #                                         max(abs(BetaSexOutlier),na.rm=T)),
    #                                  by=.(Ind, Gene)]
    #dat[, .(count = .N, var = sum(VAR)), by = MNTH]
    test.stats=(data.melted.wsexdegs.inout.withN)[,
                                             .(MedZ=median(Z, na.rm = T),
                                              Df=Df,
                                              MaxAbsBetaSex=max(abs(BetaSex),na.rm=T),
                                              MaxBetaSexOutlier=max(BetaSexOutlier,na.rm=T),
                                              MinBetaSexOutlier=min(BetaSexOutlier,na.rm=T),
                                              AbsBetaSexOutlier=max(abs(BetaSexOutlier),na.rm=T)),
                                             by=.(Ind, Gene)]
    # test.stats=data.melted.wsexdegs.inout[, .(MedZ=median(Z, na.rm = T),
    #                                         Df=.N,
    #                                         MaxBetaSexOutlier=max(BetaSexOutlier,na.rm=T),
    #                                         MinBetaSex=min(BetaSexOutlier,na.rm=T),
    #                                         AbsBetaSexmax(abs(BetaSexOutlier),na.rm=T)),
    #                                  by=.(Ind, Gene)]
    #don't report by tissue
     test.stats.noTiss<-unique(test.stats)
    
    # test.stats2=data.melted.wsexdegs[, .(MedZ=median(Z, na.rm = T)),
    #                                   by=.(Ind, Gene)]
    
    head(test.stats)
    outliers = pick.outliers(test.stats.noTiss, nphen, zthresh, medz = T)
    return(outliers)
}


#### Function to pick top outliers per gene
## Input: Data frame of outlier test statistics and the number of phenotypes required per test
## Output: Assignments of outliers and controls for genes with outliers
## Performs BF correction on the gene levels and Bh across genes
pick.outliers <- function(test.stats, nphen, zthresh,medz){
  this.test.stats = test.stats %>% filter(Df >= nphen) %>%
    group_by(Gene) %>%
    mutate(N = n())
  out = this.test.stats %>% arrange(desc(abs(MedZ))) %>%
    mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
    ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y,AbsBetaSex,MaxOutlierBetaSex,MinOutlierBetaSex,AbsOutlierBetaSex)
  out_topgene <- out %>% 
    group_by(Gene) %>% mutate(maxMedZ = max(abs(MedZ))) %>% 
    mutate(topOutlier=ifelse(abs(MedZ)==maxMedZ & abs(MedZ)>zthresh,"outlier","control")) %>%ungroup()%>% 
    select(Ind, Gene, N, Df, MedZ, topOutlier,AbsBetaSex,MaxOutlierBetaSex,MinOutlierBetaSex,AbsOutlierBetaSex)
  colnames(out_topgene)[6]<-"Y"
  return(list(out,out_topgene))
}

#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, nphen=3, zthresh,metric,sexdeg) {
  data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
  colnames(data.melted)[3:4] = c('Ind', 'Z')
  data.melted_NoNA<-data.melted %>% dplyr::filter(!is.na(Z))

  # data.melted.wsexdegs<-as.data.table(data.melted)[sexdeg, on = c('Tissue','Gene'), BetaSex := beta_sex]
  # data.melted.wsexdegs<-data.melted.wsexdegs[,BetaSexOutlier:=ifelse(abs(Z)>zthresh,BetaSex,NA)]

  data.melted.wsexdegs<-as.data.table(data.melted_NoNA)[sexdeg, on = c('Tissue','Gene'), BetaSex := beta_sex]
  data.melted.wsexdegs.inout<-data.melted.wsexdegs[,BetaSexOutlier:=ifelse(abs(Z)>zthresh,BetaSex,NA)]
  #filtering here will speed things up!
  data.melted.wsexdegs.inout.withN<-data.melted.wsexdegs[,Df:=.N,by=.(Ind, Gene)] %>% dplyr::filter(Df>=nphen)
  
  # head(data.melted)
  # test.stats = data.melted.wsexdegs %>% group_by(Ind, Gene) %>%
  #     summarise(MedZ = median(Z, na.rm = T),
  #               MaxBeta=max(BetaSex,na.rm=T),MinBeta=max(BetaSex,na.rm=T),AbsBeta=max(abs(BetaSex),na.rm=T),
  # #               Df = sum(!is.na(Z)))
  # no_na_data<-na.omit(data.melted.wsexdegs,cols="Z")
  # test.stats=no_na_data[, c("MedZ","Df","AbsBetaSex","MaxOutlierBetaSex","MinOutlierBetaSex","AbsOutlierBetaSex"):=
  #                         list(median(Z, na.rm = T),
  #                              .N,
  #                              max(abs(BetaSex),na.rm=T),
  #                              max(BetaSexOutlier,na.rm=T),
  #                              min(BetaSexOutlier,na.rm=T),
  #                              max(abs(BetaSexOutlier),na.rm=T)),
  #                       by=.(Ind, Gene)]
  test.stats=(data.melted.wsexdegs.inout.withN)[,
                                                .(MedZ=median(Z, na.rm = T),
                                                  Df=Df,
                                                  AbsBetaSex=max(abs(BetaSex),na.rm=T),
                                                  MaxOutlierBetaSex=max(BetaSexOutlier,na.rm=T),
                                                  MinOutlierBetaSex=min(BetaSexOutlier,na.rm=T),
                                                  AbsOutlierBetaSex=max(abs(BetaSexOutlier),na.rm=T)),
                                                by=.(Ind, Gene)]
  # test.stats.noTiss<-unique(test.stats[,`:=`(Tissue=NULL,Z=NULL)] )
   test.stats.uniq<-unique(test.stats)
  outliers = pick.outliers(test.stats.uniq, nphen, zthresh, medz = T)
  return(outliers)
}

##remove gloabl outliers
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

#### Function to write the outlier information to file.
## Input: data frame with data to write and output filename
## Output: None, writes directly to file.
write.outliers <- function(outliers, filename) {
    write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
}


#--------------- MAIN

##-- Read command line arguments
option_list = list(
  make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
  make_option(c('--N.PHEN'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
  make_option(c('--chrtype'), type = 'character', help = 'either aut or x'),
  make_option(c("--gene_chr_table"),type="character", default=NULL, help= "gtf with just ensg and chr"),
  make_option(c("--outfile"),type="character", default=NULL, help= "output file path"),
  make_option(c("--outfile_top"),type="character", default=NULL, help= "output file path"),
  make_option(c("--sexdeg_file"),type="character", default=NULL, help= "output file path"),
  make_option(c('--ZTHRESH'), type = 'numeric', default = 2, help = 'threshold for abs(MedZ) outliers')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

nphen = opt$N.PHEN
zthresh = opt$ZTHRESH
outfile_any=opt$outfile
outfile_top=opt$outfile_top
chrtype=opt$chrtype
chr_file=opt$gene_chr_table
zscore=opt$Z.SCORES
sexdeg_file=opt$sexdeg_file
# 
# zscore="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_x_f.txt.gz"
# globalFile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2_nphen5_globalOutliersRemoved_x_f.txt"
# # chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_gtf_protclinc_wchr.txt"
# sexdeg_file="/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/alltissues_genes_sexDEGs_beta0.111.txt.gz"
 #  chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/both_gtf_protclinc_wchr.txt"
 # chrtype="x"
 # nphen=3
 # zthresh=2.5


## Read in the normalized data
data = as.data.frame(fread(zscore))
sexdeg<-fread(sexdeg_file)
colnames(sexdeg)<-c("Tissue","chr","Gene","beta_sex")

if(chrtype!="x"){
  chrs_table=fread(chr_file,header = F)
  chrs_dic<-as.character(chrs_table$V2)
  names(chrs_dic)<-as.character(chrs_table$V1)
}
print(paste0("ZSCORES: ", zscore))
##-- Analysis


## Read in the normalized data
# data = (fread( zscore))
outliers_medz = call.outliers(data, nphen, zthresh,metric = 'medz',sexdeg)
#outliers_medz_top = pick.outliers(outliers_medz, nphen, zthresh)
outliers<-outliers_medz[[1]]
outliers_top<-outliers_medz[[2]] ### outlier is only top gene + zthresh


outliers_noglobal<-remove_global_outliers(outliers)
outliers_top_noglobal<-remove_global_outliers(outliers_top)
## Median Z-score
print('MEDZ')
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
#write.outliers(outliers.medz, paste0(dir, prefix, '.medz.txt'))
write.outliers(outliers_noglobal, outfile_any)
write.outliers(outliers_top_noglobal, outfile_top)