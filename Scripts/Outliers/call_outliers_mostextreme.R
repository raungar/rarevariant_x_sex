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
library(readr)
require(robustbase)

## Register parallel backend
registerDoMC(cores = 3)

#--------------- FUNCTIONS

#### Function to pick top outliers per gene
## Input: Data frame of outlier test statistics and the number of phenotypes required per test
## Output: Assignments of outliers and controls for genes with outliers
## Performs BF correction on the gene levels and Bh across genes
pick.outliers <- function(test.stats, nphen, zthresh,medz){
    this.test.stats = test.stats %>% filter(Df >= nphen) %>%
        group_by(Gene) %>%
        mutate(N = n())
        print(.4)
        out = this.test.stats %>% arrange(desc(abs(MedZ))) %>%
            mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
            ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y)
        print(.5)
        out_topgene <- out %>% 
          group_by(Gene) %>% mutate(maxMedZ = max(abs(MedZ))) %>% 
          mutate(topOutlier=ifelse(abs(MedZ)==maxMedZ & abs(MedZ)>zthresh,"outlier","control")) %>%ungroup()%>% 
        select(Ind, Gene, N, Df, MedZ, topOutlier)
        print(.6)
        colnames(out_topgene)[6]<-"Y"
    return(list(out,out_topgene))
}

#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, nphen, zthresh,metric) {
     data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
    #data.melted<-melt.data.table(data,id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
    print(0.1)
    colnames(data.melted)[3:4] = c('Ind', 'Z')
    print(0.2)
    head(data.melted)
    # summ_z<-as.data.table(data.melted)[,Z:=.(median(Z,na.rm=T)),by=list(Ind,Gene)]
    # summ_df<-as.data.table(data.melted)[,Df:=.(sum(!is.na(Z))),by=list(Ind,Gene)]
    # test.stats<-cbind(summ_z,summ_df$Df)
    # my.summs=function(x) list(MedZ=median(x,na.rm=T))
    test.stats = as.data.table(data.melted) %>% group_by(Ind, Gene) %>%
        summarise(MedZ = median(Z, na.rm = T),
                  Df = sum(!is.na(Z)))
    print(0.3)

    # head(test.stats)
    outliers = pick.outliers(test.stats, nphen, zthresh, medz = T)
    print(0.9)
    return(outliers)
}

#### Function to write the outlier information to file.
## Input: data frame with data to write and output filename
## Output: None, writes directly to file.
write.outliers <- function(outliers, filename) {
    write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
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
#--------------- MAIN

##-- Read command line arguments
option_list = list(
	make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
	make_option(c('--N.PHEN'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier [num tissue]'),
	make_option(c('--chrtype'), type = 'character', help = 'either aut or x'),
	make_option(c("--gene_chr_table"),type="character", default=NULL, help= "gtf with just ensg and chr"),
	make_option(c("--outfile"),type="character", default=NULL, help= "output file path"),
	make_option(c("--outfile_top"),type="character", default=NULL, help= "output file path"),
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
# 
# zscore="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_x_f.txt.gz"
# globalFile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_zthresh2_nphen5_globalOutliersRemoved_x_f.txt"
# # chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_gtf_protclinc_wchr.txt"
#  chr_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/x_proteincoding_lncrna.gtf"
# chrtype="x"
# nphen=3
# z=2.5

if(chrtype!="x"){
  chrs_table=fread(chr_file,header = F)
  chrs_dic<-as.character(chrs_table$V2)
  names(chrs_dic)<-as.character(chrs_table$V1)
}
print(paste0("ZSCORES: ", zscore))
##-- Analysis


## Read in the normalized data
data = read_tsv(zscore)
print(head(data))
print(0)
outliers_medz = call.outliers(data, nphen, zthresh,metric = 'medz')
print(1)
#outliers_medz_top = pick.outliers(outliers_medz, nphen, zthresh)
outliers<-outliers_medz[[1]]
outliers_top<-outliers_medz[[2]] ### outlier is only top gene + zthresh

outliers_noglobal<-remove_global_outliers(outliers)
print(2)
outliers_top_noglobal<-remove_global_outliers(outliers_top)
print(3)
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
