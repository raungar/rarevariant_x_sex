library(optparse)
library(glmnet)
library(data.table)
library(dplyr)
library(preprocessCore) #normalize.quantile
library(bestNormalize) #boxcox
library(coefplot)
library("ggplot2")


option_list = list(
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-k", "--sex_key"), type="character", default=NULL, help="sex_key file", metavar="character")
  
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sex_key_file<-as.character(opt$sex_key)
outfile<-as.character(opt$outfile)

 # sex_key_file="/Volumes/groups/smontgom/raungar/Sex/Output/nullshuffled_v8/sex_scrambled_key.txt"

metadata<-fread(sex_key_file,data.table=F)
colnames(metadata)<-c("SUBJID","SEX","NTISS","SexScrambled")
sex_dic<-metadata$SEX
names(sex_dic)<-metadata$SUBJID

### split into training and test
set.seed(90368)
#designed so already 50/50 m/f so if we are splitting 50/50, just need a random
random_sampling<-sample(1:nrow(metadata)) %% 2
test_train<-ifelse(random_sampling=="1","test","train") 
metadata$category<-test_train
metadata$NullSUBJID<-sample(metadata$SUBJID)
write.table(metadata,file=outfile,sep="\t",row.names=F,col.names=FALSE, quote=FALSE)

