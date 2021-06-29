library(optparse)
library(glmnet)
library(data.table)
library(dplyr)
library(preprocessCore) #normalize.quantile
library(bestNormalize) #boxcox
library(coefplot)
library("ggplot2")


option_list = list(
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character"),
  make_option(c("-f", "--filter"), type="character", default=NULL, help="filter: aut, x, or both", metavar="character"),
  make_option(c("-p", "--outfile_preds"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-a", "--outfile_accuracy"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-g", "--genes_to_chrtype_file"), type="character", default=NULL, help="genes_to_chrtype_file file", metavar="character"),
  make_option(c("-c", "--counts_file"), type="character", default=NULL, help="counts control", metavar="character"),

) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
outfile_preds<-as.character(opt$outfile_preds)
outfile_accuracy<-as.character(opt$outfile_accuracy)
outdir<-as.character(opt$outdir)
region_filter<-as.character(opt$filter)
genes_to_chrtype_file<-as.character(opt$genes_to_chrtype_file)
tissue<-as.character(opt$tissue)
counts_file<-as.character(opt$counts_file)

### read in counts file
# counts_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/Whole_Blood.log2.ztrans.both.txt"
# genes_to_chrtype_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtf_padded10kb.bed"
# metadata_file="/Volumes/groups/smontgom/raungar/Sex/Files/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# tissue<-"Whole_Blood"
# region_filter<-"both"

metadata<-fread(metadata_file)
sex_dic<-metadata$SEX
names(sex_dic)<-metadata$SUBJID

# counts_f<-as.data.frame(t(data.frame(fread(counts_file_f),row.names=1)))
# counts_f$Id<-rownames(counts_f)
# counts_m<-as.data.frame(t(data.frame(fread(counts_file_m),row.names=1)))
# counts_m$Id<-rownames(counts_m)

#counts_all<-fread(counts_file)
counts_all<-as.data.frame(t(data.frame(fread(counts_file),row.names=1)))
counts_all$Id<-rownames(counts_all)



genes_to_chrtype<-fread(genes_to_chrtype_file)
genes_to_chrtype_dic<-genes_to_chrtype$V1
names(genes_to_chrtype_dic)<-genes_to_chrtype$V4
##filter by region x, aut, or both
if(region_filter=="x"){
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
  
  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
}else if(region_filter=="aut"){
  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  
}else if(region_filter=="both"){
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
 #  counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
 # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
}else{
  stop("ERROR: INCORRECT REGION FILTER. CHOOSE aut,both, or x ONLY IN LOWERCASE")
}

#only keep genes in both males and females
# overlapping_genes<-intersect(counts_subset_f$Id,counts_subset_m$Id)
# counts_f_touse<-counts_subset_f %>% dplyr::filter(Id %in% overlapping_genes)
# counts_m_touse<-counts_subset_m %>% dplyr::filter(Id %in% overlapping_genes)


### split into training and test
set.seed(90368)
#sample half to train, half to test for male and feamle
# sampling_f<-as.data.frame(counts_f_touse[,!(names(counts_f_touse) %in% "Id")])[,sample(1:(ncol(counts_f_touse)-1),ncol(counts_f_touse)-1)]
# rownames(sampling_f)<-counts_f_touse$Id
# train_f<-sampling_f[,1:ceiling(ncol(sampling_f)*.5)]
# test_f <- sampling_f[,(ceiling(ncol(sampling_f)*.5)+1):ncol(sampling_f)]
# 
# sampling_m<-as.data.frame(counts_m_touse[,!(names(counts_m_touse) %in% "Id")])[,sample(1:(ncol(counts_m_touse)-1),ncol(counts_m_touse)-1)]
# rownames(sampling_m)<-counts_m_touse$Id
# train_m<-sampling_m[,1:ceiling(ncol(sampling_m)*0.5)]
# test_m <- sampling_m[,(ceiling(ncol(sampling_m)*0.5)+1):ncol(sampling_m)]

sampling_all<-as.data.frame(counts_subset_all[,!(names(counts_subset_all) %in% "Id")])[,sample(1:(ncol(counts_subset_all)-1),ncol(counts_subset_all)-1)]
rownames(sampling_all)<-counts_subset_all$Id
train_all<-sampling_all[,1:ceiling(ncol(sampling_all)*.5)]
train_sex_lab<-sex_dic[colnames(train_all)]
test_all <- sampling_all[,(ceiling(ncol(sampling_all)*.5)+1):ncol(sampling_all)]
test_sex_lab<-sex_dic[colnames(test_all)]


#create actual training and test tests
# train_set <- cbind(train_f, train_m)
# train_sex_lab <- c(rep(0, ncol(train_f)), rep(1, ncol(train_m)))
# test_set <- cbind(test_f, test_m)
# test_sex_lab <- c(rep(0, ncol(test_f)), rep(1, ncol(test_m)))


# #convert to matrix and then then normalize
# x_train_mtrx<-apply(as.matrix(t(train_set)), c(1,2),as.numeric)
# #x_train<-apply(x_train_mtrx, 2, function(col) yeojohnson(col+0.5)$x.t) 
# x_train <- normalize.quantiles(x_train_mtrx)
# #x_train<-yeojohnson(x_train_mtrx)$x.t
# colnames(x_train)<-colnames(x_train_mtrx)
# rownames(x_train)<-rownames(x_train_mtrx)
# #x_test_mtrx<-apply(as.matrix(t(test_set)), c(1,2),as.numeric)
#  x_test <- normalize.quantiles(apply(as.matrix(t(test_set)), c(1,2),as.numeric))
# #x_test<-yeojohnson(x_test_mtrx)$x.t

#x_test<-apply(x_test_mtrx, 2, function(col) yeojohnson(col+0.5)$x.t) 
# colnames(x_test)<-colnames(x_test_mtrx)
# rownames(x_test)<-rownames(x_test_mtrx)
# data<- apply(as.matrix(t(cbind(sampling_f,sampling_m))),c(1,2),as.numeric)
# all_labs<-c(rep(0, ncol(sampling_f)), rep(1, ncol(sampling_m)))
# cvfit_all<-cv.glmnet(data,all_labs, 
#           family="binomial", 
#           alpha=my.alpha, nfolds = 6,
#           standardize=F)
# lambdas_to_try <- 10^seq(-40, 5, length.out = 100)
# cvfit_test<-glmnet(data,all_labs, lambda = lambdas_to_try,
#                      family="binomial", 
#                      alpha=0.5)
### train model
#https://github.com/erflynn/sl_label/blob/6e81605f4c336e0c1baa07abc168c72c9a9eaceb/code/sandbox/02_sex_labeling/02_train_enet_rnaseq.R
# full_matrix<-as.data.frame(rbind(x_train,x_test)) #as.data.frame(cbind("sex"=train_sex_lab,x_train))
# mylogit <- glm(c(train_sex_lab) ~ ., data=as.data.frame(t(train_all)), family = "binomial",maxit=100)
# mylogit_pred<- predict(mylogit,data=as.data.frame(test_all), type = "response")
# print(summary(mylogit))
list.alphas <- seq(0,1,0.1)
# my.alpha=.7
# cvfit = cv.glmnet(t(train_all), train_sex_lab, 
#                   family="binomial", 
#                   alpha=my.alpha, nfolds = 6,type.measure="class",
#                   standardize=F)
# preds_train <- predict(cvfit, newx=t(train_all), s="lambda.1se", type="response") #predict sex training
# preds_class_train <- sapply(predict(cvfit, newx=t(train_all), s="lambda.1se", type="class"), as.numeric)
# train_acc <- sum(preds_class_train==train_sex_lab)/length(train_sex_lab) #accuracy
# 
# preds_test <- predict(cvfit, newx=t(test_all), s="lambda.1se", type="response") #predict sex test
# preds_class_test <- sapply(predict(cvfit, newx=t(test_all), s="lambda.1se", type="class"), as.numeric)
# test_acc <- sum(preds_class_test==test_sex_lab)/length(test_sex_lab) 
# print(paste0("alpha=",my.alpha,", training accuracy=",train_acc,", test accuracy=",test_acc))
# coefs<-extract.coef(cvfit)
# coefs$chr<-genes_to_chrtype_dic[coefs$Coefficient]
# print(nrow(coefs))
# to_plot<-as.data.frame(preds_test)
# colnames(to_plot)<-"prob"
# to_plot$sex<-as.factor(sex_dic[rownames(to_plot)])

# ggplot(to_plot,aes(x=prob,fill=sex,alpha=0.9))+geom_density() +xlim(c(0,1))+ggtitle(paste0(region_filter,": alpha=",my.alpha))

hyperparam_res2 <- lapply(list.alphas, function(my.alpha){
  cvfit = cv.glmnet(t(train_all), train_sex_lab, 
                    family="binomial", 
                    alpha=my.alpha, nfolds = 6,type.measure="class",
                    standardize=F)
  
  preds_train <- predict(cvfit, newx=t(train_all), s="lambda.1se", type="response") #predict sex training
  preds_class_train <- sapply(predict(cvfit, newx=t(train_all), s="lambda.1se", type="class"), as.numeric)
  train_acc <- sum(preds_class_train==train_sex_lab)/length(train_sex_lab) #accuracy

  preds_test <- predict(cvfit, newx=t(test_all), s="lambda.1se", type="response") #predict sex test
  preds_class_test <- sapply(predict(cvfit, newx=t(test_all), s="lambda.1se", type="class"), as.numeric)
  test_acc <- sum(preds_class_test==test_sex_lab)/length(test_sex_lab) 
  print(paste0("alpha=",my.alpha,", training accuracy=",train_acc,", test accuracy=",test_acc))
  
  coefs<-extract.coef(cvfit)
  coefs$tissue<-tissue
  coefs$alpha<-my.alpha
  coefs$chr<-genes_to_chrtype_dic[coefs$Coefficient]
  this_coefs_file<-paste0(outdir,"/",tissue,"-coefs-",region_filter,"-alpha",my.alpha,".txt")
  #write.table(coefs,file=this_coefs_file,sep="\t",row.names=T,col.names=T, quote=FALSE)
  
  
  return(list("alpha"=my.alpha, "train_acc"=train_acc, "test_acc"=test_acc))
})




