library(optparse)
library(glmnet)
library(data.table)
library(dplyr)
library(preprocessCore) #normalize.quantile
library(caret)
library(coefplot)
<<<<<<< HEAD
library(tidyverse)
=======
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
library("ggplot2")

option_list = list(
  make_option(c("-s", "--sextype"), type="character", default=NULL, help="used as colname for dictionary file. choose sextype=SexScrambled or SEX", metavar="character"),
  make_option(c("-f", "--filter"), type="character", default=NULL, help="filter: aut, x, or both", metavar="character"),
  make_option(c("-p", "--outfile_preds"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-a", "--outfile_accuracy"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-k", "--sex_key"), type="character", default=NULL, help="sex_key file", metavar="character"),
  make_option(c("-g", "--genes_to_chrtype_file"), type="character", default=NULL, help="genes_to_chrtype_file file", metavar="character"),
  make_option(c("-c", "--counts_file"), type="character", default=NULL, help="counts control", metavar="character"),
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character")

) 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
outfile_preds<-as.character(opt$outfile_preds)
outfile_accuracy<-as.character(opt$outfile_accuracy)
outdir<-as.character(opt$outdir)
region_filter<-as.character(opt$filter)
genes_to_chrtype_file<-as.character(opt$genes_to_chrtype_file)
tissue<-as.character(opt$tissue)
sextype<-as.character(opt$sextype) #used as colname for dictionary file. choose sextype=SexScrambled or SEX
counts_file<-as.character(opt$counts_file)
sex_key_file<-as.character(opt$sex_key)
# 
# ## read in counts file
<<<<<<< HEAD
# counts_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8eqtl/PEER_v8/Stomach.both_half.peer.v8ciseQTL.ztrans.txt"
=======
# counts_file="/Volumes/groups/smontgom/raungar/Sex/Output/nullshuffled_v8/PEER_v8/Artery_Aorta.log2.ztrans.both.txt"
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
# genes_to_chrtype_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtf_padded10kb.bed"
# sex_key_file="/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/Assignments/test_train.txt"
# metadata_file="/Volumes/groups/smontgom/raungar/Sex/Files/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# tissue<-"Artery_Aorta"
# region_filter<-"both"
<<<<<<< HEAD
#sex_key_file=metadata_file
=======
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
# dir="tmp"
# sextype="SexScrambled"
if(sextype != "SEX" & sextype != "SexScrambled"){
  stop("ERROR sextype must be either: SEX or SexScrambled")
}
print(paste0("filtering for: ",sextype))

metadata<-fread(sex_key_file)
if(sextype == "SexScrambled"){
  colnames(metadata)<-c("SUBJID","SEX","NTISS","SexScrambled","TestTrainStatus","NullSUBJID")
}
sex_dic<-metadata %>% pull(sextype)
names(sex_dic)<-as.character(metadata$SUBJID)
<<<<<<< HEAD
head(sex_dic)
age_dic<-metadata %>% pull(AGE)
names(age_dic)<-as.character(metadata$SUBJID)
=======
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
# group_dic<-metadata$TestTrainStatus
# names(group_dic)<-metadata$SUBJID
# counts_f<-as.data.frame(t(data.frame(fread(counts_file_f),row.names=1)))
# counts_f$Id<-rownames(counts_f)
# counts_m<-as.data.frame(t(data.frame(fread(counts_file_m),row.names=1)))
# counts_m$Id<-rownames(counts_m)

#counts_all<-fread(counts_file)
# tmp<-fread(counts_file)
# counts_f<-fread("/Volumes/groups/smontgom/raungar/Sex/Output/nullshuffled_v8/PEER_v8/Artery_Aorta.log2.ztrans.f.txt",data.table=F)
# counts_m<-fread("/Volumes/groups/smontgom/raungar/Sex/Output/nullshuffled_v8/PEER_v8/Artery_Aorta.log2.ztrans.m.txt",data.table=F)
<<<<<<< HEAD
counts_all<-as.data.frame((data.frame(fread(counts_file),row.names=1)))
colnames(counts_all) <-  str_replace(colnames(counts_all),"\\.","-")
=======
counts_all<-as.data.frame(t(data.frame(fread(counts_file),row.names=1)))
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
# counts_all<-as.data.frame(t(data.frame(fread(counts_file),row.names=1)))
counts_all$Id<-rownames(counts_all)


<<<<<<< HEAD
genes_to_chrtype<-fread(genes_to_chrtype_file)
genes_to_chrtype_dic<-genes_to_chrtype$V1
names(genes_to_chrtype_dic)<-genes_to_chrtype$V4
print("genes_to_chrtype_dic")
print(head(genes_to_chrtype_dic))
print(paste0("filtering to ",region_filter))
print(head(counts_all$Id))

print("That was head of counts all id")
##filter by region x, aut, or both
print((counts_all)[1:5,1:6])
if(region_filter=="x"){
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)]=="chrX")
  # counts_subset_all<-counts_all[,genes_to_chrtype_dic[colnames(counts_all)] =="chrX"]
  
  print("HI I AM IN X FILTER BABE")
  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
}else if(region_filter=="aut"){

  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22)))
  # counts_subset_all<-counts_all[,genes_to_chrtype_dic[colnames(counts_all)] %in% c(paste0("chr",1:22))]
  
  print("HI I AM IN aut FILTER BABE")

}else if(region_filter=="both"){
  # counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22),"chrX"))

  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22),"chrX"))
  # counts_subset_all<-counts_all[,genes_to_chrtype_dic[colnames(counts_all)] %in% c(paste0("chr",1:22),"chrX")]
  print("HI I AM IN both FILTER BABE")
  # print(head(genes_to_chrtype_dic[as.character(counts_all$Id)]))
  #  counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
=======

genes_to_chrtype<-fread(genes_to_chrtype_file)
genes_to_chrtype_dic<-genes_to_chrtype$V1
names(genes_to_chrtype_dic)<-genes_to_chrtype$V4
print(paste0("filtering to ",region_filter))
##filter by region x, aut, or both
if(region_filter=="x"){
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)]=="chrX")
  
  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id]=="chrX")
}else if(region_filter=="aut"){
  # counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22)))
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22)))
  
}else if(region_filter=="both"){
  counts_subset_all<-counts_all %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22),"chrX"))
 #  counts_subset_f<-counts_f %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
 # counts_subset_m<-counts_m %>% dplyr::filter(genes_to_chrtype_dic[Id] %in% c(paste0("chr",1:22),"chrX"))
}else{
  stop("ERROR: INCORRECT REGION FILTER. CHOOSE aut,both, or x ONLY IN LOWERCASE")
}

# print(paste0("using the folliwng chrs: ", 
#              unique(genes_to_chrtype_dic[rownames(counts_subset_all)])))

#only keep genes in both males and females
# overlapping_genes<-intersect(counts_subset_f$Id,counts_subset_m$Id)
# counts_f_touse<-counts_subset_f %>% dplyr::filter(Id %in% overlapping_genes)
# counts_m_touse<-counts_subset_m %>% dplyr::filter(Id %in% overlapping_genes)


### split into training and test
<<<<<<< HEAD
set.seed(90368)

print("iterating through models")
print((counts_subset_all)[1:5,1:6])
data_all_prefilt<-counts_subset_all %>% select(-Id)
data_all_lab_prefilt<-sex_dic[as.character(colnames(data_all_prefilt))]
f_lab<-names(data_all_lab_prefilt[data_all_lab_prefilt==2])
m_lab<-names(data_all_lab_prefilt[data_all_lab_prefilt==1])
counts_f_touse=data_all_prefilt[,f_lab]
counts_m_touse=data_all_prefilt[,m_lab]
print(paste0("NA? ",as.character(any(is.na(data_all_prefilt)))))

#sample half to train, half to test for male and feamle
# sampling_f<-as.data.frame(counts_f_touse[,!(names(counts_f_touse) %in% "Id")])[,sample(1:(ncol(counts_f_touse)-1),ncol(counts_f_touse)-1)]
# sampling_f_all<-as.data.frame(counts_f_touse[,!(names(counts_f_touse) %in% "Id")])
reorder_f<-sample(1:(ncol(counts_f_touse)),ncol(counts_f_touse))
sampling_f=as.data.frame(counts_f_touse)[,as.numeric(reorder_f)]
# rownames(sampling_f)<-counts_f_touse$Id
train_f<-sampling_f[,1:ceiling(ncol(sampling_f)*.5)]
test_f <- sampling_f[,(ceiling(ncol(sampling_f)*.5)+1):ncol(sampling_f)]

# sampling_m<-as.data.frame(counts_m_touse[,!(names(counts_m_touse) %in% "Id")])[,sample(1:(ncol(counts_m_touse)-1),ncol(counts_m_touse)-1)]
reorder_m<-sample(1:(ncol(counts_m_touse)),ncol(counts_m_touse))
sampling_m=as.data.frame(counts_m_touse)[,as.numeric(reorder_m)]
# rownames(sampling_m)<-counts_m_touse$Id
train_m<-sampling_m[,1:ceiling(ncol(sampling_m)*0.5)]
test_m <- sampling_m[,(ceiling(ncol(sampling_m)*0.5)+1):ncol(sampling_m)]

train_all<-cbind(train_f,train_m)
test_all<-cbind(test_f,test_m)
train_all_lab<-sex_dic[as.character(colnames(train_all))] 
test_all_lab<-sex_dic[as.character(colnames(test_all))] 
=======
# set.seed(90368)
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

>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
# sampling_all<-as.data.frame(counts_subset_all[,!(names(counts_subset_all) %in% "Id")])[,sample(1:(ncol(counts_subset_all)-1),ncol(counts_subset_all)-1)]
# rownames(sampling_all)<-counts_subset_all$Id
# print("separating test and train")
# in_train<-as.logical(group_dic[colnames(counts_subset_all)]=="train")
# in_train[is.na(in_train)]<-FALSE
# in_test<-as.logical(group_dic[colnames(counts_subset_all)]=="test")
# in_test[is.na(in_test)]<-FALSE

# train_all<-counts_subset_all[,group_dic[colnames(counts_subset_all)]=="train"]
# train_all<-counts_subset_all[,in_train]
# train_sex_lab<-sex_dic[colnames(train_all)]
# test_all <- counts_subset_all[,in_test]
# test_sex_lab<-sex_dic[colnames(test_all)]


# 
# print(paste0("using the folliwng chrs in train: ", 
#              unique(genes_to_chrtype_dic[rownames(train_all)])))
# print(paste0("using the folliwng chrs in test: ", 
#              unique(genes_to_chrtype_dic[rownames(test_all)])))

<<<<<<< HEAD

list.alphas <- seq(0,1,0.1)
preds=data.frame(matrix(ncol=0,nrow= ncol(test_all)))
=======
print("iterating through models")
data_all<-counts_subset_all %>% select(-Id)
data_all_lab<-sex_dic[as.character(colnames(data_all))]
print(paste0("NA? ",as.character(any(is.na(data_all)))))

list.alphas <- seq(0,1,0.1)
preds=data.frame(matrix(ncol=0,nrow= ncol(data_all)))
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
accuracy=data.frame(matrix(ncol=12,nrow=0))
colnames(accuracy)<-c("alpha","train_accuracy","test_accuracy")
for(my.alpha in list.alphas){
  print(paste0("running model with: ",my.alpha))
<<<<<<< HEAD
  # print("1")
  # print(head(data_all))
  # print(head(data_all_lab))
  # print(length(data_all_lab))
  # print(dim(data_all))
  # print(my.alpha)
  # print(colnames(data_all) %in% names(data_all_lab))
  # print(names(data_all_lab) %in%  colnames(data_all) )
  cvfit = cv.glmnet(t(train_all), train_all_lab, 
=======
  print("1")
  print(head(data_all))
  print(head(data_all_lab))
  print(length(data_all_lab))
  print(dim(data_all))
  print(my.alpha)
  print(colnames(data_all) %in% names(data_all_lab))
  print(names(data_all_lab) %in%  colnames(data_all) )
  cvfit = cv.glmnet(t(data_all), data_all_lab, 
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
                    family="binomial", 
                    alpha=my.alpha, nfolds = 6,type.measure="class",
                    standardize=F)
  print("2")
<<<<<<< HEAD
  preds_all <- predict(cvfit, newx=t(test_all), s="lambda.1se", type="response") #predict sex training
  print(head(preds_all))
  preds_class_all <- sapply(predict(cvfit, newx=t(test_all), s="lambda.1se", type="class"), as.numeric)
  print(head(preds_all))
  print("NOW ACC")
  all_acc <- sum(preds_class_all==test_all_lab)/length(test_all_lab) #accuracy
  print(head(all_acc))
  tmp_coeffs<-coef(cvfit,s="lambda.1se")
  print("COEFS")
  print(head(tmp_coeffs))
  coefs=data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  if(nrow(coefs) !=0){
    coefs$tissue<-tissue
    print(head(coefs))
    coefs$alpha<-my.alpha
    print(head(coefs))
    coefs$chr<-genes_to_chrtype_dic[as.character(coefs$name)] #NAMES to factors
    print(head(coefs))
    confuseMatrix<-confusionMatrix(as.factor(preds_class_all),as.factor(test_all_lab))
    print("HERE AND")
    acc<-data.frame("alpha"=my.alpha,t(as.data.frame(confuseMatrix$byClass)))
    
    # print(paste0("now in coefs: ",unique(coefs$chr)))
    this_coefs_file<-paste0(outdir,"/",tissue,"-coefs-",region_filter,"-alpha",my.alpha,"-",sextype,".txt")
    print(head(coefs))
    write.table(coefs,file=this_coefs_file,sep="\t",row.names=T,col.names=T, quote=FALSE)
    accuracy<-rbind(accuracy,acc)
    
  }else{    
    this_coefs_file<-paste0(outdir,"/",tissue,"-coefs-",region_filter,"-alpha",my.alpha,"-",sextype,".txt")
    write.table(coefs,file=this_coefs_file,sep="\t",row.names=T,col.names=T, quote=FALSE)
  }  
=======
  preds_all <- predict(cvfit, newx=t(data_all), s="lambda.1se", type="response") #predict sex training
  preds_class_all <- sapply(predict(cvfit, newx=t(data_all), s="lambda.1se", type="class"), as.numeric)
  all_acc <- sum(preds_class_all==data_all_lab)/length(data_all_lab) #accuracy
  
  tmp_coeffs<-coef(cvfit,s="lambda.1se")
  coefs=data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  coefs$tissue<-tissue
  coefs$alpha<-my.alpha
  coefs$chr<-genes_to_chrtype_dic[as.character(coefs$name)] #NAMES to factors
  
  confuseMatrix<-confusionMatrix(as.factor(preds_class_all),as.factor(data_all_lab))
  
  # print(paste0("now in coefs: ",unique(coefs$chr)))
  this_coefs_file<-paste0(outdir,"/",tissue,"-coefs-",region_filter,"-alpha",my.alpha,"-",sextype,".txt")
  write.table(coefs,file=this_coefs_file,sep="\t",row.names=T,col.names=T, quote=FALSE)
  acc<-data.frame("alpha"=my.alpha,t(as.data.frame(confuseMatrix$byClass)))
  accuracy<-rbind(accuracy,acc)
>>>>>>> 202c5a6ea887360d6e510761cb21611ce0d6089e
  this_preds<-data.frame(preds_all)
  colnames(this_preds)<-my.alpha
  preds<-cbind(preds,this_preds)
}
# in_train<-as.logical(group_dic[colnames(counts_subset_all)]=="train")
# in_train[is.na(in_train)]<-FALSE
# in_test<-as.logical(group_dic[colnames(counts_subset_all)]=="test")
# in_test[is.na(in_test)]<-FALSE
# 
# # train_all<-counts_subset_all[,group_dic[colnames(counts_subset_all)]=="train"]
# train_all<-counts_subset_all[,in_train]
# train_sex_lab<-sex_dic[colnames(train_all)]
# test_all <- counts_subset_all[,in_test]
# test_sex_lab<-sex_dic[colnames(test_all)]
# hyperparam_res2 <- lapply(list.alphas, function(my.alpha){
# for(my.alpha in list.alphas){
#   print(paste0("running model with: ",my.alpha))
#   
#   cvfit = cv.glmnet(t(train_all), train_sex_lab, 
#                     family="binomial", 
#                     alpha=my.alpha, nfolds = 6,type.measure="class",
#                     standardize=F)
#   
#   preds_train <- predict(cvfit, newx=t(train_all), s="lambda.1se", type="response") #predict sex training
#   preds_class_train <- sapply(predict(cvfit, newx=t(train_all), s="lambda.1se", type="class"), as.numeric)
#   train_acc <- sum(preds_class_train==train_sex_lab)/length(train_sex_lab) #accuracy
#   
#   preds_test <- predict(cvfit, newx=t(test_all), s="lambda.1se", type="response") #predict sex test
#   preds_class_test <- sapply(predict(cvfit, newx=t(test_all), s="lambda.1se", type="class"), as.numeric)
#   test_acc <- sum(preds_class_test==test_sex_lab)/length(test_sex_lab) 
#   # print(paste0("alpha=",my.alpha,", training accuracy=",train_acc,", test accuracy=",test_acc))
#   # print(paste0("following chrs in coef : ", 
#   #              unique(genes_to_chrtype_dic[rownames(coef(cvfit,s="lambda.1se"))])))
#   tmp_coeffs<-coef(cvfit,s="lambda.1se")
#   coefs=data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
#   coefs$tissue<-tissue
#   coefs$alpha<-my.alpha
#   coefs$chr<-genes_to_chrtype_dic[as.character(coefs$name)] #NAMES to factors
#   
#   confuseMatrix<-confusionMatrix(as.factor(preds_class_test),as.factor(test_sex_lab))
#   
#   # print(paste0("now in coefs: ",unique(coefs$chr)))
#   this_coefs_file<-paste0(outdir,"/",tissue,"-coefs-",region_filter,"-alpha",my.alpha,".txt")
#   write.table(coefs,file=this_coefs_file,sep="\t",row.names=T,col.names=T, quote=FALSE)
#   acc<-data.frame("alpha"=my.alpha,"train_accuracy"=train_acc,"test_accuracy"=test_acc,t(as.data.frame(confuseMatrix$byClass)))
#   accuracy<-rbind(accuracy,acc)
#   this_preds<-data.frame(preds_test)
#   colnames(this_preds)<-my.alpha
#   preds<-cbind(preds,this_preds)
# }

#   return(list(preds,accuracy))
# })



write.table(preds,file=outfile_preds,sep="\t",row.names=T,col.names=T, quote=FALSE)
write.table(accuracy,file=outfile_accuracy, sep="\t",row.names=F,col.names=T, quote=FALSE)



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
# my.alpha=.5
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
# 
# ggplot(to_plot,aes(x=prob,fill=sex,alpha=0.9))+geom_density() +xlim(c(0,1))+ggtitle(paste0(region_filter,": alpha=",my.alpha))


