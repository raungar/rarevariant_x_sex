library(data.table)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(purrr)
library(Matrix.utils) #aggregate.Matrix
library(magrittr) #set_colnames
library(glmnet)
library(SingleCellExperiment)
library("ggridges")
library("caret")
library("ggplot2")
sent_data<-readRDS("/oak/stanford/groups/smontgom/shared/TWC_ADRC_blood_scRNA/adrc_scrnaseq.rds.gz")
md=sent_data@meta.data
#number of UMI reads detected per cell (nCount_RNA)
summary(md$nCount_RNA)
count_data=sent_data@assays$RNA@counts
print(paste0("COUNTS (UMI): I will remove more than 2sd away from the median where the sd is ",sd(md$nCount_RNA), " and median is: ",median(md$nCount_RNA),
             ". So, anything greater than ", sd(md$nCount_RNA)*2+median(md$nCount_RNA), " or less than ",median(md$nCount_RNA)-sd(md$nCount_RNA)*2 ))
#number of expressed (detected) genes per same cell (nFeature_RNA).
summary(md$nFeature_RNA)
print(paste0("FEATURES: I will remove more than 2sd away from the median where the sd is ",sd(md$nFeature_RNA), " and median is: ",median(md$nFeature_RNA),
             ". So, anything greater than ", sd(md$nFeature_RNA)*2+median(md$nFeature_RNA), " or less than ",median(md$nFeature_RNA)-sd(md$nFeature_RNA)*2 ))
#percent mitochondiral
summary(md$percent.mt)
#percentage
colnames(md)[ncol(md)]
set.seed(12345)
sce <- SingleCellExperiment(assays = list(counts = count_data ), 
                            colData = md)
groups <- colData(sce)[, c("celltype", "Sample")]%>% as.data.frame%>%mutate(Sample=str_replace_all(Sample,"_","-"))
rownames(groups)<-rownames(colData(sce))
##pseudobulking goes here?
#from: https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
pseudobulked <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
# num_clusters <- length(cluster_ids <- purrr::set_names(levels(sce$celltype_cluster)))
# num_samples <- length(sample_ids <- purrr::set_names(levels(sce$SampleI)))

# pseudobulked_formatted<- split.data.frame(pseudobulked, rep(num_clusters, num_samples)) %>% 
#   lapply(function(u) set_colnames(t(u), unname(sample_ids)))
# splitf <- sapply(stringr::str_split(rownames(pseudobulked),
#                                     pattern = "_",
#                                     n = 2), `[`, 1)
# pseudobulked_renamed<-split.data.frame(pseudobulked, factor(splitf)) %>%
#   lapply(function(u)
#     set_colnames(t(u),
#                stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
# 
# 
# get_sample_ids <- function(x){pseudobulked[[x]] %>%colnames()}
# get_cluster_ids <- function(x){
#   rep(names(pseudobulked)[x],  each = length(samples_list[[x]]))
# }
# samples_list <- map(1:length(cluster_ids), get_sample_ids)
# 
# de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
#   unlist()
# ggplot(md %>% mutate(type="nFeature_RNA"),aes(x=type,y=nFeature_RNA))+geom_violin()
#cells with a percentage of mitochondrial genes below 0.05% were included. Cells with the highest (top 0.2%) or lowest (bottom 0.2%) numbers of detected genes were considered as outliers and excluded from the downstream analyses
#removed low-quality cells according to the standard, which is cells with fewer than 200 unique molecular identifiers (UMIs) or mitochondrial gene expression exceeding 60%. 
#n. Any PBMC with more than 7% of mitochondrial UMI counts was considered to be a low-quality cell [12]. PBMC GEMs with greater than 2500 genes expressed or CD8 GEMs with more than 2000 detected genes were checked in order to determine the rate of doublets.
#Any gene detected in less than three cells or a cell with less than 200 genes detected was excluded for downstream data analysis. https://link.springer.com/article/10.1186/s12967-018-1578-4#Sec2
#Cells that expressed more than 2,500 genes, more than 10,000 unique molecular identifiers (UMIs) and more than 10% mitochondrial genes were excluded. #https://www.nature.com/articles/s41586-019-1895-7#Sec2
#Genes were excluded if they were expressed in fewer than 10 cells, and cells were excluded if they expressed fewer than 200 genes. 
celltypes<-unique(md$celltype)
visits<-paste0("Y",unique(md$Visit))
for(this_celltype in celltypes){
  for(this_visit in visits){
    print(paste0("looking at ",this_celltype,"in year ",this_visit))
    # red_data <- sent_data #subset(sent_data, subset = nFeature_RNA > 200 & nFeature_RNA < (2+median(md$nFeature_RNA)) & percent.mt <= 7)
    # red_data <- pseudobulked #subset(sent_data, subset = nFeature_RNA > 200 & nFeature_RNA < (2+median(md$nFeature_RNA)) & percent.mt <= 7)
    # red_data_visit1<-subset(red_data, subset = Visit == 1)
    # red_data<-RunPCA(red_data)
    # red_data_visit1<-RunPCA(red_data_visit1)
    # DimPlot(red_data_visit1,reduction="pca",group.by = "Sex")
    rownames_celltype <- sapply(stringr::str_split(rownames(pseudobulked),pattern = "_", n = 2), `[`, 1)
    rownames_visit <- sapply(strsplit(rownames(pseudobulked),"-"), "[[", 2)
   
    red_data<-pseudobulked[ which((rownames_celltype == this_celltype & rownames_visit==this_visit)==T),]

    if(is.vector(red_data)){print(dim(red_data)); print( "THIS IS TOO SMALL");next } ## for when just one
    if(dim(red_data)[1]<30){print(dim(red_data)); print( "THIS IS TOO SMALL");next }
    print(dim(red_data))
    print(red_data[1:5,1:5])
    rnaseq=red_data
    # rnaseq_t<-t(red_data)
     rnaseq<-(red_data)
    test_rows<-sample(1:nrow(rnaseq),nrow(rnaseq)/2)
    train_rows<-(1:nrow(rnaseq))[!1:nrow(rnaseq) %in% test_rows]
    rnaseq_test<-rnaseq[test_rows,]
    rnaseq_train<-rnaseq[train_rows,]
    test_sampleids<-paste0(sapply(strsplit(rownames(rnaseq[test_rows,]),"\\_|\\-"),"[[",2),"_",this_visit)
    train_sampleids<-paste0(sapply(strsplit(rownames(rnaseq[train_rows,]),"\\_|\\-"),"[[",2),"_",this_visit)
    md_test_sex<-md%>% dplyr::filter(celltype==this_celltype  & Sample%in%test_sampleids)%>%select(Sample,Sex) %>% unique() %>%pull(Sex)
    md_train_sex<-md%>% dplyr::filter(celltype==this_celltype  & Sample%in%train_sampleids)%>%select(Sample,Sex) %>% unique() %>%pull(Sex)
    
    my.alpha=0.5
    cvfit = cv.glmnet(rnaseq_train, md_train_sex, 
                      family="binomial", 
                      alpha=my.alpha, nfolds = 6,type.measure="class",
                      standardize=F)
    preds_all <- predict(cvfit, newx=rnaseq_train, s="lambda.1se", type="response") #predict sex training
    write.table(data.frame(preds_all),paste0("/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/SingleCell/preds_test_",this_celltype,"_Y",this_visit,".txt"))
    # preds_class_all <- sapply(predict(cvfit,  newx=rnaseq_train, s="lambda.1se", type="class"), as.numeric)
    #all_acc <- sum(preds_class_all==data_all_lab)/length(data_all_lab) #accuracy
    
    tmp_coeffs<-coef(cvfit,s="lambda.1se")
    coefs=data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    coefs$alpha<-my.alpha
    write.table(coefs,file=paste0("/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/SingleCell/coefs_",this_celltype,"_Y",this_visit,".txt"),sep="\t",row.names=T,col.names=T, quote=FALSE)
  }
}
if(1 ==0){
preds_test<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/SingleCell/preds_test_Dendritic cell.txt")
colnames(preds_test)<-c("cell","sexcont")
preds_w_md<-cbind(preds_test,md[preds_test$cell,],"md_cell"=rownames(md[preds_test$cell,]))
preds_w_md_byind<-preds_w_md[,c("cell","sexcont","Patient","Diagnosis","Age","Sex","Diagnosis_path")] %>%
  group_by(Patient,Age,Sex,Diagnosis_path)%>%
  summarise(Median=median(sexcont), Mean=mean(sexcont), Min=min(sexcont), Max=max(sexcont),  StandDev=sd(sexcont),IQR=IQR(sexcont))%>%
  mutate(contsex_binary=ifelse(Median>0.5,"male","female"))
preds_w_md_byind_subtypes<-preds_w_md[,c("cell","sexcont","Patient","Diagnosis","Age","Sex","Diagnosis_path","celltype_cluster")] %>%
  group_by(Patient,Age,Sex,Diagnosis_path,celltype_cluster)%>%
  summarise(Median=median(sexcont), Mean=mean(sexcont), Min=min(sexcont), Max=max(sexcont),  StandDev=sd(sexcont),IQR=IQR(sexcont))%>%
  mutate(contsex_binary=ifelse(Median>0.5,"male","female"))

tmp=preds_w_md_byind_subtypes %>% dplyr::filter(Diagnosis_path=="Mild Cognitive Impairment->Probable Alzheimer's Disease")
ggplot(preds_w_md,aes(x=Diagnosis_path,y=sexcont,fill=Sex))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values=c("#dbab3b","#5d8596"))+ #"#926fa8",
  #geom_point()
  ggtitle("Sex scores across test set cells, CD8")+
  geom_violin()
to_plot=preds_w_md# %>% dplyr::filter(celltype_cluster=="Dend")
ggplot(to_plot,aes(x=sexcont,y=Diagnosis_path,fill=Sex,alpha=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values=c("#dbab3b","#5d8596"))+ #"#926fa8",
  xlim(c(0,1))+
  #geom_point()
  ggtitle("Sex scores across test set cells, CD8 T cell 4")+
  geom_density_ridges()

table(preds_w_md_byind$Sex==preds_w_md_byind$contsex_binary)
table(preds_w_md_byind_subtypes$Sex==preds_w_md_byind_subtypes$contsex_binary)

mismatch=preds_w_md_byind_subtypes%>% dplyr::filter(contsex_binary!=Sex)
table(mismatch$Sex)
table(mismatch$Diagnosis_path)
f_preds=preds_w_md_byind %>% dplyr::filter(Sex=="female")
f_preds_w_md_byind_cont_forlm=data.frame("Sex"=as.numeric(f_preds$Median),
                                         as.data.frame(model.matrix(~Diagnosis_path,f_preds)))
m_preds=preds_w_md_byind %>% dplyr::filter(Sex=="male")
m_preds_w_md_byind_cont_forlm=data.frame("Sex"=as.numeric(m_preds$Median),
                                         as.data.frame(model.matrix(~Diagnosis_path,m_preds)))
summary(lm(Sex~.,m_preds_w_md_byind_cont_forlm))

preds_w_md_byind_bin_forlm=data.frame("Sex"=as.numeric(preds_w_md_byind$Sex)-1,
                                           as.data.frame(model.matrix(~Diagnosis_path,preds_w_md_byind)))
preds_w_md_byind_cont_forlm=data.frame("Sex"=as.numeric(preds_w_md_byind$Median),
                                           as.data.frame(model.matrix(~Diagnosis_path,preds_w_md_byind)))
lm_binsex=lm((Sex)~.,preds_w_md_byind_bin_forlm)
lm_contsex=lm((Sex)~.,preds_w_md_byind_cont_forlm)


}
