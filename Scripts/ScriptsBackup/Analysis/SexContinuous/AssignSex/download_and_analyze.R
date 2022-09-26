library("recount3")
library("dplyr")
library("ggplot2")
library("recount")

#whole blood rna-seq of 379 samples across 3 timepoints from 157 inds with lupus
lupus=recount3::create_rse_manual(
  project = "SRP150872",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
#SRP156583
sle_bcells<-recount3::create_rse_manual(
  project = "SRP156583",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
assign_sex<-function(counts,gene_dic,region_filter){
  

  ##filter by region x, aut, or both
  if(region_filter=="x"){
    counts_subset_all<-counts %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)]=="chrX")
  }else if(region_filter=="aut"){
    counts_subset_all<-counts %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22)))
  }else if(region_filter=="both"){
    counts_subset_all<-counts %>% dplyr::filter(genes_to_chrtype_dic[as.character(Id)] %in% c(paste0("chr",1:22),"chrX"))
  }else{
    stop("ERROR: INCORRECT REGION FILTER. CHOOSE aut,both, or x ONLY IN LOWERCASE")
  }
  
  data_all<-counts_subset_all %>% select(-Id)
  data_all_lab<-sex_dic[as.character(colnames(data_all))]

  # list.alphas <- seq(0,1,0.1)
  list.alphas <- c(0,0.5,1)
  preds=data.frame(matrix(ncol=0,nrow= ncol(data_all)))
  accuracy=data.frame(matrix(ncol=12,nrow=0))
  colnames(accuracy)<-c("alpha","train_accuracy","test_accuracy")
  for(my.alpha in list.alphas){
    print(paste0("running model with: ",my.alpha))
    cvfit = cv.glmnet(t(data_all), data_all_lab, 
                      family="binomial", 
                      alpha=my.alpha, nfolds = 6,type.measure="class",
                      standardize=F)
    preds_all <- predict(cvfit, newx=t(data_all), s="lambda.1se", type="response") #predict sex training
    preds_class_all <- sapply(predict(cvfit, newx=t(data_all), s="lambda.1se", type="class"), as.numeric)
    all_acc <- sum(preds_class_all==data_all_lab)/length(data_all_lab) #accuracy
    
    tmp_coeffs<-coef(cvfit,s="lambda.1se")
    coefs=data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    coefs$tissue<-tissue
    coefs$alpha<-my.alpha
    coefs$chr<-genes_to_chrtype_dic[as.character(coefs$name)] #NAMES to factors
    
    confuseMatrix<-confusionMatrix(as.factor(preds_class_all),as.factor(data_all_lab))
    acc<-data.frame("alpha"=my.alpha,t(as.data.frame(confuseMatrix$byClass)))
    accuracy<-rbind(accuracy,acc)
    this_preds<-data.frame(preds_all)
    colnames(this_preds)<-my.alpha
    preds<-cbind(preds,this_preds)
  }
}

sleb_samplemd<-data.frame(sample=rownames(colData(sle_bcells)),md=colData(sle_bcells)$sra.sample_attributes)
sleb_samplemd_arrange<-colData(sle_bcells)$sra.sample_attributes

sleb_samplemd <- cbind(rownames(colData(sle_bcells)),data.frame(t(sapply(strsplit(colData(sle_bcells)$sra.sample_attributes,"[;|\\|]"), `[`))))
colnames(sleb_samplemd)<-c("sample_id","na1","empty1","cell_subtype","na2","empty2","cell_type",
                           "na3","empty3","sorting_markers","na4","empty4","source_name","na5","empty5","status")
sleb_sample_md<-sleb_samplemd  %>% select(sample_id,cell_subtype,cell_type,sorting_markers,source_name,status)
md_dic_status<-sleb_sample_md$status
names(md_dic_status)<-sleb_sample_md$sample_id
md_dic_cellsubtype<-sleb_sample_md$cell_subtype
names(md_dic_cellsubtype)<-sleb_sample_md$sample_id
ggplot(sleb_sample_md,aes(x=cell_subtype))+geom_bar() + 
facet_wrap(~status)+theme_bw()+  theme(axis.text.x = element_text(angle = 45,  hjust=1))

sle_bcells_rawcounts=assays(sle_bcells)$raw_counts
assays(sle_bcells)$counts<-transform_counts(sle_bcells)
assays(sle_bcells)$TPM<-transform_counts(sle_bcells)
assays(sle_bcells)$TPM <-getTPM(sle_bcells,length_var = "score")


genes_to_chrtype_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtf_padded10kb.bed"
genes_to_chrtype<-fread(genes_to_chrtype_file)
genes_to_chrtype_dic<-genes_to_chrtype$V1
names(genes_to_chrtype_dic)<-genes_to_chrtype$V4
tpms_filt<-assays(sle_bcells)$TPM[]
tpm_filt=0.5;ind_filt=0.3
indices.keep.single = (rowSums(assays(sle_bcells)$TPM > tpm_filt) >= ind_filt )
print(table(indices.keep.single))
red_tpm=assays(sle_bcells)$TPM[indices.keep.single,]

pca_tpms<-prcomp(t(red_tpm))
to_plot<-cbind.data.frame(md_dic_status[rownames(pca_tpms$x)],md_dic_cellsubtype[rownames(pca_tpms$x)],pca_tpms$x)
colnames(to_plot)[1:2]<-c("status","cell_subtype")
pc_eigenvalues <- tibble(PC = factor(1:length(pca_tpms$sdev^2)), 
                         variance = pca_tpms$sdev^2) %>%   mutate(pct = variance/sum(variance)*100) %>%  mutate(pct_cum = cumsum(pct))
ggplot(pc_eigenvalues[1:50,],aes(x = PC)) +geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")
ggplot(to_plot,aes(x=PC1,y=PC3,color=cell_subtype,shape=status,size=2,alpha=0.8))+geom_point()+theme_bw()
assign_sex_res_x<-assign_sex(sle_bcells_rawcounts,genes_to_chrtype_dic,"x")


