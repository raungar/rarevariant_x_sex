#!/usr/bin/env Rscript

rm(list = ls())

require(data.table)
require(plyr)
library(stringr)
require(dplyr)
library(optparse)
library(reshape2)


#--- OPTION PARSER
option_list = list(
  make_option(c("-e", "--counts_f"), type="character", default=NULL, help="counts file", metavar="character"),
  make_option(c("-c", "--pc_file"), type="character", default=NULL, help="pcs", metavar="character"),
  make_option(c("-f", "--factors_file"), type="character", default=NULL, help="[dir]/factors.tsv or [dir]/sex_regress.factors.tsv etc ", metavar="character"),
  make_option(c("-s", "--incl_sex"), type="character", default=NULL, help="incl_sex or regress_sex or no_sex", metavar="character"),
  make_option(c("-n", "--is_null"), type="logical", default=NULL, help="T or F incl sex", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-m", "--metadata_file"), type="character", default=NULL, help="metadata file", metavar="character"),
  make_option(c("-p", "--peer_file"), type="character", default=NULL, help="metadata file", metavar="character"),
  make_option(c("-b", "--sex_continuous"), type="logical", default=NULL, help="Sex continuous? T or F", metavar="character")
  
) 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
counts_f<-as.character(opt$counts_f)
pc_file<-as.character(opt$pc_file)
incl_sex<-as.logical(opt$character)
is_null<-as.logical(opt$incl_sex) #null distr or no
tissue_dir<-as.character(opt$tissue_dir)
out_file<-as.character(opt$out_file)
metadata_file<-as.character(opt$metadata_file)
peer_file<-as.character(opt$peer_file)
sex_continuous<-as.logical(opt$sex_continuous)

# ## Define arguments
# expr_file = args[1]
# covs_file = args[2]
# tissue_dir = args[3] #ex: /preprocessing/PEER_v8/Whole_Blood_Factors35
# out_file = args[4]
# metadata_file=args[5]
# factors_type=args[6] #either factors or factors_sexregress
# incl_sex=args[7]
# sex_continuous=args[8]

print(out_file)
 #counts_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/Adipose_Subcutaneous.log2.ztrans.both_half.txt"
# pc_file<-"/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt"
# factors_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/Adipose_Subcutaneous_Factors30_m/factors.tsv"
#metadata_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/Regressions/md_Adipose_Subcutaneous-preprocessing_v8-0.5-aut.txt"
# metadata_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/Assignments/test_train.txt"
## Read in expression and covariate matrices
counts<-data.frame(fread(counts_file),row.names=1)
pcs = data.frame(fread(pc_file, header = T, sep = '\t'),row.names=1)
rownames(pcs)<-sapply(str_split((rownames(pcs)),"-"),function(x){paste0(x[[1]],"-",x[[2]])})
## Reorder and subset rows in covariates file to match countsession matrix rows
## Also only keep first 3 genotype PCs and sex (if applicable)
pcs_red = pcs[rownames(counts),c("PC1","PC2","PC3") ]



## Read in PEER factors, fix subject names, and make column order match expression rows
peer = t(data.frame(fread(factors_file, header = T, sep = '\t'), row.names = 1))
rownames(peer) = gsub('\\.', '-', rownames(peer))
peer = peer[rownames(counts), ]
## Combine covariates and PEER factors

if(sex_continuous=="F"){
  print("treating sex as a binary variable")
		md_sex_full<-data.frame(fread(metadata_file,sep="\t"),row.names=1)
		colnames(md_sex_full)<-c("SEX","ntiss","SexScrambled","test_or_train","swap_id")
		
    if(is_null=="F"){
  		  md_sex=md_sex_full%>%dplyr::select(SEX)
  		  md_sex$SEX<-md_sex$SEX-1
    }else if (is_null=="T"){
        md_sex=md_sex_full%>%dplyr::select(SexScrambled)
        md_sex$SexScrambled<-md_sex$SexScrambled-1
    }else{
        stop("ERROR: IS_NULL IS INCORRECTLY ASSIGNED")
    }

}else if (sex_continuous=="T"){
  print("treating sex as a continuous variable")
  
	  md_sex<-data.frame(fread(metadata_file,sep="\t"),row.names=1)
}else{
	    stop("ERROR: INCORRECT ASSIGNMETN FOR SEX_CONTINUOUS")
}

if(sex == "incl_sex"){
  
  peer=cbind(peer,"SEX"=md_sex[rownames(covs),])
  print(head(peer,1))
  print("INCLUDING SEX")
  
}else if(sex=="regress_sex"){
  ##THIS REGRESSES SEX FROM PER TO PROTECT IT FROM BEING REGRESSED OUT
  tmp_peer=cbind(peer,"SEX"=md_sex[rownames(covs),])
  lm_form<-as.formula(paste0("cbind(",paste0(colnames(peer),collapse=","),") ~ SEX"))
  lm_fit<-lm(lm_form,data=as.data.frame(tmp_peer))
  peer<-(residuals(lm_fit))
}else if(sex=="no_sex"){
  peer=peer
}else{
  stop("ERROR: INVALID SEX VARIABLE")
}
covs = cbind(pcs_red, peer)

## Remove individuals with missing covariates
print("removing missing individuals")
inds_to_keep = rowSums(is.na(covs)) == 0
covs = covs[inds_to_keep, ]
counts = counts[inds_to_keep, ]

## For each gene in the countsession file, perform a linear regression 
## Keep residuals
resids = matrix(, ncol = ncol(counts), nrow = nrow(counts))
rownames(resids) = rownames(counts)
colnames(resids) = colnames(counts)
print("performing regression")
#will get counts regresesd on covariates
for(i in 1:ncol(counts)){
	gene = names(counts)[i]
	data = as.data.frame(cbind(counts[, i], covs))
	colnames(data) = c('RPKM', colnames(covs))
	model = lm(RPKM ~ ., data = data)
	resids[, i] = model$residuals
}

print("center and scale")
# Center and scale, then transpose
resids = t(scale(resids))

print("write")
# Write out the residuals
#outfile=${prefix}.${sex}.peer.ztrans.txt
#resids are counts regressed on metadata
write.table(matrix(c('Id', colnames(resids)), nrow = 1), out_file, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(resids, out_file, row.names = T, col.names = F, quote = F, sep = '\t', append = T)

print("DONE")

  inds_all<-colnames(resids)[-1] #colnames are the individuals

  #set dependent variable
  # resids_forlm<-as.matrix(resids[,-1])
  # rownames(resids_forlm)<-resids[,1]
  
  #combine these variables into one dataframe for lm()
  lm_matrix<-cbind.data.frame(t(resids),md_sex[colnames(resids),])
  #print(head(lm_matrix))
  #formula where it is essentially factors ~ sex, specifically cbind(Factor1, Factor2, ..., FactorN) ~ Sex
  lm_form<-as.formula(paste0("cbind(",paste0(rownames(resids),collapse=","),") ~ ind_sex"))
  #calculate the lm fit and residuals
  #print(lm_form)
  lm_fit<-lm(lm_form,data=as.data.frame(lm_matrix))
  lm_resid<-residuals(lm_fit)
  lm_resid_t<-(t(lm_resid))
  lm_resid_t<-cbind("Id"=rownames(lm_resid_t),lm_resid_t)
  #write this to the new folder
  #write.table(lm_resid,file=paste0(resid_dir_path,"/residuals_sex.tsv"),sep="\t",quote=F)
  write.table(lm_resid_t,file=fileout,sep="\t",quote=F,row.names=F)



