#!/usr/bin/env Rscript

rm(list = ls())

require(data.table)
require(plyr)
require(dplyr)
library(stringr)

#-------------- MAIN
print("ENTERING R SCRIPT...")

args = commandArgs(trailingOnly = T)
if (length(args) < 6) {
  cat("Usage: Rscript calculate_PEER_residuals.R RPKM COV PEER EQTLCALLS EQTLGENOS OUT\n", file = stderr())
  quit(status = 2)
}

## Define arguments
expr_file = args[1]
covs_file = args[2]
tissue_dir = args[3] #ex: /preprocessing/PEER_v8/Whole_Blood_Factors35
eqtl_call_file = args[4]
eqtl_geno_file = args[5]
out_file = args[6]
metadata_file=args[7]
factors_type=args[8] #either factors or factors_sexregress
incl_sex=args[9]

print(out_file)

#expr_file = '/srv/scratch/restricted/GOATs/preprocessing/PEER_v7/Whole_Blood.rpkm.log2.ztrans.txt' 
#covs_file = '/srv/scratch/restricted/GOATs/preprocessing/PEER_v7/covariates.txt'
#peer_file = '/srv/scratch/restricted/GOATs/preprocessing/PEER_v7/Whole_Blood_Factors35/factors.tsv'
#eqtl_call_file = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_FOR_QC_ONLY/Whole_Blood_Analysis.v6p.FOR_QC_ONLY.egenes.txt.gz'
#eqtl_geno_file = '/srv/scratch/restricted/GOATs/preprocessing/gtex_2016-01-15_v7_genotypes_v6p_cis_eQTLs_012_processed.txt'
#out_file = '/srv/scratch/restricted/GOATs/preprocessing/PEER_v7/Whole_Blood.peer.v6pciseQTLs.ztrans.txt'

## Read in expression and covariate matrices
expr = read.table(expr_file, header = T, sep = '\t', row.names = 1)
covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)
#covs=read.csv(covs_file, header = T, stringsAsFactors = F, sep = '\t')
#print(head(covs,1))
##if(as.numeric(str_count(rownames(covs),"-"))>1){
#	before_rn<-rownames(covs)
#	before_rn_split<-strsplit(before_rn,"-")
#	print("NEW RN")
#	new_rn<-paste0(sapply(before_rn_split,"[[",1),"-",sapply(before_rn_split,"[[",2))
#	print(new_rn)
#	rownames(covs)<-new_rn
#}
#print(head(rownames(covs)))

#global_outliers = fread('/users/nferraro/data/goats_data/v8_data/gtexV8_global_outliers_medz3_iqr.txt', header=F)
#kinds = which(!(rownames(expr) %in% global_outliers$V1))
#expr = expr[kinds,]

## Reorder and subset rows in covariates file to match expression matrix rows
## Also only keep first 3 genotype PCs and sex (if applicable)
print("HEAD EXPR RN")
print(head(rownames(expr)))
covs = covs[rownames(expr), ]
#print("FIRST HEAD: COVS")
#print(head(covs))
#remove.cols = paste0('C', 4:20)
#covs = covs[, !(colnames(covs) %in% remove.cols)]
covs = covs[, c("PC1","PC2","PC3")]


## Read in PEER factors, fix subject names, and make column order match expression rows
peer_file<-paste0(tissue_dir,"/",factors_type,".tsv")
print("PEER FILE: ")
print(peer_file)
peer = t(read.table(peer_file, header = T, sep = '\t', row.names = 1))
rownames(peer) = gsub('\\.', '-', rownames(peer))
#print("PEER THEN EXPR")
#print(head(peer))
#print(head(rownames(expr)))
peer = peer[rownames(expr), ]

## Combine covariates and PEER factors
#num_keep = round(0.5*ncol(peer))
#peer = peer[,1:num_keep]
covs = cbind(covs, peer)

print("HEAD COVS 1")
print(head(covs,1))

#print("NROW COVS")
#print(nrow(covs))


if(incl_sex == "T"){
	print("INCLUDING SEX")
	md_sex<-fread(metadata_file,sep="\t",header=T)
	print(head(md_sex[,c("SEX","SUBJID")]))
	ind_dic<-as.numeric(md_sex$SEX)-1
	names(ind_dic)<-md_sex$SUBJID
	print(head(ind_dic))

	before_rn_split<-strsplit(rownames(covs),"-")
	new_rn<-paste0(sapply(before_rn_split,"[[",1),"-",sapply(before_rn_split,"[[",2))
	print("NEW RN")
	print(head(new_rn))

	
	sex_vals=ind_dic[new_rn]
	print(sex_vals)
	print("sex_vals")	
	covs=cbind(covs,"SEX"=sex_vals)
	print(head(covs,1))
}

#print("COVS post sex")
#print(head(covs))



## Remove individuals with missing covariates
inds_to_keep = rowSums(is.na(covs)) == 0
covs = covs[inds_to_keep, ]
expr = expr[inds_to_keep, ]

## Read in eQTL data for this tissue
## Restrict to the individuals with expression data for this tissue
eqtl_calls = read.table(eqtl_call_file, sep = '\t', header = T) %>% select(Gene = gene_id, Chrom = chr, Pos = variant_pos, Qval = qval)
eqtl_genos = as.data.frame(fread(eqtl_geno_file,header=T))

#print("eqtl genos")
#head(eqtl_genos[1:5,1:5])
#head(rownames(expr))
#print("rownames expr in colnames genos")
#print(rownames(expr) %in% colnames(eqtl_genos))
#print("colnames genos in rownames expr")
#print(colnames(eqtl_genos) %in% rownames(expr))
#print(head(eqtl_calls))

#print("DFSKLFJLDSKFJ")
#eqtl_genos = eqtl_genos %>% dplyr::select_(c("Chrom", "Pos", rownames(expr))) %>%
#			merge(., eqtl_calls)
eqtl_genos = eqtl_genos[,c("Chrom", "Pos", rownames(expr))] %>%
			merge(., eqtl_calls)
 

## For each gene in the expression file, perform a linear regression 
## Keep residuals
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)

#print("EXPR")
#head(expr[,1])
#print("COVS")
#print(head(covs,2))

for(i in 1:ncol(expr)){
	gene = names(expr)[i]
	data = as.data.frame(cbind(expr[, i], covs))
	colnames(data) = c('RPKM', colnames(covs))
	#print("DATA")
	#print(head(data,1))
	if(gene %in% eqtl_genos$Gene){
		genos_temp = eqtl_genos %>% filter(Gene == gene) # %>% select_(rownames(expr)) %>% t()		
		genos = genos_temp[,rownames(expr)] %>% t()		
		#print("GENOS")
		#print(head(genos,1))
		#genos = eqtl_genos %>% filter(Gene == gene) %>% select_(rownames(expr)) %>% t()		
		mean.geno = mean(genos, na.rm = T)
		#print("MEAN GENO")
		#print(head(mean.geno))
		genos = ifelse(is.na(genos), mean.geno, genos)
		data = cbind(data, genos)
		colnames(data)[ncol(data)] = 'EQTL'
		#print("FINAL DATA")
		#print(head(data,1))
	}
	model = lm(RPKM ~ ., data = data)
	#print("MODEL SUMMARY")
	#print(summary(model))
	resids[, i] = model$residuals
}

print("center and scale")
# Center and scale, then transpose
resids = t(scale(resids))

print("write")
# Write out the residuals
write.table(matrix(c('Id', colnames(resids)), nrow = 1), out_file, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(resids, out_file, row.names = T, col.names = F, quote = F, sep = '\t', append = T)

print("DONE")
