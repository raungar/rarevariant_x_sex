#!/usr/bin/env Rscript

rm(list = ls())

require(data.table)
require(plyr)
require(dplyr)

#-------------- MAIN
print("ENTERING R SCRIPT...")
print(date())
args = commandArgs(trailingOnly = T)
if (length(args) < 6) {
  cat("Usage: Rscript calculate_PEER_residuals.R RPKM COV PEER EQTLCALLS EQTLGENOS OUT\n", file = stderr())
  quit(status = 2)
}

## Define arguments
expr_file = args[1]
covs_file = args[2]
tissue_dir = args[3] #ex: /preprocessing/PEER_v8/Whole_Blood_Factors35
out_file = args[4]
metadata_file=args[5]
factors_type=args[6] #either factors or factors_sexregress
incl_sex=args[7]
sex_continuous=args[8]

print(out_file)

## Read in expression and covariate matrices
expr = read.table(expr_file, header = T, sep = '\t', row.names = 1)
covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)

## Reorder and subset rows in covariates file to match expression matrix rows
## Also only keep first 3 genotype PCs and sex (if applicable)
covs = covs[rownames(expr), ]
covs = covs[, c("PC1","PC2","PC3")]


## Read in PEER factors, fix subject names, and make column order match expression rows
peer_file<-paste0(tissue_dir,"/",factors_type,".tsv")
print("PEER FILE: ")
print(peer_file)
peer = t(read.table(peer_file, header = T, sep = '\t', row.names = 1))
rownames(peer) = gsub('\\.', '-', rownames(peer))
peer = peer[rownames(expr), ]
## Combine covariates and PEER factors
covs = cbind(covs, peer)

if(incl_sex == "T"){
	if(sex_continuous=="F"){
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
}

## Remove individuals with missing covariates
print("removing missing individuals")
inds_to_keep = rowSums(is.na(covs)) == 0
covs = covs[inds_to_keep, ]
expr = expr[inds_to_keep, ]

## For each gene in the expression file, perform a linear regression 
## Keep residuals
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)
print("performing regression")
for(i in 1:ncol(expr)){
	gene = names(expr)[i]
	data = as.data.frame(cbind(expr[, i], covs))
	colnames(data) = c('RPKM', colnames(covs))
	model = lm(RPKM ~ ., data = data)
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
