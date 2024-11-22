#!/usr/bin/env Rscript

rm(list = ls())

require(data.table)
require(plyr)
require(dplyr)
library(tidyverse)
library("readr")


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
# tissue_dir = args[3] #ex: /preprocessing/PEER_v8/Whole_Blood_Factors35
eqtl_call_file=args[3]
eqtl_geno_file=args[4]
out_file = args[5]
metadata_file=args[6]
peer_factors=args[7] #either factors or factors_sexregress
incl_sex=args[8]
sex_continuous=args[9]

print("OUTFILE HERE")
print(out_file)

if(grepl("pinal_cord",eqtl_geno_file,fixed=T)){
  print("Then new eqtl_gene_file")
  eqtl_geno_file="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Brain_Spinal_cord_cervical_c-1.v8.egenes.txt.gz"

}

if(grepl("pinal_cord",eqtl_call_file,fixed=T)){
  print("Then new eqtl_gene_file")
  eqtl_call_file="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Brain_Spinal_cord_cervical_c-1.v8.egenes.txt.gz"
}

## Read in expression and covariate matrices
# expr = read.table(expr_file, header = T, sep = '\t', row.names = 1)
# covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)
expr = fread(expr_file)%>%as.data.frame()
covs = fread(covs_file)%>%as.data.frame()
rownames(expr)=expr[,1]
expr=expr[,-1]
print(head(expr[1:5,1:5]))
#print("expr^^, covs below")
rownames(covs)<-sapply(strsplit(covs[,1],"-"),function(x){paste0(x[1],'-',x[2])})
#print(head(covs[1:5,1:5]))
## Reorder and subset rows in covariates file to match expression matrix rows
## Also only keep first 3 genotype PCs and sex (if applicable)
covs = covs[rownames(expr), ]
covs = covs[, c("PC1","PC2","PC3")]


## Read in PEER factors, fix subject names, and make column order match expression rows
# peer_file<-paste0(tissue_dir,"/",factors_type,".tsv")
# print("PEER FILE: ")
# print(peer_file)
# peer = t(read.table(peer_factors, header = T, sep = '\t', row.names = 1))
peer = t(fread(peer_factors))
colnames(peer)<-paste0("Factor",1:ncol(peer))
# peer = t(read.table(peer_file, header = T, sep = '\t', row.names = 1))
print(rownames(expr))
rownames(peer) = gsub('\\.', '-', rownames(peer))
print(rownames(peer))
print(head(peer))

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


###Read in eQTL data for this tissue
## Restrict to the individuals with expression data for this tissue
print(paste0("reading eqtl_calls: ",eqtl_call_file))
eqtl_calls = read.table(eqtl_call_file, sep = '\t', header = T) %>% select(Gene = gene_id, Chrom = chr, Pos = variant_pos, Qval = qval)
print(paste0("reading eqtl geno file: ",eqtl_geno_file))
#eqtl_genos = as.data.frame(fread(eqtl_geno_file,header=T))
#gc() #memory issues trying this
#eqtl_genos = read.csv(eqtl_geno_file,header=T)
eqtl_genos=read_delim(eqtl_geno_file,'\t',col_names=TRUE)
print(head(eqtl_genos))
#eqtl_genos = readRDS(eqtl_geno_file)
#eqtl_genos = read.csv(eqtl_geno_file,header=T)
# gc() #memory issues trying this
# eqtl_genos = fread(eqtl_geno_file,header=T)
#eqtl_genos = (read.table(eqtl_geno_file))
print(0)
eqtl_genos = eqtl_genos[,c("Chrom", "Pos", rownames(expr))] %>%
			merge(., eqtl_calls)
print(1)

## For each gene in the expression file, perform a linear regression 
## Keep residuals
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
print(2)
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)
print("performing regression")
for(i in 1:ncol(expr)){
	gene = names(expr)[i]
	data = as.data.frame(cbind(expr[, i], covs))
	colnames(data) = c('RPKM', colnames(covs))
	if(gene %in% eqtl_genos$Gene){
		genos_temp = eqtl_genos %>% filter(Gene == gene)
		genos = genos_temp[,rownames(expr)] %>% t()		
		mean.geno = mean(genos, na.rm = T)
		genos = ifelse(is.na(genos), mean.geno, genos)
		data = cbind(data, genos)
		colnames(data)[ncol(data)] = 'EQTL'
	}
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
