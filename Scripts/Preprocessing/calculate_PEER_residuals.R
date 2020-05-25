#!/usr/bin/env Rscript

rm(list = ls())

require(data.table)
require(plyr)
require(dplyr)

#-------------- MAIN

args = commandArgs(trailingOnly = T)
if (length(args) != 6) {
  cat("Usage: Rscript calculate_PEER_residuals.R RPKM COV PEER EQTLCALLS EQTLGENOS OUT\n", file = stderr())
  quit(status = 2)
}

## Define arguments
expr_file = args[1]
covs_file = args[2]
peer_file = args[3]
eqtl_call_file = args[4]
eqtl_geno_file = args[5]
out_file = args[6]

## Read in expression and covariate matrices
expr = read.table(expr_file, header = T, sep = '\t', row.names = 1)
covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)
print(head(covs,1))

## Reorder and subset rows in covariates file to match expression matrix rows
## Also only keep first 3 genotype PCs and sex (if applicable)
covs = covs[rownames(expr), ]
covs = covs[, c("PC1","PC2","PC3")]




## Read in PEER factors, fix subject names, and make column order match expression rows
peer = t(read.table(peer_file, header = T, sep = '\t', row.names = 1))
rownames(peer) = gsub('\\.', '-', rownames(peer))
peer = peer[rownames(expr), ]

## Combine covariates and PEER factors
#num_keep = round(0.5*ncol(peer))
#peer = peer[,1:num_keep]
covs = cbind(covs, peer)



## Remove individuals with missing covariates
inds_to_keep = rowSums(is.na(covs)) == 0
covs = covs[inds_to_keep, ]
expr = expr[inds_to_keep, ]

## Read in eQTL data for this tissue
## Restrict to the individuals with expression data for this tissue
eqtl_calls = read.table(eqtl_call_file, sep = '\t', header = T) %>% select(Gene = gene_id, Chrom = chr, Pos = variant_pos, Qval = qval)
eqtl_genos = as.data.frame(fread(eqtl_geno_file,header=T))

eqtl_genos = eqtl_genos[,c("Chrom", "Pos", rownames(expr))] %>%
			merge(., eqtl_calls)
 

## For each gene in the expression file, perform a linear regression 
## Keep residuals
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)

for(i in 1:ncol(expr)){
	gene = names(expr)[i]
	data = as.data.frame(cbind(expr[, i], covs))
	colnames(data) = c('RPKM', colnames(covs))
	if(gene %in% eqtl_genos$Gene){
		genos_temp = eqtl_genos %>% filter(Gene == gene) # %>% select_(rownames(expr)) %>% t()		
		genos = genos_temp[,rownames(expr)] %>% t()		
		mean.geno = mean(genos, na.rm = T)
		genos = ifelse(is.na(genos), mean.geno, genos)
		data = cbind(data, genos)
		colnames(data)[ncol(data)] = 'EQTL'
	}
	model = lm(RPKM ~ ., data = data)
	print("MODEL SUMMARY")
	print(summary(model))
	resids[, i] = model$residuals
}

print("center and scale")
# Center and scale, then transpose
resids = t(scale(resids))

print("write")
# Write out the residuals
write.table(matrix(c('Id', colnames(resids)), nrow = 1), out_file, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(resids, out_file, row.names = T, col.names = F, quote = F, sep = '\t', append = T)
