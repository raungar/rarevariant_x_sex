#!/usr/bin/env

## Merge outliers with variant annotations 

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)


parser = ArgumentParser()
parser$add_argument('--method', help = 'Method') #medz or splicing or ase
parser$add_argument('--outfile', help = 'outfile')
parser$add_argument('--sexdeg_file', help = ' sex deg file')
parser$add_argument('--variant_file', help='variant file') #/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv
parser$add_argument('--infile', help = 'input file for medz, splicing, or ase')
args = parser$parse_args()

method = args$method
outfile = as.character(args$outfile)
infile = as.character(args$infile)
sexdeg_file=as.character(args$sexdeg_file)
this_variant_file<- as.character(args$variant_file)

get_bin <- function(pval) {
  if (pval < 1e-13) {
    bin = '0~1e-13'
  } else if (pval < 1e-11) {
    bin = '1e-13~1e-11'
  } else if (pval < 1e-09) {
    bin = '1e-11~1e-09'
  } else if (pval < 1e-07) {
    bin = '1e-09~1e-07'
  } else if (pval < 1e-05) {
    bin = '1e-07~1e-05'
  } else if (pval < 1e-04) {
    bin = '1e-05~1e-04'
  } else if (pval < 1e-03) {
    bin = '1e-04~1e-03'
  } else if (pval < 1e-02) {
    bin = '1e-03~1e-02'
  } else if (pval < 5e-02) {
    bin = '1e-02~5e-02'
  } else {
    bin = 'nonOutlier'
  }
  return(bin)
}

if (method == 'medz') {
  #medz_outliers = fread(paste0(RAREDIR, '/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'))
  medz_outliers = fread(infile)
  #colnames(medz_outliers)[1:2] = c('ind','ensg')
  colnames(medz_outliers)[1:2] = c('ind','ensg')
  #medz_outliers$ind = sapply(medz_outliers$ind, function(x) strsplit(x,'GTEX-')[[1]][2])
} else if (method == 'splicing') {
  #splicing_outliers = fread(paste0(RAREDIR, '/data_v8/splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
  splicing_outliers = fread(infile)
  splicing_melted = melt(splicing_outliers)
  splicing_melted = filter(splicing_melted,!is.nan(value))
  splicing_bins = sapply(splicing_melted$value, function(x) get_bin(x))
  splicing_melted$pval_bin = splicing_bins
  #colnames(splicing_melted)[1:2] = c('ensg','ind')
  colnames(splicing_melted)[1:2] = c('ensg','ind')
  #splicing_melted$indiv_id = as.character(splicing_melted$indiv_id)
  #splicing_melted$indiv_id = sapply(splicing_melted$indiv_id, function(x) strsplit(x,'GTEX-')[[1]][2])
  splicing_melted$ind = as.character(splicing_melted$ind)
  splicing_melted$ind = sapply(splicing_melted$ind, function(x) strsplit(x,'GTEX-')[[1]][2])
  outlier_genes = filter(splicing_melted,pval_bin!='nonOutlier')$ensg
  splicing_melted = filter(splicing_melted,ensg %in% outlier_genes)
} else if (method == 'ase') {
  #ase_data = fread(paste0(RAREDIR, '/data_v8/ase/combined.ad.scores.in.MEDIAN_4_10_update.tsv'))
  #ase_data = fread('/users/nferraro/data/goats_data/v8_data/ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.tsv', stringsAsFactors=F)
  ase_data=fread(infile,stringsAsFactors=F)
  ase_outliers = melt(ase_data)
  ase_outliers$pval_bin = sapply(ase_outliers$value, function(x) ifelse(is.na(x), NA, get_bin(x)))
  colnames(ase_outliers)[1:2] = c('ensg','ind')
  ase_outliers$ind = as.character(ase_outliers$ind)
  ase_outliers$ind = sapply(ase_outliers$ind, function(x) strsplit(x,'GTEX-')[[1]][2])
  gene_map = fread('/users/nferraro/data/goats_data/v8_data/gencode.v26.GRCh38.genes.bed',header=F)
  gene_map$ensg = sapply(gene_map$V4, function(x) strsplit(x, '[.]')[[1]][1])
  gene_map = gene_map %>% select(V4,ensg)
  colnames(gene_map) = c('Gene', 'ensg')
  ase_outliers = inner_join(ase_outliers, gene_map, by='ensg')
  ase_outliers = ase_outliers %>% select(Gene,ind,value,pval_bin)
  colnames(ase_outliers)[1] = 'ensg'
  
} else {
  #cor_data = fread(paste0(RAREDIR, '/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.newParams.knn.txt'))
  cor_data = fread(infile)
  cor_outliers = filter(cor_data,Y=='outlier',FDR < 0.01)
  cor_controls = filter(cor_data, Gene %in% cor_outliers$Gene, Y=='control')
  cor_bins = sapply(cor_outliers$FDR, function(x) get_bin(x))
  cor_outliers$pval_bin = cor_bins
  cor_controls$pval_bin = 'nonOutlier'
  cor_outliers = rbind(cor_outliers %>% select(Ind,Gene,FDR,pval_bin,Y),cor_controls %>% select(Ind,Gene,FDR,pval_bin,Y))
  #colnames(cor_outliers)[1:2] = c('ind','ensg')
  colnames(cor_outliers)[1:2] = c('ind','ensg')
  cor_outliers$ind = sapply(cor_outliers$ind, function(x) strsplit(x,'GTEX-')[[1]][2])
}

#variant_file = read.table('/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv',header=T,fill=TRUE)
variant_file = fread(paste0("zcat -f ",this_variant_file),header=F,fill=TRUE)
#variant_file = read.table(this_variant_file,header=T,fill=NA)
colnames(variant_file)<-c("chr","pos","ensg","vartype","ind","sex","ref","alt","both","gtex_maf","gnomad_both","gnomad_m","gnomad_f","genetype","mafdiff","numrv")
variant_file$ensg = as.character(variant_file$ensg)
sexdegs<-fread(sexdeg_file,header=F)
colnames(sexdegs)<-c("chr","ensg","beta")
print("VARIANT FILE")
head(variant_file)
print("MEDZ OUTLIERS")
head(medz_outliers)
print("SEXDEGS")
print(head(sexdegs))

if (method == 'medz') {
  medz_variants = merge(medz_outliers,variant_file,by=c('ind','ensg'),all.x=T)
  medz_variants_sexdeg=merge(medz_variants,sexdegs,by=c('ensg'),all.x=T)
  #write.table(medz_variants,file=paste0(RAREDIR,'/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.lncRNA.variantAnnotations.withPromoter.txt'),sep='\t',quote=F,row.names=F)
  write.table(medz_variants_sexdeg,file=outfile,sep='\t',quote=F,row.names=F)
} 
