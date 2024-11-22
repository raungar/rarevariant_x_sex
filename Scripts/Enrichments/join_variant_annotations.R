#!/usr/bin/env

## Merge outliers with variant annotations 

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)
library(stringr)
library(readr)

parser = ArgumentParser()
parser$add_argument('--method', help = 'Method') #medz or splicing or ase
parser$add_argument('--outfile', help = 'outfile')
parser$add_argument('--chr', help = 'chrom')
parser$add_argument('--variant_file', help='variant file') #/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv
parser$add_argument('--infile', help = 'input file for medz, splicing, or ase')
args = parser$parse_args()

# 
# infile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_f_maxoutliers3.txt.gz"
# this_variant_file="/Volumes/groups/smontgom/raungar/Sex/Output/features_v8/Collapsed/collapsed_chr16_CADDtypesALL_cadd15_window10000.txt.gz"
# method="medz"
#  this_chr="chrX"

method = args$method
this_variant_file<- as.character(args$variant_file)
outfile = as.character(args$outfile)
infile = as.character(args$infile)
this_chr = as.character(args$chr)
print(paste0("THIS CHR: ",this_chr))
if(!startsWith(this_chr,"chr") & this_chr !="all" & this_chr !="aut"){
  this_chr=paste0("chr",this_chr)
}
print(this_chr)
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
  #infile="/Volumes/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliersMuscle-Skeletal_noglobal_medz_zthresh2.5_nphen1_x_f_maxoutliers3.txt.gz"
  medz_outliers = fread(infile)
   print("OUTLIERS")
  print(head(medz_outliers))
  #colnames(medz_outliers)[1:2] = c('ind','ensg')
  colnames(medz_outliers)[1:2] = c('ind','ensg')
  #medz_outliers$chr<-str_replace(medz_outliers$chr,"chr","")
  #if(nchar(medz_outliers$chr[1])<=2){
  if(TRUE %in% head(nchar(medz_outliers$chr)<=2)){
    medz_outliers$chr<-paste0("chr",toupper(medz_outliers$chr))
  }
  if(this_chr=="aut"){
    medz_outliers<-medz_outliers%>%dplyr::filter(chr!="chrX")
    print(head(medz_outliers))
  }else if(this_chr=="all"){
    medz_outliers<-medz_outliers
  }else{
    print("HEREe in ELSE")
    print(head(medz_outliers));print(this_chr)
   medz_outliers<-medz_outliers%>%dplyr::filter(chr==this_chr)
  }
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
print(0)
#variant_file = read.table('/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv',header=T,fill=TRUE)
#this_variant_file="/Volumes/groups/smontgom/raungar/Sex/Output/features_v8/Collapsed/collapsed_x_CADDtypesGQ5BlacklistRemovedALL_cadd15_window10000.txt.gz"
# variant_file = read_tsv(this_variant_file,col_names=FALSE) #,header=F,fill=TRUE)
variant_file = fread(this_variant_file,header=F,fill=TRUE)
print("ORIGINAL...")
print(head(variant_file))
if(ncol(variant_file)==16){
	colnames(variant_file)<-c("chr","start","end","ensg","vartype","ind","sex","var_location","gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","geno","numrv")
}else if(ncol(variant_file)==19){
  colnames(variant_file)<-c("chr","start","end","gtex_maf","gnomad_maf","use_maf","ind","vartype","ensg","genetype","sex","var_location","cadd_raw","cadd_phred","geno","TSS_distance","TSE_distance","vep_type","vep_consq")
}else{
  colnames(variant_file)<-c("chr","start","end","ensg","vartype","ind","sex","var_location","gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","geno","TSS_distance","TSE_distance","numrv")
  
	# colnames(variant_file)<-c("chr","start","end","gtex_maf","gnomad_maf","use_maf","ind","vartype","genetype","sex","cadd_raw","cadd_phred","geno","num_rv")
	#variant_file=select(variant_file,chr,start,end,ensg,vartype,ind,sex,gtex_maf,gnomad_maf,use_maf,genetype,cadd_raw,cadd_phred,geno,TSS_distance,TSE_distance,num_rv)
}
print(head(variant_file))

#variant_file = read.table(this_variant_file,header=T,fill=NA)
variant_file$ensg = as.character(variant_file$ensg)

print(head(variant_file))
if(this_chr=="chrX" ){
  medz_outliers$chr<-paste0("chr",toupper(medz_outliers$chr))
  variant_file=variant_file %>% dplyr::filter(chr == "chrX")
  
}else if(this_chr=="aut"){
  variant_file=variant_file %>% dplyr::filter(chr != "chrX")  %>% dplyr::filter(chr !="NA")
  head(variant_file)
}else if(this_chr=="all"){
  variant_file=variant_file %>% dplyr::filter(chr !="NA")
}else{
 print("IN HERE")
variant_file=variant_file %>% dplyr::filter(chr == this_chr)
}
print("VARIANT FILE")
head(variant_file)
print("MEDZ OUTLIERS")
head(medz_outliers)

if (method == 'medz') {
  print(0)
  # medz_outliers_rename<-medz_outliers %>% mutate("ensg"=sapply(strsplit(ensg,"\\."),"[[",1))
  print(1)
    # variant_file_rename<-variant_file %>% mutate("ensg"=sapply(strsplit(ensg,"\\."),"[[",1))
  # medz_variants = merge(medz_outliers_rename,variant_file_rename,by=c('ind','ensg'),all.x=T)
  medz_variants = merge(medz_outliers,variant_file,by=c('ind','ensg'),all.x=T)
  print(2)
  
  medz_variants=medz_variants%>%mutate("ensg"=sapply(strsplit(ensg,"\\."),"[[",1))
  print(3)
  print(head(medz_variants))
  #write.table(medz_variants,file=paste0(RAREDIR,'/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.lncRNA.variantAnnotations.withPromoter.txt'),sep='\t',quote=F,row.names=F)
  write.table(medz_variants,file=outfile,sep='\t',quote=F,row.names=F)
} else if (method == 'splicing') {
  splice_variants = merge(splicing_melted,variant_file,by=c('ind','ensg'),all.x=T)
  #write.table(splice_variants,file=paste0(RAREDIR,'/data_v8/splicing/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.splice.variantAnnotations.withPromoter.txt'),sep='\t',quote=F,row.names=F)
  write.table(splice_variants,file=outfile,sep='\t',quote=F,row.names=F)

} else if (method == 'ase') {
  print(ase_outliers[1,])
  print(variant_file[1,])
  ase_variants = merge(ase_outliers,variant_file,by=c('ind','ensg'),all.x=T)
  #write.table(ase_variants,file=paste0(RAREDIR,'/data_v8/ase/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.v8.vg.txt'),sep='\t',quote=F,row.names=F)
  write.table(ase_variants,file=outfile,sep='\t',quote=F,row.names=F)
} else {
  cor_variants = merge(cor_outliers,variant_file,by.x=c('indiv_id','gene_id'),by.y=c('ind','ensg'),all.x=T)
  #write.table(cor_variants,file=paste0(RAREDIR,'/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.knn.variantAnnotations.txt'),sep='\t',quote=F,row.names=F)
  write.table(cor_variants,file=outfile,sep='\t',quote=F,row.names=F)

}

