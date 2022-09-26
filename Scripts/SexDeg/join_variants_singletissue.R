#!/usr/bin/env

## Merge outliers with variant annotations 

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)


parser = ArgumentParser()
parser$add_argument('--method', help = 'Method') #ztrans or splicing or ase
parser$add_argument('--outfile', help = 'outfile')
parser$add_argument('--sexdeg_file', help = ' sex deg file')
parser$add_argument('--chr', help = 'chr')
parser$add_argument('--linc_prot_only_dir', help = ' linc_prot_only_dir')
parser$add_argument('--variant_file', help='variant file') #/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv
parser$add_argument('--infile', help = 'input file for ztrans, splicing, or ase')
args = parser$parse_args()
print(args)
method = args$method
outfile = as.character(args$outfile)
infile = as.character(args$infile)
sexdeg_file=as.character(args$sexdeg_file)
mychr=as.character(args$chr)
if(!(startsWith(mychr,"chr")) & mychr != "aut" & mychr !="all"){
  mychr=paste0("chr",mychr)
}
this_variant_file<- as.character(args$variant_file)
linc_prot_only_dir<-as.character(args$linc_prot_only_dir)
# 
# infile="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/Liver.f.peer.ztrans.txt"
# this_variant_file="/Volumes/groups/smontgom/raungar/Sex/Output/features_v8/Collapsed/collapsed_chr7_CADDtypesALL_cadd15.txt.gz"
# sexdeg_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/all_genes_sexDEGs_beta0.111_LIVER.txt"
 #linc_prot_only_dir<-"/Volumes/groups/smontgom/raungar/Sex/Output/features_v8/GTF_lnc_protcod_only"
linc_prot_genes<-c()
for (f in list.files(linc_prot_only_dir,full.names = T)){
  this_f=fread(f,select=8)
  linc_prot_genes<-c(linc_prot_genes,data.frame(this_f)[,1])
}


this_ztrans = fread(infile) %>%dplyr::filter(Id %in% linc_prot_genes)
this_ztrans_melt<-melt(this_ztrans)
colnames(this_ztrans_melt) = c("ensg","ind","z")
  
variant_file = fread(this_variant_file,header=F,fill=TRUE)
print(head(variant_file))
#chr7	100791968	100791968	ENSG00000146839.18	SNPs	GTEX-14XAO	female	0.074930	0.0765954	0.0765954	protein_coding	2.90618	23.2	9
colnames(variant_file)<-c("chr","pos","pos2","ensg","vartype","ind","sex","gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","numrv")
#colnames(variant_file)<-c("chr","pos","ensg","vartype","ind","sex","ref","alt","both","gtex_maf","gnomad_both","gnomad_m","gnomad_f","genetype","mafdiff","numrv")
variant_file$ensg = as.character(variant_file$ensg)
sexdegs<-fread(sexdeg_file,header=F,data.table = F)
colnames(sexdegs)<-c("chr","ensg","beta")
print("VARIANT FILE")
head(variant_file)
print("this ztrans")
head(this_ztrans_melt)
print("SEXDEGS")
print(head(sexdegs))




  ztrans_sexdegs = merge(this_ztrans_melt,sexdegs,by=c('ensg'),all.x=T) #%>%dplyr::filter(chr != "chrX")
  if(mychr=="all"){
    sexdegs_filtered<-ztrans_sexdegs
  }else if(mychr=="aut"){
    sexdegs_filtered<-ztrans_sexdegs %>% dplyr::filter(chr !="chrX")
    
  }else{
    sexdegs_filtered<-ztrans_sexdegs %>% dplyr::filter(chr ==mychr)
    print( head(sexdegs_filtered))
    
  }
  
  ztrans_variants_sexdeg=merge(sexdegs_filtered,variant_file,by=c('ind','ensg'),all.x=T)
  write.table(ztrans_variants_sexdeg,file=outfile,sep='\t',quote=F,row.names=F)

