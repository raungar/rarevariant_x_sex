library(data.table)
library(dplyr)
library(tidyr)
library(optparse)


parse_my_args<-function(){
  option_list = list(
    make_option(c( "--combined_file"), type="character", default=NULL, help="z hg19", metavar="character"),
    make_option(c( "--outlierfile"), type="character", default=NULL, help="z hg19", metavar="character"),
    make_option(c("--maf_min"), type="character", default=NULL, help="outfile", metavar="character"),
    make_option(c("--cadd_min"), type="character", default=NULL, help="outfile", metavar="character"),
    make_option(c("--chr"), type="character", default=NULL, help="outfile", metavar="character"),
    make_option(c("--maf_max"), type="character", default=NULL, help="outfile", metavar="character"),
    make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt)
}

 #combined_file="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/Combined/x_all_rvs_inds_CADDtypesGQ5BlacklistRemovedALL_linc_prot_window5000.txt.gz"
 #outlierfile="/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_f_maxoutliers3.txt.gz"
combine_variant_outliers<-function(combined_file,outlierfile,chr,maf_max,maf_min){
  if(chr!='x'){this_chr=paste0('chr',chr)}else{this_chr=chr}
  print(this_chr)
  variants=fread(combined_file)
  colnames(variants)<-c('chr','pos1','pos2','maf1','maf2','maf_use','Ind','bps','Gene','genetype','sex','gene_loc',
                        'cadd_raw','cadd_phred','x','tss1','tss2','vartype','type2'  )
  variants_filt=variants%>%dplyr::filter(maf_use<maf_max & maf_use>=maf_min)
  exp=fread(outlierfile)%>%dplyr::filter(chr==this_chr)
  exp_variants=merge(exp,variants_filt,by=c('Gene','Ind'),all.x=T)
  print(head(exp_variants))
  exp_variants_long=exp_variants%>%mutate(vartype = strsplit(as.character(vartype), ",")) %>% 
    unnest(vartype)
  exp_variants_combined=exp_variants_long[,c('Gene','Ind','MedZ','Y','vartype','chr.x')]%>%rename(chr=chr.x)
  print(head(exp_variants_combined))
  return(exp_variants_combined)

}
get_prop_outliers<-function(exp_variants_combined){
  
  prop_outliers=exp_variants_combined%>%group_by(Y,vartype,chr)%>%mutate(count_status=n())%>%
    ungroup()%>%group_by(chr,Y)%>%mutate(num_genes_inds=n_distinct(Gene,Ind))%>%ungroup()%>%
    group_by(Y,vartype,chr,count_status,num_genes_inds)%>%
    summarise(prop_vars=count_status/num_genes_inds)%>%mutate(direction="all") %>%unique()
  
  prop_outliers_over_under=exp_variants_combined%>%mutate(direction=ifelse(sign(MedZ)==1,"over","under"))%>%
    group_by(Y,vartype,chr,direction)%>%mutate(count_status=n())%>%
    ungroup()%>%group_by(chr,Y,direction)%>%mutate(num_genes_inds=n_distinct(Gene,Ind))%>%ungroup()%>%
    group_by(Y,vartype,chr,count_status,num_genes_inds,direction)%>%
    summarise(prop_vars=count_status/num_genes_inds)%>%unique()

  prop_outliers_all=rbind(prop_outliers,prop_outliers_over_under)
  return(prop_outliers_all)
}

main<-function(){
  opt=parse_my_args()
  exp_variants_combined=combine_variant_outliers(opt$combined_file,opt$outlierfile,opt$chr,opt$maf_max,opt$maf_min)
  prop_outliers=get_prop_outliers(exp_variants_combined)
  print(head(prop_outliers))
  print(paste0("writing to ",opt$outfile))
  fwrite(prop_outliers,file=opt$outfile,col.names = F,row.names = F)
  print("wrote")
}
main()
#outliers 




