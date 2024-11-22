#!/usr/bin/env Rscript

## R script to split the entire GTEx matrices into files for each tissue
## then process RPKM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.
## output is peer_dir/tissue_[m/f/both].ztrans.txt

set.seed(1234) #IMPORTANT: RANDOM SUBSETTING OF INDIVIDUALS 
## Load required packages
require(data.table)
require(ggplot2)
require(stringr)
library(tidyverse)

library("optparse") #for passing args
##------------- FUNCTIONS
get_opt_parser<-function(){
        option_list = list(
                make_option(c("-d", "--RAREDIR"), type="character", default=NULL, help="RAREDIR directory"),
                 make_option(c("-s", "--SUBJECTSv8"), type="character", default=NULL, help="GTEX_SUBJECTSv8"),
                 make_option(c("-m", "--samples_tissues_file"), type="character", default="gtex_2017-06-05_v8_samples_tissues.txt", help="output file samples"),
                 make_option(c("-e", "--euro_file"), type="character", default=NULL, help="unambiguous european list"),
                make_option(c("-k", "--sex_new_old_key"), type="character", default="F", help="file with original sex and key to new sex"),
                make_option(c("-n", "--new_sex_file"), type="character", default="F", help="replace original samples files to get sex")
        ) 
        opt_parser = OptionParser(option_list=option_list)
        return(opt_parser)
}


##------------- MAIN
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)

dir = as.character(args$RAREDIR) #Sys.getenv('RAREDIR')
peer.dir = paste0(dir, '/preprocessing_v8/PEER_v8/')
subject.file = as.character(args$SUBJECTSv8) #Sys.getenv('GTEX_SUBJECTSv8')
samples_tissues_file = as.character(args$samples_tissues_file) #Sys.getenv('GTEX_SUBJECTSv8')
euro_file=as.character(args$euro_file)
sex_key_file=as.character(args$sex_new_old_key)
new_sex_file=as.character(args$new_sex_file)
# 
# dir = "/Volumes/groups/smontgom/raungar/Sex/Output"
# peer_dir_preprocessed = paste0(dir, '/preprocessing_v8/PEER_v8/')
# peer_dir_null = paste0(dir, '/nullshuffled_v8/PEER_v8/')
# pc.file = "/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt"
# subject.file = "/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
# sample.file = "/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt, /Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/rna_seq"
# samples_tissues_file = "/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt"
# GTEX_RNAv8="/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/rna_seq"
# do_tissues<-as.logical("T")
# tpm_dif_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/tpm_dif_file.txt"
# euro_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids_notambiguous.txt"

euro_df<-read_tsv(euro_file,col_names = F)
euro<-euro_df$X1

samples_tissues<-fread(samples_tissues_file,header=F)
# colnames(samples_tissues)<-c("SUBJID","Tissue")
## Read in list of tissues and keep those with more than 50 samples
tissue.list = read.table(samples_tissues_file, header = F, stringsAsFactors = F)[,2]
tissues = names(table(tissue.list)[table(tissue.list) > 50]) 
tissues = tissues[tissues != 'Cells_Leukemia_cell_line_CML'] # exclude K-652 samples

# Read in data for covariates and build the covariate matrix
subj = read.csv(subject.file, header = T, stringsAsFactors = F, sep = '\t')
print(head(subj))
sex = subj[,c(1,3)]
#IF I WANT TO RANDOMIZE GENERALLY sex$SEXsample(sex$SEX)

get_ind<-function(tissues, dir,euro){
  ind_list<-c()
  for(tissue in tissues){
    this_header = fread(paste0(dir, tissue, '.tpm.txt'),nrows=1,header=F,sep="\t")
    ind_list<-c(ind_list,as.character(unlist(this_header)))
  }
  ind_table<-table(ind_list)
  ind_table_nogene<-rev(sort(ind_table[-1]))
  ind_table_nogene_euro<-ind_table_nogene[names(ind_table) %in% euro]
  return(ind_table_nogene_euro)
  
}

  # print("Do tissues")
  # if(tpm_dif_file != "F"){
  #   header<-data.frame("tissue","m_only","f_only")
  #   write.table(header, file=tpm_dif_file,sep="\t",quote = F,row.names = F,col.names = F)
  # }
  #get the number of tissues each individual has (artificially randomly select for most number of tissues)
  ind_table<-get_ind(tissues, peer.dir,euro)
  sex_and_ntiss<-sex
  sex_and_ntiss$ntiss<-ind_table[sex_and_ntiss$SUBJID]
  sex_and_ntiss<-arrange((sex_and_ntiss),-ntiss)
  sex_and_ntiss_m<-sex_and_ntiss %>% filter(SEX==1)%>% filter(ntiss!="NA")
  sex_and_ntiss_f<-sex_and_ntiss %>% filter(SEX==2)%>% filter(ntiss!="NA")
  min_sex_amt<- min(nrow(sex_and_ntiss_f),nrow(sex_and_ntiss_m))
  sex_and_ntiss_f_subset<-sex_and_ntiss_f[1:min_sex_amt,]
  sex_and_ntiss_m_subset<-sex_and_ntiss_m[1:min_sex_amt,]
  sex_and_ntiss_sub<-rbind(sex_and_ntiss_m_subset,sex_and_ntiss_f)
  sex_and_ntiss_sub$SampledSex<-sample(sex_and_ntiss_sub$SEX) #SAMPLED
  new_sex_df=sex_and_ntiss_sub[ ,c("SUBJID","SampledSex")]
  colnames(new_sex_df)<-colnames(sex)
 
  write.table(sex_and_ntiss_sub,file=sex_key_file,sep="\t",row.names=F,col.names=FALSE, quote=FALSE)
  write.table(new_sex_df,file=new_sex_file,sep="\t",row.names=F,col.names=FALSE, quote=FALSE)
                      