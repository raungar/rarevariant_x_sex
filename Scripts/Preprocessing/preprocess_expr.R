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
                 make_option(c("-i", "--indir"), type="character", default=NULL, help="indir"),
                 make_option(c("-o", "--outdir"), type="character", default=NULL, help="outdir"),
                 make_option(c("-w", "--tpm_reads_dir"), type="character", default=NULL, help="where to read in tpm dir from"),
                 make_option(c("-d", "--randomize_sex"), type="character", default='', help="if you do not want to randomize sex leave as an empty string."),
                 make_option(c("-p", "--GTEX_PCv8"), type="character", default=NULL, help="GTEX_PCv8"),
                 make_option(c("-s", "--SUBJECTSv8"), type="character", default=NULL, help="GTEX_SUBJECTSv8"),
                 make_option(c("-j", "--SAMPLESv8"), type="character", default=NULL, help="GTEX_SAMPLESv8"),
                 make_option(c("-m", "--map_file"), type="character", default="gtex_2017-06-05_v8_samples_tissues.txt", help="output file samples"),
                 make_option(c("-t", "--do_tissues"), type="logical", default=NULL, help="T (first time run) or F (second time run)"),
                 make_option(c("-r", "--GTEX_RNAv8"), type="character", default=NULL, help="GTEX_RNAv8"),
                 make_option(c("-e", "--euro_file"), type="character", default=NULL, help="unambiguous european list"),
                 make_option(c("-g", "--TPM_DIF_FILE"), type="character", default="F", help="tpm differences otuptu file")
        ) 
        opt_parser = OptionParser(option_list=option_list)
        return(opt_parser)
}


##get most frequent individuals
#return type is a table
get_ind<-function(tissues, dir,euro){
  ind_list<-c()
  for(tissue in tissues){
    this_header = fread(paste0(dir,"/", tissue, '.tpm.txt'),nrows=1,header=F,sep="\t")
    ind_list<-c(ind_list,as.character(unlist(this_header)))
  }
  ind_table<-table(ind_list)
  ind_table_nogene<-rev(sort(ind_table[-1]))
  ind_table_nogene_euro<-ind_table_nogene[names(ind_table) %in% euro]
  return(ind_table_nogene_euro)

}

## For each tissue, read in the TPM and read files from the given directory.
## Subset and reorder columns of TPM file to match the given covariates.
## Filter for genes with > 20% individuals with TPM > 0.1 and read count > 6.
## RAU edit: in M or F (>20% of M OR >20% of F with TPM >0.1 and read count >6)
## Log2(tpm + 2) transform the data, then z-transform.
## Finally, output the transformed matrix to a new file.
ztrans.tissue = function(tissue, indir,outdir, covs, read.filt = 6, tpm.filt = 0.1,tpm_dif_file=F,inds,randomize_sex) {
    print("in ztrans")

    tpm = fread(paste0(indir,"/", tissue, '.tpm.txt'))
    reads = fread(paste0(indir,"/", tissue, '.reads.txt'))
    ## sanity checks
    stopifnot(sum(colnames(tpm) != colnames(reads)) == 0)
    print("0")
    if(sum(tpm$Gene != reads$Gene)){
      	print("WARNING: HAVING TO SUBSET GENES")
      	print(paste0("there are ", length(tpm$Gene[!(tpm$Gene %in% reads$Gene)]), " genes in TPM not in reads and there are ", 
      		length(tpm$Gene[!(reads$Gene %in% tpm$Gene)]), " genes in reads but not in TPM. "))	
      	#tpm<-tpm$Gene[(tpm$Gene %in% reads$Gene),]
      	tpm_red<-tpm %>% dplyr::filter(Gene %in% reads$Gene) %>%  dplyr::arrange(match(Gene,reads$Gene))
      	print(head(tpm_red$Gene))
      	tpm <- as.data.table(tpm_red )
      	print(head(tpm$Gene))
    }
    stopifnot(sum(tpm$Gene != reads$Gene) == 0)
    genes = tpm$Gene
    print(1)

    ### MALES AND FEMALES SEPARATELY NOW: RAU
    covs.subset = covs$SUBJID[covs$SUBJID %in% colnames(tpm)]
    covs.subset_full=covs %>% dplyr::filter(SUBJID %in% colnames(tpm))
    sex_table<-table(covs.subset_full$SEX)
    print(2)
    #check to see if single sex tissue
    if(length(sex_table)==1){
      print(3)
        print(paste0(tissue,": no subsetting, this is a sex-specific tissue"))
        tpm.single = tpm[, c(covs.subset), with = F]
        reads.single = reads[, c(covs.subset), with = F]
        ind.filt.single = round(0.2*ncol(tpm.single)) #20% of people
        #how many people pass 20% of tpm filtering and min reads numbering
        indices.keep.single = (rowSums(tpm.single > tpm.filt & reads.single > read.filt) >= ind.filt.single )
        tpm.cut.single = tpm[indices.keep.single, -1]
        tpm.out.single = scale(t(log2(tpm.cut.single + 2))) #log and z transform
        colnames(tpm.out.single) = genes[indices.keep.single]
        
        if(names(sex_table)==1){
          this_sex="m"
        }else if(names(sex_table)==2){
          this_sex="f"
        }else{
          stop("ERROR, sex does not make sense (not 1/2)")
        }
        
        write.table(tpm.out.single, paste0(outdir,"/", tissue, ".log2.ztrans.",this_sex,".txt"), quote = F, sep = '\t', row.names = T, col.names = T)
        
        ###DO SPECIAL THINGY HERE
        return()	
    }else{
      
      print(4)
      #subset individuals to males and females and both
      #keep 
      min_sex<-as.numeric(names(sex_table)[sex_table %in% min(sex_table)])
      max_sex<-as.numeric(names(sex_table)[sex_table %in% max(sex_table)])
      if(min_sex==1){
        mix_sex_fm="male"
        covs.m.all<-covs$SUBJID[covs$SEX==min_sex]
        covs.m<-covs.mf.all[covs.m.all %in% colnames(tpm)]
	      half_m<-floor(length(covs.m)/2)
	      covs.m.half<-names(rev(sort(inds[covs.m]))[1:half_m])
        covs.m.all<-covs$SUBJID[covs$SEX==max_sex]
        covs.f.sub<-covs.m.all[covs.f.all %in% colnames(tpm)]
        #sample to individuals that are in the most tissues
        covs.f<-names(rev(sort(inds[covs.f.sub]))[1:length(covs.m)]) 
        covs.f.half<-names(rev(sort(inds[covs.f.sub]))[1:half_F]) 
      }else if(min_sex==2){
        min_sex_fm="female"
        covs.f.all<-covs$SUBJID[covs$SEX==min_sex]
        covs.f<-covs.f.all[covs.f.all %in% colnames(tpm)]
      	half_f<-floor(length(covs.f)/2)
      	covs.f.half<-names(rev(sort(inds[covs.f]))[1:half_f])
         covs.m.all<-covs$SUBJID[covs$SEX==max_sex]
         covs.m.sub<-covs.m.all[covs.m.all %in% colnames(tpm)]
              #sample to individuals that are in the most tissues
         covs.m<-names(rev(sort(inds[covs.m.sub]))[1:length(covs.f)]) 
        covs.m.half<-names(rev(sort(inds[covs.m.sub]))[1:half_f]) 
      }else {
        min_sex_fm=stop("ERROR")
      }
      print(5)
      covs.both<-c(covs.m,covs.f)
      covs.both.half<-c(covs.m.half,covs.f.half)
      
      ####SUBSET, LOG TRANSFORM, AND Z TRANSFORM
      #M
      tpm.m = tpm[, c(covs.m), with = F]
      reads.m = reads[, c(covs.m), with = F]
      ind.filt.m = round(0.2*ncol(tpm.m)) #20% of people
      #how many people pass 20% of tpm filtering and min reads numbering
      indices.keep.m = (rowSums(tpm.m > tpm.filt & reads.m > read.filt) >= ind.filt.m )
      tpm.cut.m= tpm.m[indices.keep.m, -1]
      tpm.out.m = scale(t(log2(tpm.cut.m + 2))) #log and z transform
      colnames(tpm.out.m) = genes[indices.keep.m]
      write.table(tpm.out.m, paste0(outdir,"/", tissue, '.log2.ztrans.m.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
      
      print(6)
      #F
      tpm.f = tpm[, c(covs.f), with = F]
      print(6.2)
      reads.f = reads[, c(covs.f), with = F]
      print(6.3)
      ind.filt.f = round(0.2*ncol(tpm.f)) #20% of people
      print(7)
      #how many people pass 20% of tpm filtering and min reads numbering
      indices.keep.f = (rowSums(tpm.f > tpm.filt & reads.f > read.filt) >= ind.filt.f )
      print(7.5)
      tpm.cut.f = tpm.f[indices.keep.f, -1]
      print(8)
      tpm.out.f = scale(t(log2(tpm.cut.f + 2))) #log and z transform
      colnames(tpm.out.f) = genes[indices.keep.f]
      write.table(tpm.out.f, paste0(outdir,"/", tissue, '.log2.ztrans.f.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
      
      #Both
      print(paste0("WRITING BOTH:",outdir ,"/",tissue, '.log2.ztrans.allboth.txt'))
      tpm.both = tpm[, c(covs.both), with = F]
      reads.covs = reads[, c(covs.both), with = F]
      ind.filt.both = round(0.2*ncol(tpm.both)) #20% of people
      #how many people pass 20% of tpm filtering and min reads numbering
      #this passes female or male
      indices.keep.both = (rowSums(tpm.m > tpm.filt & reads.m > read.filt) >= ind.filt.m ) | (rowSums(tpm.f > tpm.filt & reads.f > read.filt) >= ind.filt.f)
      tpm.cut.both = tpm.both[indices.keep.both, -1]
      tpm.out.both = scale(t(log2(tpm.cut.both + 2))) #log and z transform
      colnames(tpm.out.both) = genes[indices.keep.both]
      write.table(tpm.out.both, paste0(outdir,"/", tissue, '.log2.ztrans.allboth.txt'), quote = F, sep = '\t', row.names = T, col.names = T)


      #Both half
      print(paste0("WRITING BOTH HALF:",outdir, tissue, '.log2.ztrans.both.txt'))
      tpm.both.half = tpm[, c(covs.both.half), with = F]
      reads.covs = reads[, c(covs.both.half), with = F]
      ind.filt.both.half = round(0.2*ncol(tpm.both.half)) #20% of people

      tpm.m.half = tpm[, c(covs.m.half), with = F]
      reads.m.half = reads[, c(covs.m.half), with = F]
      ind.filt.m.half = round(0.2*ncol(tpm.m.half)) #20% of people
      tpm.f.half = tpm[, c(covs.f.half), with = F]
      reads.f.half = reads[, c(covs.f.half), with = F]
      ind.filt.f.half = round(0.2*ncol(tpm.f.half)) #20% of people

      indices.keep.both.half = (rowSums(tpm.m.half > tpm.filt & reads.m.half > read.filt) >= ind.filt.m.half ) | (rowSums(tpm.f.half > tpm.filt & reads.f.half > read.filt) >= ind.filt.f.half)
      tpm.cut.both.half = tpm.both.half[indices.keep.both.half, -1]
      tpm.out.both.half = scale(t(log2(tpm.cut.both.half + 2))) #log and z transform
      colnames(tpm.out.both.half) = genes[indices.keep.both.half]
      write.table(tpm.out.both.half, paste0(outdir,"/", tissue, '.log2.ztrans.both.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
      
      ##OPTIONAL DEPENDING ON INPUT
      #write to a file if indicated how the tpm subsetting changed number of genes
      if(tpm_dif_file!= "F"){
        m_specific<-as.numeric(table(!(genes[indices.keep.m] %in% genes[indices.keep.f]))["TRUE"])
        f_specific<-as.numeric(table(!(genes[indices.keep.f] %in% genes[indices.keep.m]))["TRUE"])
        
        this_line<-data.frame(tissue,m_specific,f_specific)
        write.table(this_line, file=tpm_dif_file,sep="\t",quote = F,row.names = F,col.names = F,append = T)
      }
      ###RAU Commented out: code to compare different filtering
  
      
      # indices.keep.m = (rowSums(tpm.m > tpm.filt & reads.m > read.filt) >= ind.filt.m )
      # indices.keep.f = (rowSums(tpm.f > tpm.filt & reads.f > read.filt) >= ind.filt.f)
      # zero.keep.m = (rowSums(tpm.m > 0) >= zero.filt) 
      # zero.keep.f = (rowSums(tpm.f > 0) >= zero.filt)
      # indices.zero.keep.m = indices.keep.m[which(indices.keep.m %in% zero.keep.m)]
      # indices.zero.keep.f = indices.keep.f[which(indices.keep.f %in% zero.keep.f)]
      # genes.kept.m=as.character(t(tpm[indices.zero.keep.m, 1]))
      # genes.kept.f=as.character(t(tpm[indices.zero.keep.f, 1]))
      # 
      # library(ggVennDiagram)
      # my_venn_list<-list("mand2"=genes.kept.orig,
      #               "m"=genes.kept.m,
      #               "f"=genes.kept.f, 
      #               "mor2"=genes[indices.zero.keep])
      # my_venn_diagram<-ggVennDiagram(my_venn_list)+ scale_fill_gradient(low="blue",high = "purple")
      
      return()
    }
}

##------------- MAIN
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)

#read in files

indir = paste0(as.character(args$indir)) 
pc.file = as.character(args$GTEX_PCv8)
subject.file = as.character(args$SUBJECTSv8) 
sample.file = as.character(args$SAMPLESv8) 
map_file = as.character(args$map_file) 
GTEX_RNAv8=as.character(args$GTEX_RNAv8)
do_tissues=as.logical(args$do_tissues)
print(paste0("DO TISSUES: ",do_tissues))
tpm_reads_dir=as.character(args$tpm_reads_dir)
tpm_dif_file=as.character(args$TPM_DIF_FILE)
euro_file=as.character(args$euro_file)
euro_df<-read_tsv(euro_file,col_names = F)
euro<-euro_df$X1 
randomize_sex=as.character(args$randomize_sex)
outdir = paste0(as.character(args$indir),randomize_sex,"/PEER_v8") 
# randomize_sex=''
# indir = "/oak/stanford/groups/smontgom/raungar/Sex/Output"
# outdir=paste0(indir,"/preprocessing_v8",randomize_sex)
# tpm_reads_dir = paste0(dir, '/preprocessing_v8/PEER_v8/')
# pc.file = "/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt"
# subject.file = "/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
# sample.file = "/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt, /oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/rna_seq"
# map_file = "/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt"
# GTEX_RNAv8="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/rna_seq"
# do_tissues<-as.logical("T")
# tpm_dif_file="/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/tpm_dif_file.txt"
# euro_file="/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids_notambiguous.txt" #
# 

## Make output directory if it doesn't exist
#print(paste('mkdir -p', peer.dir))
system(paste('mkdir -p', outdir))

## Generate file mapping sample identifiers to tissues
## Restrict to samples that pass RNA-seq QC (marked as RNASEQ in column 28 of the sample file)
#map.file = paste0(dir, '/preprocessing_v8/',map_file_prefix) #gtex_2017-06-05_v8_samples_tissues.txt')
#command = paste("cat $GTEX_SAMPLESv8 | tail -n+2",
if(!do_tissues){
command = paste("cat ", sample.file ," | tail -n+2",
                "| cut -f1,14,28 | sed 's/ - /_/' | sed 's/ /_/g'",
                "| sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/'",
                "| awk '$3==\"RNASEQ\" {print $1\"\t\"$2}' >",map_file)
system(command)
}

  pcs = read.table(pc.file, header = T, stringsAsFactors = F)
  ## shorten names: e.g., only keep GTEX-1117F of GTEX-1117F-0003-SM-6WBT7
  pcs$SUBJID = apply(str_split_fixed(pcs$IID, '-', 5)[, c(1,2)], 1, paste, collapse = '-')
  pcs = pcs[, -c(1,2)]

if(do_tissues){
  ## Split the RPKM and read matrices into matrices for each tissue
  # DONE SEPARATELY
  

  ## Read in list of tissues and keep those with more than 50 samples
  tissue.list = read.table(map_file, header = F, stringsAsFactors = F)[,2]
  tissues = names(table(tissue.list)[table(tissue.list) > 50]) 
  tissues = tissues[tissues != 'Cells_Leukemia_cell_line_CML'] # exclude K-652 samples

  print(tpm_reads_dir)
  ind_table<-get_ind(tissues, tpm_reads_dir,euro)
  # Read in data for covariates and build the covariate matrix

  sex = read.csv(subject.file, header = T, stringsAsFactors = F, sep = '\t')[,c("SUBJID","SEX")]
  sex=sex%>%dplyr::filter(SUBJID%in%names(ind_table))
  print(paste0("random sex value: ",randomize_sex))
  #if we want to randomize the sex labels
  if(randomize_sex!=""){
    set.seed(NULL) #unset seed
    Sys.sleep(runif(1)/10) #wait a small random amoutn of time in case jobs are submitted at the same time
    op <- options(digits.secs = 6) #larger number of digits to make sure it's not randomized at the same second
    new_seed=as.numeric(Sys.time()) #ensure random seed by setting it to the system time
    print(paste0("NEW SEED: ",new_seed))
    set.seed(new_seed) #set new seed
  
    original_sex_labels=sex$SEX
    shuffled_sex_labels=sample(original_sex_labels)
    sex$SEX=shuffled_sex_labels
  }else{
    print("not randomizing sex")
  }
  
  
  
  covariates = merge(pcs, sex, by = 'SUBJID')



  print("Do tissues")
  if(tpm_dif_file != "F"){
    header<-data.frame("tissue","m_only","f_only")
    write.table(header, file=tpm_dif_file,sep="\t",quote = F,row.names = F,col.names = F)
  }
  #get the number of tissues each individual has (artificially randomly select for most number of tissues)
	sapply(tissues, ztrans.tissue, indir = tpm_reads_dir,outdir, covs = covariates,
	       tpm_dif_file=tpm_dif_file,inds=ind_table,randomize_sex)
} else{
	print("don't do tisues")
	#DONE SEPARATELY (RAU: ????)
	###system(make.split.command('tpm', map.file, peer.dir, GTEX_RNAv8))
	###system(make.split.command('reads', map.file, peer.dir, GTEX_RNAv8))

	# # Output the covariates matrix
	# write.table(covariates, paste0(outdir, '/covariates.txt'),
  #       	    col.names = T, row.names = F, quote = F, sep = '\t')
  	write.table(pcs, paste0(outdir, '/covariates.txt'),
        	    col.names = T, row.names = F, quote = F, sep = '\t')
}
