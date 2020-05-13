#!/usr/bin/env Rscript

## R script to split the entire GTEx matrices into files for each tissue
## then process RPKM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.

## Load required packages
require(data.table)
require(ggplot2)
require(stringr)

library("optparse") #for passing args
##------------- FUNCTIONS
get_opt_parser<-function(){
        option_list = list(
                 make_option(c("-d", "--RAREDIR"), type="character", default=NULL, help="RAREDIR directory"),
                 make_option(c("-p", "--GTEX_PCv8"), type="character", default=NULL, help="GTEX_PCv8"),
                 make_option(c("-s", "--SUBJECTSv8"), type="character", default=NULL, help="GTEX_SUBJECTSv8"),
                 make_option(c("-j", "--SAMPLESv8"), type="character", default=NULL, help="GTEX_SAMPLESv8"),
                 make_option(c("-m", "--map_file_prefix"), type="character", default="gtex_2017-06-05_v8_samples_tissues.txt", help="output file samples"),
                 make_option(c("-t", "--do_tissues"), type="logical", default=NULL, help="T (first time run) or F (second time run)"),
                 make_option(c("-r", "--GTEX_RNAv8"), type="character", default=NULL, help="GTEX_RNAv8")
                 ) 
        opt_parser = OptionParser(option_list=option_list)
        return(opt_parser)
}


## Function to put together the command to split the given matrix by tissue
##    datatype should be either 'read' or 'rpkm'
##    mapfile should be the file with the sampleID-tissue correspondance
##    outdir should be the directory to write the new files to
make.split.command = function(datatype, mapfile, outdir, GTEX_RNAv8) {
   # comm = paste0('python preprocessing/correction/preprocess_expr_split_by_tissues.py ',
    comm = paste0('python preprocessing/split_expr_by_tissues.py ',
                  '--gtex ', GTEX_RNAv8, '/GTEx_Analysis_v8_RNA-seq_RNA-SeQCv1.1.8_gene_',
                  datatype, '.gct.gz ',
                  '--sample ', mapfile, ' --out ', outdir,
                  ' --end .', datatype, '.txt')
    return(comm)
}

## For each tissue, read in the TPM and read files from the given directory.
## Subset and reorder columns of TPM file to match the given covariates.
## Filter for genes with > 20% individuals with TPM > 0.1 and read count > 6.
## Log2(tpm + 2) transform the data, then z-transform.
## Finally, output the transformed matrix to a new file.
ztrans.tissue = function(tissue, dir, covs, read.filt = 6, tpm.filt = 0.1) {
    tpm = fread(paste0(dir, tissue, '.tpm.txt'))
    reads = fread(paste0(dir, tissue, '.reads.txt'))
    ## sanity checks
    stopifnot(sum(colnames(tpm) != colnames(reads)) == 0)
    
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
    print("did not stop")
    covariates.subset = covs$SUBJID[covs$SUBJID %in% colnames(tpm)]
    genes = tpm$Gene
    tpm = tpm[, c(covariates.subset), with = F]
    reads = reads[, c(covariates.subset), with = F]
    ind.filt = round(0.2*ncol(tpm))
    zero.filt = round(0.05*ncol(tpm))
    indices.keep = rowSums(tpm > tpm.filt & reads > read.filt) >= ind.filt
    zero.keep = rowSums(tpm > 0) >= zero.filt
    indices.keep = indices.keep[which(indices.keep %in% zero.keep)]
    tpm = tpm[indices.keep, ]
    tpm.out = scale(t(log2(tpm + 2))) 
    colnames(tpm.out) = genes[indices.keep]
    write.table(tpm.out, paste0(dir, tissue, '.log2.ztrans.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
    return()
}

##------------- MAIN
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)

dir = as.character(args$RAREDIR) #Sys.getenv('RAREDIR')
peer.dir = paste0(dir, '/preprocessing_v8/PEER_v8/')
pc.file = as.character(args$GTEX_PCv8) #Sys.getenv('GTEX_PCv8')
subject.file = as.character(args$SUBJECTSv8) #Sys.getenv('GTEX_SUBJECTSv8')
sample.file = as.character(args$SAMPLESv8) #Sys.getenv('GTEX_SUBJECTSv8')
map_file_prefix = as.character(args$map_file_prefix) #Sys.getenv('GTEX_SUBJECTSv8')
GTEX_RNAv8=as.character(args$GTEX_RNAv8)
do_tissues=as.logical(args$do_tissues)

print("args start ---")
print(peer.dir)
print(pc.file)
print(subject.file)
print(map_file_prefix)
print(GTEX_RNAv8)
print(do_tissues)
print("--- end args")
## Make output directory if it doesn't exist
#print(paste('mkdir -p', peer.dir))
system(paste('mkdir -p', peer.dir))

## Generate file mapping sample identifiers to tissues
## Restrict to samples that pass RNA-seq QC (marked as RNASEQ in column 28 of the sample file)
map.file = paste0(dir, '/preprocessing_v8/',map_file_prefix) #gtex_2017-06-05_v8_samples_tissues.txt')
#command = paste("cat $GTEX_SAMPLESv8 | tail -n+2",
command = paste("cat ", sample.file ," | tail -n+2",
                "| cut -f1,14,28 | sed 's/ - /_/' | sed 's/ /_/g'",
                "| sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/'",
                "| awk '$3==\"RNASEQ\" {print $1\"\t\"$2}' >", map.file)
print(command)
system(command)

## Split the RPKM and read matrices into matrices for each tissue
# DONE SEPARATELY

## Read in list of tissues and keep those with more than 50 samples
tissue.list = read.table(map.file, header = F, stringsAsFactors = F)[,2]
tissues = names(table(tissue.list)[table(tissue.list) > 50]) 
tissues = tissues[tissues != 'Cells_Leukemia_cell_line_CML'] # exclude K-652 samples

# Read in data for covariates and build the covariate matrix
pcs = read.table(pc.file, header = T, stringsAsFactors = F)
## shorten names: e.g., only keep GTEX-1117F of GTEX-1117F-0003-SM-6WBT7
pcs$SUBJID = apply(str_split_fixed(pcs$IID, '-', 5)[, c(1,2)], 1, paste, collapse = '-')
pcs = pcs[, -c(1,2)]
sex = read.csv(subject.file, header = T, stringsAsFactors = F, sep = '\t')[, c(1,3)]
covariates = merge(pcs, sex, by = 'SUBJID')


if(do_tissues){
	sapply(tissues, ztrans.tissue, dir = peer.dir, covs = covariates)
} else{
	#DONE SEPARATELY (RAU: ????)
	#system(make.split.command('tpm', map.file, peer.dir, GTEX_RNAv8))
	#system(make.split.command('reads', map.file, peer.dir, GTEX_RNAv8))

	# Output the covariates matrix
	write.table(covariates, paste0(peer.dir, '/covariates.txt'),
        	    col.names = T, row.names = F, quote = F, sep = '\t')
}
