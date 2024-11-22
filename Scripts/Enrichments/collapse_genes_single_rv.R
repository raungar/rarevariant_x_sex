library(dict)
#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--outfile"), type = 'character', default = NULL, help = "path of output file whcih is an RData file")
        )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
outfile <- as.character(opt$outfile)

# zscore<-2
#  infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen5_noglobal_medz_varAnnot_x_m.txt"
#  #out_rdata_relative<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_z3_x_f.xci.RData"
#  xinact_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact.tsv"
# zscore<-3

### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
exp_data$gene_id_red<-sapply(strsplit(exp_data$ensg,"\\."), "[[",1)


#######reduce to unique genes only
seen_genes<-dict()
i=1
while (i<=nrow((exp_data))){
   # print(i)
   this_row<-exp_data[i,]
   if (!is.na(seen_genes$get(c(this_row$ind,this_row$gene_id_red),NA))){
      #don't even need to check bc no RVs in this
      if(is.na(this_row$num_rvs)){
         i=i+1
         next
      }
      if(is.na( seen_genes[[c(this_row$ind,this_row$gene_id_red)]][1])){
         seen_genes[this_row[,1]]=c(this_row$num_rvs,i)
         i=i+1
         next
      }
      #if the max genes are recorded, continue
      if( seen_genes[[c(this_row$ind,this_row$gene_id_red)]][1]>=this_row$num_rvs){
         i=i+1
         next 
      }else{
         # else, then max
         # seen_genes[this_row[,1]]=c(this_row$num_rvs,i)
         seen_genes[[c(this_row$ind,this_row$gene_id_red)]]=c(this_row$num_rvs,i)
         i=i+1
         next
         #saves as first value the number of rare variatns, the second value is the row number
      }
   }else{
      seen_genes[[c(this_row$ind,this_row$gene_id_red)]]=c(this_row$num_rvs,i)
   }
   i=i+1
}


saveRDS(seen_genes,outfile)
