## R script to pick tissues and individuals for imputation and outlier calling
## Run from upper level directory of the repo (or adjust path to tissue colors)
print("RUN FILTER TISSUES")
library(stringr)
library(data.table)
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c('--RAREDIR'), type='character', default=NULL, help='RAREDIR'),
                  make_option(c("--outfile_x"), type="character", default=NULL, help="output file for x only"),
                  make_option(c("--outfile_a"), type="character", default=NULL, help="output file for autosome only"),
                  make_option(c("--x_gtf_file"), type="character", default=NULL, help="path for x gtf file"),
                  make_option(c("--a_gtf_file"), type="character", default=NULL, help="path for autosomal gtf file"),
                  make_option(c("--group"), type="character", default=NULL, help="group: both, m, f"),
                  make_option(c("--sample_file"), type="character", default=NULL, help="sample file like gtex_2017-06-05_v8_samples_tissues.txt"),
                  make_option(c("--norm_expr_file"), type="character", default=NULL, help="norm expr file like gtex normalized")
)

# option_list = list(
#   make_option(c('--RAREDIR'), type='character', default=NULL, help='RAREDIR')
# )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("NOW?")
# print(opt_parser)
print("PARSED")
dir<-as.character(opt$RAREDIR)
print(dir)
outfile_x<-as.character(opt$outfile_x)
outfile_a<-as.character(opt$outfile_a)
my_group<-as.character(opt$group)
print(my_group)
sample_file<-as.character(opt$sample_file)
print(paste0(sample_file,"IS THE SAMPLE FILE"))
norm_expr_file<-as.character(opt$norm_expr_file)
x_gtf_file<-as.character(opt$x_gtf_file)
a_gtf_file<-as.character(opt$a_gtf_file)

# if(TRUE){stop("STOP ME")}

print(norm_expr_file)

# RAREDIR="/Volumes/groups/smontgom/raungar/Sex/Output"
# sample_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt"
# norm_expr_file="/Volumes/groups/smontgom/raungar/Sex/Output/expression_v8/Counts/normalized_expression_gathered-f-protect_sex-nullshuffled_v8-both-binary.txt.gz"
# x_gtf_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/x_proteincoding_lncrna.gtf"
# a_gtf_file="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_proteincoding_lncrna.gtf"
# outfile_x="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subset.x.both.txt.gz"
# outfile_a="/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subset.aut.both.txt.gz"
# my_group="m"
# 


#baseDir = Sys.getenv('RAREDIR')
#filesDir = Sys.getenv('FILESDIR')

#dir = paste0(baseDir, '/preprocessing_v8')

##---------------- FUNCTIONS
plot.ind.miss <- function(design, title = '', thresh = NULL) {
    ind.miss = data.frame(ind = colnames(design), perc = 1 - colMeans(design))
    ind.miss = ind.miss[order(ind.miss$perc), ]
    ind.miss.bar = ggplot(data = ind.miss, aes(x = c(1:nrow(ind.miss)), y = perc)) +
        geom_bar(stat = 'identity', fill = 'dodgerblue3')  +
        xlab('Individuals') + ylab('Missingness') + ggtitle(title)
    if (!is.null(thresh)) {
        ind.miss.bar = ind.miss.bar + geom_hline(yintercept = thresh)
        keep = ind.miss$ind[ind.miss$perc <= thresh]
    } else {
        keep = ind.miss$ind
    }
    return(list(plot = ind.miss.bar, keep = keep))
}

plot.tissue.miss <- function(design, title = '', thresh = NULL) {
    tiss.miss = data.frame(tissue = rownames(design), perc = 1 - rowMeans(design))
    tiss.miss = tiss.miss[order(tiss.miss$perc), ]
    tiss.miss$tissue = factor(tiss.miss$tissue, levels = tiss.miss$tissue)
    tiss.miss.bar = ggplot(data = tiss.miss, aes(x = tissue, y = perc, fill = tissue)) +
        geom_bar(stat = 'identity') + coord_flip() + guides(fill = F) +
        xlab('') + ylab('Missingness') + ggtitle(title) #+
       # scale_fill_manual(values = gtex.colors) 
    if (!is.null(thresh)) {
        tiss.miss.bar = tiss.miss.bar + geom_hline(yintercept = thresh)
        keep = tiss.miss$tissue[tiss.miss$perc <= thresh]
    } else {
        keep = tiss.miss$tissue
    }
    return(list(plot = tiss.miss.bar, keep = keep))
}



##---------------- MAIN
## Read in list of EA samples
#eas.wgs = scan(paste0(filesDir, '/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt'), what = character())

## Read in sample to tissue correspondence and turn it into a individual to tissue correspondence
#meta = read.table(sample_file, header = F,
print("SAMPLE FILE")
print(sample_file)
meta = fread(sample_file, header = F,
                  stringsAsFactors = F, col.names = c('Sample', 'Tissue'))
meta$Id = apply(str_split_fixed(meta$Sample, '-', 6)[, c(1:2)], 1, paste, collapse = '-')

## Read in normalized expression data and subset above correspondance to individuals with corrected data
## (this is a subset because it only includes individuals that were genotyped)
#expr = read.table(gzfile(paste0(dir,'/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.txt.gz')), header=T)
#expr = read.table(gzfile(norm_expr_file), header=T)
print("EXPR") #Adipose_Subcu
expr = fread(norm_expr_file, header=T,fill = TRUE)
#expr2 = fread(norm_expr_file, header=T)
colnames(expr)[1:2]<-c('Gene','Tissue')

print(head(expr[1:4,1:4]))
#expr=expr%>%select(-Id)
tissue_dic<-as.character(sort(unique(meta$Tissue)))
print("tissue dic pre becoming a dictionary so it's just vals")
print(tissue_dic)
print("len meta than expr tissue")
print(length(unique(meta$Tissue)))
print(length(unique(expr$tissue)))
expr$Tissue<-sapply(strsplit(expr$Tissue,"\\."),"[[",1)
print(head(expr$Tissue))

names(tissue_dic)<-as.character(sort(unique(expr$Tissue)))
print("tissue_dic!!")
print(tissue_dic)
print("now printing unique(meta$Tissue)")
print(unique(meta$Tissue))
print("now printing unique(expr$Tissue)")
print(unique(expr$Tissue))
######expr$Tissue<-tissue_dic[expr$Tissue]
#print(head(expr))
print("expr tissues....")
print(unique(expr$Tissue))

colnames(expr) = gsub("[.]", "-", colnames(expr))
#print(head(colnames(expr)))
meta = meta[meta$Id %in% colnames(expr), ]

## Get tissue and sample names
tissues = sort(unique(meta$Tissue))
print("TISSUES FROM META")
print(tissues)
brain.indices = grep('Brain', tissues)
individuals = unique(meta$Id)

## Make matrix of samples by tissues with expression data
exp.design = matrix(0, ncol = length(individuals), nrow = length(tissues),
                    dimnames = list(tissues, individuals))
for (t in tissues) {
    inds = meta$Id[meta$Tissue == t]
    exp.design[t, inds] = 1
}


#missingness.ind = plot.ind.miss(exp.design, 'Missingness by individual (original)', thresh = NULL)
#missingness.tiss = plot.tissue.miss(exp.design, 'Missingness by tissue (original)', thresh = NULL)
#
#print(missingness.ind$plot)
#print(missingness.tiss$plot)
#
### Filter for tissues with at most 75% missingness
##EDIT PREVENT FILTERING
#tiss.miss = data.frame(tissue = rownames(exp.design), perc = 1 - rowMeans(exp.design))
#tiss.miss = tiss.miss[order(tiss.miss$perc), ]
#tiss.miss$tissue = factor(tiss.miss$tissue, levels = tiss.miss$tissue)
#tissues.final=as.character(tiss.miss$tissue)
tissues.final=rownames(exp.design)
print("TISSUES FINAL")
print(tissues.final)

#tissues.final = as.character(missingness.tiss$keep)
exp.design = exp.design[tissues.final, ]


#print("begin ind.miss")
## Filter for individuals with at most 75% missingness (based on subset set of tissues)
#inds.final = as.character(plot.ind.miss(exp.design, thresh = NULL)$keep)
ind.miss = data.frame(ind = colnames(exp.design), perc = 1 - colMeans(exp.design))
ind.miss = ind.miss[order(ind.miss$perc), ]
inds.final = as.character(ind.miss$ind)
exp.design = exp.design[, inds.final]
print("INDS.FINAL")
#print(head(inds.final))

### Make heatmap and barplots with filtered tissues and individuals
cat('Average missingness', mean(exp.design), '\n')

## Write out selected tissues and individuals
write.table(sort(inds.final), paste0(dir, '/gtex_2017-06-05_v8_individuals_passed_',my_group,'.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(sort(tissues.final), paste0(dir, '/gtex_2017-06-05_v8_tissues_passed_',my_group,'.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(exp.design, paste0(dir, '/gtex_2017-06-05_v8_design_passed_',my_group,'.txt'),
            sep = '\t', quote = F, col.names = T, row.names = T)

## Subset the normalized expression file to the individuals and tissues selected
print("EXPR TISSUE IN FINAL")
#print(expr$Tissue)
print(head(expr))
print(unique(expr$Tissue))
#print(unique(tissues.final))
print(head(expr$Tissue))
print(head(tissues.final))
expr.subset = expr[which(expr$Tissue %in% tissues.final),]

print("expr subset tissue")
print(head(expr.subset))
print(unique(expr.subset$Tissue))
#head(expr.subset)
# rm(expr)
#expr.subset = expr.subset[, c('Tissue', 'Gene', sort(get(inds.final)))]
#print("EXPR>SUBSET")
#print(head(expr.subset))
## Then subset to the genes expressed in each tissue that are either protein coding or lincRNA
## first reading in the list of autosomal protein-coding  & lincRNA genes

print("Read gtf")

autosomal.df = fread(a_gtf_file,
                          header = F, stringsAsFactors = F)
colnames(autosomal.df)<-c("gene","type")
#print("autosomal.df")
#print(head(autosomal.df))
#autosomal.selected = autosomal.df[autosomal.df[, 2] %in% c('lincRNA', 'protein_coding'), 1]
#print("autosomal.selected")
#print(head(autosomal.selected))
autosomal.selected<-autosomal.df
##autosomal.selected$gene<-sapply(strsplit(autosomal.selected$gene,"\\."),"[[",1)
#print(head(autosomal.selected))

print("read x")

x.df = fread(x_gtf_file,
                          header = F, stringsAsFactors = F)
colnames(x.df)<-c("gene","type")
#print("x.df")
#print(head(x.df))
x.selected<-unique(x.df)
####x.selected$gene<-sapply(strsplit(x.selected$gene,"\\."),"[[",1)
print(head(x.selected))
#x.selected = x.df[x.df[, 2] %in% c('lincRNA', 'protein_coding'), 1]
#x.selected<-sapply(strsplit(x.selected,"\\."),"[[",1)


#print(length(unique(expr.subset$Gene)))
genes.counts = table(expr.subset$Gene)
#print(genes.counts)
#genes.keep = names(genes.counts)[which(genes.counts == length(tissues.final))]
print(" head expr subset plz")
print(head(expr.subset[1:4,1:5]))
genes.keep=expr.subset$Gene

print("GENES KEEP")
print(head(genes.keep))
print("selected")
print(head(x.selected$gene))
#print(length(genes.keep))
#print(length(genes.keep))
#print("AUTOSOMAL SELEcTED")
#print((unique(autosomal.selected)))
#print("x selected")
#print((unique(x.selected)))
#print("table autosome then x")
#print(table(genes.keep %in% autosomal.selected$gene))
#print(table(autosomal.selected %in% genes.keep))
#print(table(genes.keep %in% x.selected$gene))
print(table(x.selected$gene %in% genes.keep))
#print(paste0(head(genes.keep), " IN a: ", head(autosomal.selected$gene)))
#print(paste0(head(genes.keep), " IN x: ", head(x.selected$gene)))

genes.keep.a = genes.keep[which(genes.keep %in% autosomal.selected$gene)]
print(length(genes.keep.a))
print("subset 1")
print(expr.subset[1:4,1:4])
expr.subset.a = expr.subset[which(expr.subset$Gene %in% genes.keep.a), ]
print("subset 2")
print(expr.subset.a[1:4,1:4])
## finally restandardize and output subsetted expression matrix
expr.subset.a = as.data.frame(expr.subset.a)
print("Subset 3")
print(expr.subset.a[1:6,1:6])
expr.subset.a[, 3:ncol(expr.subset.a)] = t(scale(t(expr.subset.a[, 3:ncol(expr.subset.a)])))
print("HIIIIII")
print("HERE")
print(length(x.selected$gene))
genes.keep.x = genes.keep[which(genes.keep %in% x.selected$gene)]
print(length(genes.keep.x))
expr.subset.x = expr.subset[which(expr.subset$Gene %in% genes.keep.x), ]
## finally restandardize and output subsetted expression matrix
print("now subset")
expr.subset.x = as.data.frame(expr.subset.x)
print("3:ncol")
head(expr.subset.x)
expr.subset.x[, 3:ncol(expr.subset.x)] = t(scale(t(expr.subset.x[, 3:ncol(expr.subset.x)])))


print("writing aut")
#gzout = gzfile(paste0(dir, '/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.subset.txt.gz'), 'w')
gzout_a=gzfile(outfile_a,'w')
write.table(expr.subset.a, gzout_a, sep = '\t', quote = F, col.names = T, row.names = F)
close(gzout_a)

print("writing x")
gzout_x=gzfile(outfile_x,'w')
write.table(expr.subset.x, gzout_x, sep = '\t', quote = F, col.names = T, row.names = F)
close(gzout_x)

print("complete")



