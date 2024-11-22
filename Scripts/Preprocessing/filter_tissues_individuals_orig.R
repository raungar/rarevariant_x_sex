#!/usr/bin/env Rscript

## R script to pick tissues and individuals for imputation and outlier calling
## Run from upper level directory of the repo (or adjust path to tissue colors)

require(ggplot2)
require(gplots)
require(stringr)
require(data.table)

get_opt_parser<-function(){
        option_list = list(
                 make_option(c("-r", "--RAREDIR"), type="character", default=NULL, help="RAREDIR"),
                 make_option(c("-f", "--FILESDIR"), type="character", default=NULL, help="FILESDIR"),
                 make_option(c("-o", "--outfile_x"), type="character", default=NULL, help="output file for x only"),
                 make_option(c("-o", "--outfile_a"), type="character", default=NULL, help="output file for autosome only"),
                 make_option(c("-x", "--x_gtf_file"), type="character", default=NULL, help="path for x gtf file""),
                 make_option(c("-a", "--a_gtf_file"), type="character", default=NULL, help="path for autosomal gtf file"),
                 make_option(c("-s", "--sample_file"), type="character", default=NULL, help="sample file (gtex_2017-06-05_v8_samples_tissues.txt)"),
                 make_option(c("-n", "--norm_expr_file"), type="character", default=NULL, help="norm expr file (gtex_2017-06-05_normalized_expression.txt.gz)")
                 )

        opt_parser = OptionParser(option_list=option_list)
        return(opt_parser)
}


opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
baseDir<-as.character(args$RAREDIR)
filesDir<-as.character(args$FILESDIR)
outfile_x<-as.character(args$outfile_x)
outfile_a<-as.character(args$outfile_a)
sample_file<-as.character(args$sample_file)
norm_expr_file<-as.character(args$norm_expr_file)
x_gtf_file<-as.character(args$x_gtf_file)
a_gtf_file<-as.character(args$a_gtf_file)

#baseDir = Sys.getenv('RAREDIR')
#filesDir = Sys.getenv('FILESDIR')

dir = paste0(baseDir, '/preprocessing_v8')


##---------------- MAIN
## Read in list of EA samples
#eas.wgs = scan(paste0(filesDir, '/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt'), what = character())

## Read in sample to tissue correspondence and turn it into a individual to tissue correspondence
meta = read.table(sample_file, header = F,
                  stringsAsFactors = F, col.names = c('Sample', 'Tissue'))
meta$Id = apply(str_split_fixed(meta$Sample, '-', 6)[, c(1:2)], 1, paste, collapse = '-')

## Read in normalized expression data and subset above correspondance to individuals with corrected data
## (this is a subset because it only includes individuals that were genotyped)
#expr = read.table(gzfile(paste0(dir,'/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.txt.gz')), header=T)
expr = read.table(gzfile(norm_expr_file), header=T)
colnames(expr) = gsub("[.]", "-", colnames(expr))
meta = meta[meta$Id %in% colnames(expr), ]

## Get tissue and sample names
tissues = sort(unique(meta$Tissue))
brain.indices = grep('Brain', tissues)
individuals = unique(meta$Id)

## Make matrix of samples by tissues with expression data
exp.design = matrix(0, ncol = length(individuals), nrow = length(tissues),
                    dimnames = list(tissues, individuals))
for (t in tissues) {
    inds = meta$Id[meta$Tissue == t]
    exp.design[t, inds] = 1
}

missingness.ind = plot.ind.miss(exp.design, 'Missingness by individual (original)', thresh = 0.75)
missingness.tiss = plot.tissue.miss(exp.design, 'Missingness by tissue (original)', thresh = 0.75)
#
print(missingness.ind$plot)
print(missingness.tiss$plot)
#
### Filter for tissues with at most 75% missingness
tissues.final = as.character(missingness.tiss$keep)
exp.design = exp.design[tissues.final, ]

## Filter for individuals with at most 75% missingness (based on subset set of tissues)
inds.final = as.character(plot.ind.miss(exp.design, thresh = 0.75)$keep)
exp.design = exp.design[, inds.final]

### Make heatmap and barplots with filtered tissues and individuals
col.ind = ifelse(colnames(exp.design) %in% eas.wgs, 'navyblue', 'white') # color column if individual has WGS
#heatmap.2(exp.design, Rowv = T, Colv = T, dendrogram = 'n', labCol = '', labRow = '',
#          xlab = 'Individuals', ylab = 'Tissues', margins = c(2,2), lwid=c(0.1,4), lhei=c(0.7,4),
#          ColSideColors = col.ind,
#          RowSideColors = gtex.colors[tissues.final], col = c('white', 'dodgerblue3'), 
#          trace = 'none', cexRow = .5, cexCol = .3, key = F, main = 'GTEx design matrix, filtered')

#print(plot.ind.miss(exp.design, 'Missingness by individual (filtered)')$plot)
#print(plot.tissue.miss(exp.design, 'Missingness by tissue (filtered)')$plot)

#dev.off()

## Calculate overall missingess in the filtered design matrix
cat('Average missingness', mean(exp.design), '\n')

## Write out selected tissues and individuals
write.table(sort(inds.final), paste0(dir, '/gtex_2017-06-05_v8_individuals_passed.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(sort(tissues.final), paste0(dir, '/gtex_2017-06-05_v8_tissues_passed.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(exp.design, paste0(dir, '/gtex_2017-06-05_v8_design_passed.txt'),
            sep = '\t', quote = F, col.names = T, row.names = T)

## Subset the normalized expression file to the individuals and tissues selected
expr.subset = expr[which(expr$Tissue %in% tissues.final),]
rm(expr)
expr.subset = expr.subset[, c('Tissue', 'Gene', sort(inds.final))]

## Then subset to the genes expressed in each tissue that are either protein coding or lincRNA
## first reading in the list of autosomal protein-coding  & lincRNA genes
#autosomal.df = read.table(paste0(filesDir, '/preprocessing_v8/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt'),
autosomal.df = read.table(a_gtf_file,
                          header = F, stringsAsFactors = F)
autosomal.selected = autosomal.df[autosomal.df[, 2] %in% c('lincRNA', 'protein_coding'), 1]
x.df = read.table(x_gtf_file,
                          header = F, stringsAsFactors = F)
x.selected = x.df[x.df[, 2] %in% c('lincRNA', 'protein_coding'), 1]

print(length(unique(expr.subset$Gene)))
genes.counts = table(expr.subset$Gene)
genes.keep = names(genes.counts)[which(genes.counts == length(tissues.final))]

genes.keep.a = genes.keep[which(genes.keep %in% autosomal.selected)]
print(length(genes.keep.a))
expr.subset.a = expr.subset[which(expr.subset$Gene %in% genes.keep.a), ]
## finally restandardize and output subsetted expression matrix
expr.subset.a = as.data.frame(expr.subset.a)
expr.subset.a[, 3:ncol(expr.subset.a)] = t(scale(t(expr.subset.a[, 3:ncol(expr.subset.a)])))


genes.keep.x = genes.keep[which(genes.keep %in% x.selected)]
print(length(genes.keep.x))
expr.subset.x = expr.subset[which(expr.subset$Gene %in% genes.keep.x), ]
## finally restandardize and output subsetted expression matrix
expr.subset.x = as.data.frame(expr.subset.x)
expr.subset.x[, 3:ncol(expr.subset.x)] = t(scale(t(expr.subset.x[, 3:ncol(expr.subset.x)])))


#gzout = gzfile(paste0(dir, '/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.subset.txt.gz'), 'w')
gzout_a=gzfile(outfile_a,'w')
write.table(expr.subset.a, gzout_a, sep = '\t', quote = F, col.names = T, row.names = F)
close(gzout_a)

gzout_x=gzfile(outfile_x,'w')
write.table(expr.subset.x, gzout_x, sep = '\t', quote = F, col.names = T, row.names = F)
close(gzout_x)
