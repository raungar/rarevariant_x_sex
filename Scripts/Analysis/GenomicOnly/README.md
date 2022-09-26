# Genomic Analysis

This section is for analyses that do not involve RNA-seq    

This will look at the number of rare variants on the X-chromosome and autosomes and investigate the patterns of variant type for males and females and how certain regions of the X matter

### Snakefile Rules   
##### rule get_chr_rvs
purpose: generate for each chr1:22,X the number of rare variants and significance of difference between males vs. females for each MAF and variant type          
input: enrichment_v8/[x/aut]_all_rvs_inds_types.txt.gz (chr/start/end/maf/vartype/ind), eurofile, chr len file, file of samples sex      
output: number of RVs for MAF x chr x variant type x subtype x sex, then significance of differences for this combination between M vs. F                     

##### rule  get_x_subtypes
purpose: Get subtypes (xar, xcr1, xcr2, xtr, par1, par2, nonpar) and patterns across them. This will return two files, the number of rvs and significance of difference between males and females. The significance of difference is computed by a BH-adjusted wilcox.test()      
input: x_all_rvs_inds_types.txt.gz, eurofile, file of sample's sex        
output: number of rvs for MAF x subregion x variant type x subtype x sex, then significance of difference of this combination for males and females     


##### rule  get_maf_distribution
purpose: get a chosen MAF where there is no summary, all individuals at that MAF w specified sex to plot distributions later
input: enrichment_v8/[x/aut]_all_rvs_inds_types.txt.gz (chr/start/end/maf/vartype/ind), eurofile, chr len file, file of samples sex, maf
output: for a maf, sex, all individuals, all vartypes in a txt.gz file


### Plotting Files    
#### plots_x_summaries.R
purpose: combines previous analysis, manually edit to plot #RVs/bp*10000 vs MAF and -log(p) vs MAF       
