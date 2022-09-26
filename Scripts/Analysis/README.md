# Analysis Files
These files are not part of the pipeline, but are rather a series of scrips used to generate figures and other one-off analyses.

#### rule get_chr_rvs
purpose: generate for each chr1:22,X the number of rare variants and significance of difference between males vs. females for each MAF and variant type          
input: enrichment_v8/[x/aut]_all_rvs_inds_types.txt.gz (chr/start/end/maf/vartype/ind), eurofile, chr len file, file of samples sex      
output: number of RVs for MAF x chr x variant type x subtype x sex, then significance of differences for this combination between M vs. F                     

#### rule get_chr_rvs_bins
purpose: same as above, just get the bins      


#### rule  get_x_subtypes
purpose: Get subtypes (xar, xcr1, xcr2, xtr, par1, par2, nonpar) and patterns across them. This will return two files, the number of rvs and significance of difference between males and females. The significance of difference is computed by a BH-adjusted wilcox.test()      
input: x_all_rvs_inds_types.txt.gz, eurofile, file of sample's sex        
output: number of rvs for MAF x subregion x variant type x subtype x sex, then significance of difference of this combination for males and females      

#### rule get_x_subtypes_bins
purpose: same as above, just get the bins    

#### rule get_maf_distribution 
purpose: prepares file to plot a maf distribution
input: enrichments_v8/{type}_all_rvs_inds_types.txt.gz, chr len file, euro inds, metadata subjects fiel       
params: maf, chrtype       
output: analysis_v8/genomic_only/{type}_maf_{maf}.txt.gz      




other scripts not in rules are to make plots or are archives     
