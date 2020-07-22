# Genomic Analysis

This section is for analyses that do not involve RNA-seq    

This will look at the number of rare variants on the X-chromosome and autosomes and investigate the patterns of variant type for males and females and how certain regions of the X matter


#### get_chr_rvs.R
purpose: output cumulative rv across MAF (numrv), test if sig difference between males and females (sigdif)  
input: all_rvs_inds_types.txt.gz from enrichment script    
output: [x/aut]_all_numrv.txt.gz, [x/aut]_all_sigdif.txt.gz    

#### get_xsubtypes_rv.R 
purpose: output cumulative rv across MAF (numrv) for all subtypes (par1, par2, nonpar, xcr1, xcr2, xtr, xar), test if sig difference between males and females (sigdif)  
input: x_all_rvs_inds_types.txt.gz from enrichment script    
output: x_subtype_all_numrv.txt.gz, x_subtype_all_sigdif.txt.gz    
