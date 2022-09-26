# sexDEG analysis
The point of this analysis is to evaluate the interplay of rare variants with sexDEGs. The sexDEG files come from the 2020 GTEx Sex paper (full files obtained from Daniel Cotter), and were separated to create one file per sheet. These are in the Files/ directory for further use in this analysis.       


#### rule get_genes
purpose: get all the genes that were used     
input: preprocessing_v8 final output file (both sexes)     
output: sexdeg_v8/all_genes.txt      

#### rule get_chr
purpose: annotate these genes with the correct chromosome     
script: Scripts/SexDeg/gtf_to_genes_chr.py     
input: sexdeg_v8/all_genes.txt, list of gtex genes     
output: sexdeg_v8/all_genes_chr.txt   

#### rule overlap_genes_sexdegs
purpose: combine the all genes file with the tissue sexDEG file (this snakefile will run in parallel at this point for each tissue)      
script: Scripts/SexDeg/combine_allgenes_sexdegs.sh    
input: sexdeg_v8/all_genes_chr.txt Files/SexDEGs/{tissue}.csv      
output: sexdeg_v8/all_genes_sexDEGs_{tissue}.txt        

#### rule add_RVs     
purpose: finally annotate with rare variants. this collapses       
input: sexdeg_v8/all_genes_sexDEGs_{tissue}.txt, from features final output features_v8/collapsed_maf_both_linc_prot_x.tsv.gz     
script: Scripts/SexDeg/intersect_sexdegs_rvs.py    
params: min and maf MAF to consider (for bins)   
output: sexdeg_v8/all_genes_sexDEGs_rvs_{tissue}.txt.gz     
log: this contains the genes that did not find an intersection (RVs with genes that were not in this list) /sexdeg_v8/all_genes_sexDEGs_rvs_{tissue}.log     

#### rule get_frac_file
purpose: count the number of sexDEGs with/without RVs, number of non-sexDEGs with/without RVs to get RR in the following script       
input: sexdeg_v8/all_genes_sexDEGs_rvs_{tissue}.txt.gz     
script: SexDeg/get_rr.sh    
output: sexdeg_v8/risktable_{tissue}.txt    

#### rule calc_risk   
purpose: use R to get upper/lower bounds and RR    
script: Scripts/SexDeg/get_rr_epitab.R   
input: sexdeg_v8/risktable_{tissue}.txt    
output: sexdeg_v8/rr_{tissue}.txt      


Finally, you can use the Rscript Scripts/SexDeg/plot_risks.R to make the plots!        
