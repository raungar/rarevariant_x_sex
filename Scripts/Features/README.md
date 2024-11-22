# Features -- Processing RVs
Independently can be run, focuses only on processing the genetic data into downstream files. Use  SnakefileGnomad.      

### rule vcf2bed    
Goes from GTEx VCF file to a bed file for each individual with that individual's variant sites and their allele frequency      
 - input: scripts dir, euro ids (list of individuals who are european to use ), eventual outdir     
 - output: a log file, and all the files in the bed directory specified in the input   
 - scripts:  Scripts/Features/vcf2bedfiles_gnomad.sh    

### rule vcf_anno_run_chrX      
This takes the exact variants from the sample file, and see if it exists in gnomadv3. if so, this information is recorded in a column. if there is no match to an exisiting variant, then "NO_MATCH" is recorded.     
 - input:  gnomad v3 file for annotation (individually by chr for parallelization), sample bed file from previous rule      
 - output: Output/features_v8/bySiteAnnoX for individual file     
 - scripts: none!  but bedtools module loading is necessary    

### rule vcf_anno_run_aut      
This rule is the same as above, just collects input from the autosomal directory and puts output into bySiteAnnoAut       

### rule gtf_split_by_gene:    
Split gtf to separate chromosomes for parallelization    
 - input: gtex gtf     
 - params: gtf directory to output the split gtf file into       
 - output: features_v8/GTF/combine_genes.log (and not via snakefile, but split files as features_v8/GTF/chr{chrN}.gtf)     

### rule gtf_split_by_gene_lincRNA_protcod_only
get only lincRNA and protcod genes
 - input: features_v8/GTF/chr{chrN}.gtf, features_v8/GTF/combine_genes.log (for snakemake rule logic only)     
 - output: features_v8/GTF_lnc_protcod_only/chr{chrN}.gtf    

### rule vcf_add_genes_x     
This using the correct gtf annotates the actual genes that are within a 10000 bp window. if there is no gene within 10kb, annotations it with NAs but keep variants        
 - input: individual anno file from previous rule output, chr specific gtf file (choose only lincRNA or prot coding)     
 - output: plops output into GenesAnno directory within this bySiteAnnoX directory with useful information      


### rule vcf_add_genes_aut     
Again same as above, just in autosomal specific directories     


### rule collapse_vars_to_genes_x:
This collapse the variant information to just one information per gene. This does so by taking the minimum MAF and reporting only that variant per gene. This does this for the nef_both MAF, as well as the nef_m and nef_f MAF (output into separate file). Just for fun, also writes to a file for any M vs. F MAF difference greate than cutoff_mafdiff parameter     
 - input: gene anno directory from above to read from, file wtih information about the sex of the individual       
 - params:  maf min to report it was a MAF, cutoff for mafdiff between the sexes      
 - output: maff_diff file, collapsed file for both/m/f (4 total)           
 - scripts: Scripts/Features/collapse_variants_genes.py     

### rule collapse_vars_to_genes_aut
again, just aut version (has specific script collapse_variants_genes_aut.py)         

### rule combine_rvs_inds_x
This takes the final gene anno files that are by individual, and put them into a combined file in a more manageable useful way      
final headers are: chr,start,end,maf,sample,vartype,ensg,genetype       
 - input: genes anno dir 
 - output: combined file (x_all_rvs_inds_types.txt.gz)

### rule combine_rvs_inds_aut 
same as above, aut specific directories/columns       
