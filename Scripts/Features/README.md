# Features -- Processing RVs
Independently can be run, focuses only on processing the genetic data

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

### rule vcf_add_genes_x     
This using the correct gtf annotates the actual genes that are within a 10000 bp window
 - input: individual anno file from previous rule output, chr specific gtf file      
 - output: plops output into GenesAnno directory within this bySiteAnnoX directory with useful information      


### rule vcf_add_genes_aut     
Agains ame as above, just in autosomal specific directories    
