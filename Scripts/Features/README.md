# Features -- Processing RVs
Independently can be run, focuses only on processing the genetic data

### rule vcf2bed    
Goes from GTEx VCF file to a bed file for each individual with that individual's variant sites and their allele frequency      
 - input: scripts dir, euro ids (list of individuals who are european to use ), eventual outdir     
 - output: a log file, and all the files in the bed directory specified in the input   
 - scripts:  Scripts/Features/vcf2bedfiles_gnomad.sh    
