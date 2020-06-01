# Preprocessing Scripts

### rule get_samples_tissues
This creates a covariates file of sex and pcs    
 - output: outdir/preprocessing_v8/covariates.txt    
 - scripts: Scripts/Preprocessing/preprocess_expr.R    

###  rule split_tissues_reads   
This splits the tissues into reads files   
  - output: outdir/preprocessing_v8/peerdir/[tissue].reads.txt     
  - scripts: Scripts/Preprocessing/split_expr_by_tissues.py  

### split_tissues_tpm
This splits the tissues into tpm files   
  - output: otudir/preprocessing_v8/peerdir/[tissue].tpm.reads.txt    
  - scripts: Scripts/Preprocessing/split_expr_by_tissues.py  

### rule split_tissues_reads
This subsets the files into males/females/both into equal numbers based on individuals in the most tissues   
Then TPM restrictions are applied and reads are transformed    
 - output: outdir/preprocessing_v8/peerdir/[tissue].[m/f/both].log2.ztrans.txt, outdir/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt   
 - scripts: Scripts/Preprocessing/preprocess_expr.R   

### rule get_eqtl_genotypes  
This prepares a smattering of proper files for process_eqtl_genotypes    
 - output: outdir/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs, outdir/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf,    
		outdir/preprocessing_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_sorted.vcf.gz
 - scripts: Scripts/Preprocessing/get_eqtl_genotypes.sh    


### rule process_eqtl_genotypes
This creates an eQTL genotypes file for regression in PEER   
  - output: outdir/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt     
  - scripts: Scripts/Preprocessing/process_gtex_v8_cis_eqtl_genotypes.R   
 

### rule get_peer_factors
This calculates PEER factors for all tissues and groups, creating a new PEER directory for each   
  - output: outdir/preprocessing_v8/peerdir/[tissue]/factors.tsv (and log files)      
  - scripts: Scripts/Preprocessing/calculate_PEER_factors.sh, Scripts/Preprocessing/calculate_PEER_factors.R    

### rule get_resids_ztrans
This uses eqtl files,metadata file, covariates (pcs) files, and peer factors to calculate residuals.      
Parameters include factors_type which allows for factors_sexregress or factors to choose which PEER .tsv file to use    
as well as provide a logfile parameter to output logfile in tissue PEER directory   
and includes a inclsex parameter to choose whether to include sex in the PEER regression. 
At last, the residuals are calculated from a linear model including SEX if desired, PC1:3, PEER.tsv, and eQTL file.   
These residuals are then z transformed for use for outlier analysis.     
  - output: log2.ztrans.[m/f/both/both.sex] and associated logfiles
  - scripts: Scripts/Preprocessing/calculate_PEER_residuals.sh, Scripts/Preprocessing/calculate_PEER_residuals.R

### peer_correct_sex
This regresses out sex from peer factors   
  - output: outdir/preprocessing_v8/peerdir/[tissue]Both/factors_sexregress.tsv    
  - scripts: Scripts/Preprocessing/get_residuals_peer_by_sex.sh, Scripts/Preprocessing/get_residuals_peer_by_sex.R

### get_resids_ztrans_sex
This uses same structure and process as get_resids_ztrans with different parameters    
Residuals but include sex to regres out (same way in ferraro et al paper)   
  - output: log2.ztrans.both.sex.txt and associated logfiles    
  - scripts: Scripts/Preprocessing/calculate_PEER_residuals.sh, Scripts/Preprocessing/calculate_PEER_residuals.R   

### get_resids_ztrans_regresssex
This uses same structure and process as get_resids_ztrans with different parameters   
Residuals both of PEER factors where sex is regressed out (note: different than above!)    
  - output: log2.ztrans.both,sexregress.txt
  - scripts: Scripts/Preprocessing/calculate_PEER_residuals.sh, Scripts/Preprocessing/calculate_PEER_residuals.R    

### get_tissue_by_individual
makes a file that has each tissue and each individual in this tissue   
  - output: gtex_2017-06-05_tissues_all_normalized_samples_[group].txt, gtex_2017-06-05_individuals_all_normalized_samples_[group].txt   
  - scripts: Scripts/Preprocessing/get_tissue_by_individual.sh   

### generate_gtf
gets gtf files for autosomes and x chromosome separately only protein-coding and lncRNA   
  - output: outdir/preprocessing_v8/[x/autosomal]_proteincoding_lncrna.gtf   
  - scripts: Scripts/Preprocessing/subset_gtf.sh    

### gather_filter_norm_expr
combines all the individual tissue file into one file for each group    
  - output: outdir/preprocessing_v8/gtex_2017-06-05_normalized_expression_{group}.txt.gz    
  - scripts: Preprocessing/gather_filter_normalized_expression.py    

### filter_tissues_individuals
takes above file, and subsets to the genes that are in generated gtf files   
  - output: outdir/preprocessing_v8/gtex_2017-06-05_normalized_expression_subset_[aut/x]_[group].txt.gz
  - scripts: Scripts/Preprocessing/filter_tissues_individuals.R
