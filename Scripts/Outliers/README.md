# Call Outliers
Run this Snakefile after the Preprocessing Snakefile. 

### rule call_outliers   
This rule looks at the zscores, and based on a number of tissues for it to be a multi-tissue outlier, and the minimum zscore needed to be an outlier, calls multitissue outliers
 - input: outdir/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_{type}_{group}.txt.gz (final file of the preprocessing snakefile for type=variant type and group=m/f/both); zthresh is the minimum z score to be called an outlier, nphen is the minimum number of tissues to be called a multi-tissue outlier; outdir is the output directory for the file      
 - output: outdir./outliers_v8/outliers_both_zthresh(z_min)_nphen(nphen_min)_noglobal_medz_{type}_{group}.txt.gz     
 - scripts: Scripts/Outliers/call_outliers.R      

### rule identify_global_outliers   
This scripts identifies global outliers and removes them      
 - input: output from previous rule      
 - params: method is proportion (legacy)   
 - output: /outliers_v8/outliers_zthresh(z_min)_nphen(nphen_min)_globalOutliersRemoved_{type}_{group}.txt     
 - scripts: Scripts/Outliers/identify_global_outliers.R       


