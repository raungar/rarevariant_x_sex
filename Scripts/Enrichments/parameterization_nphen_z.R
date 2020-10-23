library("ggplot2")
library("dplyr")
library("data.table")

#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8"
groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress
chr_types=c("x") #,"aut") #aut
nphen=c(2,3,4,5)
z=c(2,2.5,3)
maf_min=0.01

risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(),chr=character(), sex=character(),zmin=numeric(), nphen_min=numeric(),num_outliers=numeric())


#read in all combination of CSV files, and assign to the appropriately named variables
for (this_sex in groups){
  for(this_chr in chr_types){
    for(this_nphen in nphen){
     for(this_z in z){
        
        file=paste0(mydir,"/outliers_zthresh",this_z,"_nphen",this_nphen,"_noglobal_medz_varAnnot_",
                    this_chr,"_",this_sex,".txt")
       varname<-paste0(this_chr,"_",this_sex,"_z",this_z,"_nphen",this_nphen)
       print(varname)
        combined_genes_vars=fread(file,sep="\t") %>% dplyr::filter(vartype=="SNP" | is.na(vartype)) ###limit to only a SNP analysis
        assign(varname,combined_genes_vars)
        
         all_outliers<-combined_genes_vars %>% dplyr::filter(Y=="outlier")
         all_controls<-combined_genes_vars %>% dplyr::filter(Y=="control")
         
         number_outliers<-nrow(all_outliers) # this gets the number of outliers
         
         exp_nn_all=nrow(all_controls  %>% dplyr::filter(is.na(gnomad_maf_both) | gnomad_maf_both > maf_min)) # controls without RV
         exp_ny_all=nrow(all_controls %>% dplyr::filter(!is.na(gnomad_maf_both) & gnomad_maf_both < maf_min)) # control w RV
         exp_yn_all=nrow(all_outliers %>% dplyr::filter(is.na(gnomad_maf_both) | gnomad_maf_both > maf_min)) ### outliers without RV
         exp_yy_all=nrow(all_outliers %>% dplyr::filter(!is.na(gnomad_maf_both) & gnomad_maf_both < maf_min)) # outliers w RV
         exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
         err_all = epitab(exptable_all, method = 'riskratio')
         
         risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                         Lower = err_all$tab[2,6],
                                         Upper = err_all$tab[2,7],
                                         Pval = err_all$tab[2,8],
                                         Cat = varname,
                                         chr=this_chr,
                                         sex=this_sex,
                                         zmin=this_z,
                                         nphen_min=this_nphen,
                                         num_outliers=number_outliers))  
        
      }
    }
  }
}

#filter to outliers and controls

risks$CI<-risks$Upper-risks$Lower


risks = risks %>% arrange(by=Risk) 

