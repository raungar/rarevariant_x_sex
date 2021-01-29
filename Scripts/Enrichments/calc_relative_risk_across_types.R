#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--out_rdata_relative"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--out_rdata_continuous"), type = 'character', default = NULL, help = "path of output file (RDATA) continuous risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
                make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
                make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
out_rdata_relative <- as.character(opt$out_rdata_relative)
out_rdata_continuous <- as.character(opt$out_rdata_continuous)
zscore <- as.numeric(opt$zscore)
min_maf<-as.numeric(opt$min_maf)
max_maf<-as.numeric(opt$max_maf)

# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outlier_noglobal_medz_varAnnot_zthresh2_nphen2_aut_m_linc_prot.txt.gz"
# #infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen5_noglobal_medz_varAnnot_x_m.txt"
# zscore<-3
# min_maf<-0
# max_maf<-0.01

if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
#exp_data = exp_data %>% select(indiv_id,gene_id,MedZ,Y,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
#exp_data = exp_data %>% dplyr::select(ind,ensg,MedZ,Y,af,isOutlier,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
#exp_data = dplyr::filter(exp_data,sv_v7==1)
#new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$variant_cat[x] == 'splice', exp_data$tier2[x], exp_data$variant_cat[x]))
#exp_data$variant_cat = ind	ensg	N	Df	MedZ	Y	chr	pos	vartype	sex	ref	alt	varswitch	gtex_maf	gnomad_maf_both	gnomad_maf_m	gnomad_maf_f	genetype	gnomad_maf_diff	num_rvs
colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","chr","pos","vartype","sex","ref","alt","bothvars",
                     "gtex_maf","gnomad_maf_both","gnomad_maf_m","gnomad_maf_f","genetype","mafdiff","numrv")
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
  this_gnomad_both_maf<-as.numeric(x[which(colnames(exp_data)=="gnomad_maf_both")])
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_gnomad_both_maf)){"common"}
  else if (this_gnomad_both_maf>=min_maf & this_gnomad_both_maf<max_maf){"rare"}
  else if (this_gnomad_both_maf>=max_maf){"common"}
  else (stop("ERROR: VARIANT NOT MAKING SENSE"))
})
  


   
   
print("filter for z score")
### get relative risk per category
exp_data$has_variant<-has_variant
exp_data$variant_cat<-exp_data$genetype ##########CHANGE THIS LINE FOR WHAT U WANT THE RR TO BE....

exp_outliers = dplyr::filter(exp_data, abs(MedZ) >= zscore)
exp_controls = dplyr::filter(exp_data, abs(MedZ) < zscore)




print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
vcats = na.omit(unique(exp_data$variant_cat))
for (vcat in (vcats)) {
   print(vcat)
   exp_nn = nrow(exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter( has_variant != "rare"))
   exp_ny = nrow(exp_controls %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))
   exp_yn = nrow(exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant != "rare"))
   exp_yy = nrow(exp_outliers %>% dplyr::filter(variant_cat == vcat) %>% dplyr::filter(has_variant == "rare"))
   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
   print(exptable)
   err = epitab(exptable, method = 'riskratio')
   risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   Cat = vcat,
                                   Type = 'Total expression'))   
}

exp_nn_all=nrow(dplyr::filter(exp_controls,has_variant != "rare"))
exp_ny_all=nrow(dplyr::filter(exp_controls,has_variant == "rare"))
exp_yn_all=nrow(dplyr::filter(exp_outliers,has_variant != "rare"))
exp_yy_all=nrow(dplyr::filter(exp_outliers,has_variant == "rare"))
exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
err_all = epitab(exptable_all, method = 'riskratio')
risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
                                Lower = err_all$tab[2,6],
                                Upper = err_all$tab[2,7],
                                Pval = err_all$tab[2,8],
                                Cat = "all",
                                Type = 'Total expression'))   
risks = risks %>% arrange(by=Risk) 
risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))

save(risks,  file=out_rdata_relative)
print(paste0("COMPLETE: ",out_rdata_relative))
# save(all_coefs, file=out_rdata_continuous)


# 
# print("continuous risk")
# ### Continuous risk
# all_coefs = data.frame(Beta = numeric(), SE = numeric(), Zval = numeric(), Pval = numeric(), Category = character(), Variant = character())
# 
# for (vcat in vcats) {
#   print(vcat)
#   
#   exp_data = exp_data %>% mutate(HasCAT = ifelse(variant_cat == vcat, 1, 0))
#   exp_lm = as.data.frame(summary(glm(HasCAT ~ OutlierValue, data = exp_data, family = binomial))$coefficients)[2,]
#   colnames(exp_lm) = c('Beta', 'SE', 'Zval', 'Pval')
#   new_coefs=exp_lm %>% mutate(Category = 'eOutliers', Variant = vcat)
#   all_coefs = rbind(all_coefs, new_coefs)
# }
# 

#pcols = exp_data %>% filter(!(medz_bin %in% c(NA,"0.000~0.002"))) %>%
#  group_by(variant_cat) %>% sample_n(1) %>%
 # ungroup() %>% mutate(CatCol = rgb(color_R, color_G, color_B)) %>%
#  select(variant_cat, CatCol)
#plot_cols = c(pcols$CatCol[1:11],'#3e6690','#33669a',pcols$CatCol[14:16])
#names(plot_cols) = unique(pcols$variant_cat)


