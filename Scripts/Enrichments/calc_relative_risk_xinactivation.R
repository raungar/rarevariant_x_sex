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
                make_option(c("--xinact_file"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--out_rdata_relative"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score")
        )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
xinact_file <- as.character(opt$xinact_file)
out_rdata_relative <- as.character(opt$out_rdata_relative)
zscore <- as.numeric(opt$zscore)

# zscore<-2
#  infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/outliers_zthresh3_nphen5_noglobal_medz_varAnnot_x_m.txt"
#  #out_rdata_relative<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_z3_x_f.xci.RData"
#  xinact_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact.tsv"
# zscore<-3

if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

   
### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
exp_data$gene_id_red<-sapply(strsplit(exp_data$ensg,"\\."), "[[",1)

#exp_data = exp_data %>% select(indiv_id,gene_id,MedZ.x,Y,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
# exp_data = exp_data %>% dplyr::select(indiv_id,gene_id,MedZ.x,Y,af,isOutlier,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
# exp_data_final = dplyr::filter(exp_data,sv_v7==1)
#new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$variant_cat[x] == 'splice', exp_data$tier2[x], exp_data$variant_cat[x]))
#exp_data_final$variant_cat = new_cats
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

xinact<-fread(xinact_file,data.table=F)
xinact$GeneID_red<-sapply(strsplit(xinact$`Gene ID`,"\\."), "[[",1)

exp_data<-merge(exp_data,xinact[c("Combined XCI status", "GeneID_red")], by.x="gene_id_red",by.y="GeneID_red")
colnames(exp_data)[colnames(exp_data)=="Combined XCI status"]<-"XCI_STATUS"

print("filter")
### get relative risk per category
exp_outliers = exp_data %>% dplyr::filter(abs(MedZ) >= zscore)  %>% dplyr::filter(as.numeric(Df) > 1) 
exp_controls = dplyr::filter(exp_data, abs(MedZ) < zscore)

print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
xci_status = unique(exp_data$XCI_STATUS)
for (xci in xci_status) {
   print(xci)
   exp_nn = nrow(dplyr::filter(exp_controls, XCI_STATUS != xci))
   exp_ny = nrow(dplyr::filter(exp_controls, XCI_STATUS == xci))
   exp_yn = nrow(dplyr::filter(exp_outliers, XCI_STATUS != xci))
   exp_yy = nrow(dplyr::filter(exp_outliers, XCI_STATUS == xci))
   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
   err = epitab(exptable, method = 'riskratio')
   risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   XCI_STATUS = xci,
                                   Type = 'Total expression'))   
}

# exp_nn_all=nrow(dplyr::filter(exp_controls,!(XCI_STATUS %in% xci_status)))
# exp_ny_all=nrow(dplyr::filter(exp_controls,XCI_STATUS %in% xci_status))
# exp_yn_all=nrow(dplyr::filter(exp_outliers,!(XCI_STATUS %in% xci_status)))
# exp_yy_all=nrow(dplyr::filter(exp_outliers,XCI_STATUS %in% xci_status))
# exptable_all = rbind(c(exp_nn_all,exp_ny_all),c(exp_yn_all,exp_yy_all))
# err_all = epitab(exptable_all, method = 'riskratio')
# risks = rbind(risks, data.frame(Risk = err_all$tab[2,5],
#                                 Lower = err_all$tab[2,6],
#                                 Upper = err_all$tab[2,7],
#                                 Pval = err_all$tab[2,8],
#                                 XCI_STATUS = "all",
#                                 Type = 'Total expression'))   
risks = risks %>% arrange(by=Risk) 
risks$XCI_STATUS = factor(risks$XCI_STATUS, levels=unique(risks$XCI_STATUS))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))



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

write.csv(risks,  file=out_rdata_relative)
print(paste0("COMPLETE: ",out_rdata_relative))
# save(all_coefs, file=out_rdata_continuous)

