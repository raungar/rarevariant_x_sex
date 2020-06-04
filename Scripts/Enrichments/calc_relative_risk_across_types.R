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
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score")
        )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
out_rdata_relative <- as.character(opt$out_rdata_relative)
out_rdata_continuous <- as.character(opt$out_rdata_continuous)
zscore <- as.numeric(opt$zscore)

### expression
print('Reading expression')
exp_data = fread(infile,data.table=F)
#exp_data = exp_data %>% select(indiv_id,gene_id,MedZ.x,Y,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
exp_data = exp_data %>% select(indiv_id,gene_id,MedZ.x,Y,af,isOutlier,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
exp_data = filter(exp_data,sv_v7==1)
#new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$variant_cat[x] == 'splice', exp_data$tier2[x], exp_data$variant_cat[x]))
#exp_data$variant_cat = new_cats
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ.x)))

print("filter")
### get relative risk per category
exp_outliers = filter(exp_data, abs(MedZ.x) >= zscore)
exp_controls = filter(exp_data, abs(MedZ.x) < zscore)

print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
vcats = unique(exp_data$variant_cat)
for (vcat in vcats) {
   print(vcat)
   exp_nn = nrow(filter(exp_controls, variant_cat != vcat))
   exp_ny = nrow(filter(exp_controls, variant_cat == vcat))
   exp_yn = nrow(filter(exp_outliers, variant_cat != vcat))
   exp_yy = nrow(filter(exp_outliers, variant_cat == vcat))
   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
   err = epitab(exptable, method = 'riskratio')
   risks = rbind(risks, data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   Cat = vcat,
                                   Type = 'Total expression'))   
 }
risks = risks %>% arrange(by=Risk) 
risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))




print("continuous risk")
### Continuous risk
all_coefs = data.frame(Beta = numeric(), SE = numeric(), Zval = numeric(), Pval = numeric(), Category = character(), Variant = character())

for (vcat in vcats) {
  print(vcat)
  
  exp_data = exp_data %>% mutate(HasCAT = ifelse(variant_cat == vcat, 1, 0))
  exp_lm = as.data.frame(summary(glm(HasCAT ~ OutlierValue, data = exp_data, family = binomial))$coefficients)[2,]
  colnames(exp_lm) = c('Beta', 'SE', 'Zval', 'Pval')
  new_coefs=exp_lm %>% mutate(Category = 'eOutliers', Variant = vcat)
  all_coefs = rbind(all_coefs, new_coefs)
}


#pcols = exp_data %>% filter(!(medz_bin %in% c(NA,"0.000~0.002"))) %>%
#  group_by(variant_cat) %>% sample_n(1) %>%
 # ungroup() %>% mutate(CatCol = rgb(color_R, color_G, color_B)) %>%
#  select(variant_cat, CatCol)
#plot_cols = c(pcols$CatCol[1:11],'#3e6690','#33669a',pcols$CatCol[14:16])
#names(plot_cols) = unique(pcols$variant_cat)

save(risks,  file=out_rdata_relative)
save(all_coefs, file=out_rdata_continuous)

