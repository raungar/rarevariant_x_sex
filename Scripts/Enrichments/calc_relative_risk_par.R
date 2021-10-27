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
                make_option(c("--gtf_code_file"), type = 'character', default = NULL, help = "gtf_code_file preprocessing"),
                make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--par_file"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--max_maf"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--min_maf"), type = 'character', default = NULL, help = "path of x inactivation file"),
                make_option(c("--nphen"), type = 'numeric', default = NULL, help = "min tissues"),
                make_option(c("--out_rdata_relative"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
                make_option(c("--cadd_min"), type = 'numeric', default = NULL, help = "min cadd score")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
par_file <- as.character(opt$par_file)
out_rdata_relative <- as.character(opt$out_rdata_relative)
my_zscore <- as.numeric(opt$zscore)
cadd_min <- as.numeric(opt$cadd_min)
max_maf <- as.numeric(opt$max_maf)
min_maf <- as.numeric(opt$min_maf)
sex <- as.character(opt$sex)
nphen <- as.numeric(opt$nphen)
gtf_code_file<-as.character(opt$gtf_code_file)
# # #zscore<-3
#infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/x_outlier_noglobal_medz_varAnnot_zthresh2.5_nphen3_m_CADDtypesGQ5SeenTwice_linc_prot.txt.gz"
# out_rdata_relative<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_z3_x_f.xci.RData"
#  par_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/par_table.txt"
#  gtf_code_file="/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/autosomal_proteincoding_lncrna.gtf"
# # zscore<-2

if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

   
### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","chr","start","end","vartype","sex",
                      "gtex_maf","gnomad_maf","use_maf","genetype","numrv","cadd_raw","cadd_phred")
#exp_data = exp_data %>% select(indiv_id,gene_id,MedZ.x,Y,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
# exp_data = exp_data %>% dplyr::select(indiv_id,gene_id,MedZ.x,Y,af,isOutlier,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
# exp_data_final = dplyr::filter(exp_data,sv_v7==1)
#new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$variant_cat[x] == 'splice', exp_data$tier2[x], exp_data$variant_cat[x]))
#exp_data_final$variant_cat = new_cats
exp_data$gene_id_red<-sapply(strsplit(exp_data$ensg,"\\."), "[[",1)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
   this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
   this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
   #no variants found w/in 10kb of gene, so not rare
   if(is.na(this_maf)){"common"}
   else if (this_maf>=min_maf & this_maf<max_maf){
      if(this_cadd_phred>=cadd_min){"rare"}
      else if(this_cadd_phred<cadd_min){"common"}
      else (stop("ERROR: VARIANT NOT MAKING SENSE"))
   }
   else if (this_maf>=max_maf){"common"}
   else (stop("ERROR: VARIANT NOT MAKING SENSE"))
})
exp_data$has_variant<-has_variant

gtf_code<-fread(gtf_code_file,header=F)
genetype_dic=gtf_code$V2
names(genetype_dic)<-gtf_code$V1
exp_data$variant_cat<-genetype_dic[exp_data$ensg] 

par_df<-fread(par_file,data.table=F,header = F)
colnames(par_df)<-c("ensg","subregion")
#par_df$GeneID_red<-sapply(strsplit(par_df$`Gene ID`,"\\."), "[[",1)
# exp_data<-merge(exp_data,par_df[c("PAR_BINARY", "GeneID_red")], by.x="gene_id_red",by.y="GeneID_red")
#colnames(exp_data)[colnames(exp_data)=="Combined XCI status"]<-"XCI_STATUS"

print("filter")
### get relative risk per category

exp_outliers = exp_data %>% dplyr::filter(abs(MedZ) >= my_zscore)  %>% dplyr::filter(as.numeric(Df)>=nphen)
exp_controls = dplyr::filter(exp_data, (abs(MedZ) < my_zscore) | (abs(MedZ) >= my_zscore & as.numeric(Df)<nphen))

print("relative risk")
### Relative risk
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Subregion = character(),
                   exp_nn=numeric(),exp_ny=numeric(),exp_yn=numeric(),exp_yy=numeric(),num_outliers=numeric(),
                   sex=character(),z=numeric(),nphen=numeric(),cadd_min=numeric(),Type=character())
par_status = unique(par_df$subregion)
for (par in par_status) {
  this_par_subregion<-par_df%>%dplyr::filter(subregion==par)
   print(par)

   exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% this_par_subregion$ensg) %>% dplyr::filter( has_variant != "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% this_par_subregion$ensg) %>% dplyr::filter(has_variant == "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_yn = nrow(exp_outliers %>% dplyr::filter(ensg %in% this_par_subregion$ensg) %>% dplyr::filter(has_variant != "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    exp_yy = nrow(exp_outliers %>% dplyr::filter(ensg %in% this_par_subregion$ensg) %>% dplyr::filter(has_variant == "rare")) #%>% dplyr::filter(cadd_phred>cadd_min))
    number_of_outliers=nrow(exp_outliers%>% dplyr::filter(ensg %in% this_par_subregion$ensg))
   # exp_nn = nrow(dplyr::filter(exp_controls, PAR_BINARY != par))
   # exp_ny = nrow(dplyr::filter(exp_controls, PAR_BINARY == par))
   # exp_yn = nrow(dplyr::filter(exp_outliers, PAR_BINARY != par))
   # exp_yy = nrow(dplyr::filter(exp_outliers, PAR_BINARY == par))
   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
   err = epitab(exptable, method = 'riskratio')
   print(err)
   print(paste(sex,my_zscore,nphen,cadd_min,par,sep=","))
   new_df<-data.frame(Risk = err$tab[2,5],
                    Lower = err$tab[2,6],
                    Upper = err$tab[2,7],
                    Pval = err$tab[2,8],
                    Subregion = par,
                    exp_nn=exp_nn,
                    exp_ny=exp_ny,
                    exp_yn=exp_yn,
                    exp_yy=exp_yy,
                    num_outliers=number_of_outliers,
                    sex=sex,z=my_zscore,nphen=nphen,cadd=cadd_min,
                    Type = 'Subregion')
   print(new_df)
   risks = rbind(risks, new_df)  
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
risks$Subregion = factor(risks$Subregion, levels=unique(risks$Subregion))
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


write.csv(risks,  file=out_rdata_relative,quote=F,sep="\t")
print(paste0("COMPLETE: ",out_rdata_relative))
# save(all_coefs, file=out_rdata_continuous)

