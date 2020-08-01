library("ggplot2")
library("dplyr")
library("mltools") #ecdf
library("data.table")

infile_x_numrv<-"Output/analysis_v8/genomic_only/x_all_numrv.txt.gz"
infile_x_sigdif<-"Output/analysis_v8/genomic_only/x_all_sigdif.txt.gz"
infile_aut_numrv<-"Output/analysis_v8/genomic_only/aut_all_numrv.txt.gz"
infile_aut_sigdif<-"Output/analysis_v8/genomic_only/aut_all_sigdif.txt.gz"
infile_x_subtypes_numrv<-"Output/analysis_v8/genomic_only/x_subtypes_all_numrv.txt.gz"
infile_x_subtypes_sigdif<-"Output/analysis_v8/genomic_only/x_subtypes_all_sigdif.txt.gz"


x_numrv<-fread(infile_x_numrv,data.table=F, header=T)
x_sigdif<-fread(infile_x_sigdif,data.table=F)

aut_numrv<-fread(infile_aut_numrv,data.table=F)
aut_numrv$chr<-factor(aut_numrv$chr, levels=c(paste0("chr",1:22)))
aut_sigdif<-fread(infile_aut_sigdif,data.table=F)

x_subtypes_numrv<-fread(infile_x_subtypes_numrv,data.table=F)
x_subtypes_sigdif<-fread(infile_x_subtypes_sigdif,data.table=F)
#split by variant types
# cum_maf_summ_plot_m_and_f_snps<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="SNPs")
# cum_maf_summ_plot_m_and_f_indels<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="indels")
# cum_maf_summ_plot_m_and_f_sv<-cum_maf_summ_plot_m_and_f %>% dplyr::filter(vartype=="SV")

x_numrv_snps<-x_numrv %>% dplyr::filter(vartype=="SNPs")%>% dplyr::filter(sex=="male" | sex== "female")
x_numrv_indels<-x_numrv %>% dplyr::filter(vartype=="indels") %>% dplyr::filter(sex=="male" | sex== "female")
x_numrv_sv<-x_numrv %>% dplyr::filter(vartype=="SV") %>% dplyr::filter(sex=="male" | sex== "female")

x_numrv_subtypes_snps<-x_subtypes_numrv %>% dplyr::filter(vartype=="SNPs")%>% dplyr::filter(sex=="male" | sex== "female")
x_numrv_subtypes_indels<-x_subtypes_numrv %>% dplyr::filter(vartype=="indels") %>% dplyr::filter(sex=="male" | sex== "female")
x_numrv_subtypes_sv<-x_subtypes_numrv %>% dplyr::filter(vartype=="SV") %>% dplyr::filter(sex=="male" | sex== "female")
x_subtypes_sigdif_indels<-x_subtypes_sigdif %>% dplyr::filter(vartype=="indels")
x_subtypes_sigdif_snps<-x_subtypes_sigdif %>% dplyr::filter(vartype=="SNPs")
x_subtypes_sigdif_sv<-x_subtypes_sigdif %>% dplyr::filter(vartype=="SV")

aut_numrv_snps<-aut_numrv %>% dplyr::filter(vartype=="SNPs")%>% dplyr::filter(sex=="male" | sex== "female")
aut_numrv_indels<-aut_numrv %>% dplyr::filter(vartype=="indels") %>% dplyr::filter(sex=="male" | sex== "female")
aut_numrv_sv<-aut_numrv %>% dplyr::filter(vartype=="SV") %>% dplyr::filter(sex=="male" | sex== "female")

aut_comb_numrv<-as.data.table(aut_numrv[,c("vartype","maf","chr","this_mean","sex")])[,
                                                       Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                                       by=.(vartype,maf,sex),.SDcols=c("this_mean")]
aut_comb_numrv_reformatted<-cbind(as.data.frame(aut_comb_numrv[,c(1,2)]), "chr"="aut", as.data.frame(aut_comb_numrv[,c(4:9)]), "this_sd"="NA","upper"="NA","lower"="NA","sex"=aut_comb_numrv$sex)
colnames(aut_comb_numrv_reformatted)<-colnames(aut_numrv_snps)
aut_comb_numrv_snps<-aut_comb_numrv_reformatted %>% dplyr::filter(vartype=="SNPs")%>% dplyr::filter(sex=="male" | sex== "female") %>% dplyr::mutate(chr="aut")
aut_comb_numrv_indels<-aut_comb_numrv_reformatted %>% dplyr::filter(vartype=="indels") %>% dplyr::filter(sex=="male" | sex== "female") %>% dplyr::mutate(chr="aut")
aut_comb_numrv_sv<-aut_comb_numrv_reformatted %>% dplyr::filter(vartype=="SV") %>% dplyr::filter(sex=="male" | sex== "female") %>% dplyr::mutate(chr="aut")

###CHR!!
all_numrv_snps<-rbind(aut_numrv_snps,x_numrv_snps)
all_numrv_indels<-rbind(aut_numrv_indels,x_numrv_indels)
all_numrv_sv<-rbind(aut_numrv_sv,x_numrv_sv)


aut_sigdif_indels<-aut_sigdif %>% dplyr::filter(vartype=="indels")
aut_sigdif_snps<-aut_sigdif %>% dplyr::filter(vartype=="SNPs")
aut_sigdif_sv<-aut_sigdif %>% dplyr::filter(vartype=="SV")

#and plot!!!
# 95% CI
to_plot<-all_numrv_sv
vartype<-"SVs"
chrtype<-"X and autosomes"

ggplot(to_plot,aes(x=maf,y=this_mean,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  geom_ribbon(aes(ymin=lower,ymax=upper, fill=sex), alpha=0.2)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype, " [95% CI]"))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean))

# IQR
ggplot(to_plot,aes(x=maf,y=this_mean,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  #xlim(0,0.01)+ylim(0,0.2)+
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean))
# IQR zoom
ggplot(to_plot,aes(x=maf,y=this_mean,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  xlim(0,0.01)+ylim(0,0.0004)+
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean))

## do the chromosomes
ggplot(to_plot,aes(x=maf,y=this_mean,color=chr))+
  # scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  #scale_fill_manual(values=rep("gray",length(unique(all_numrv_snps$chr))*2))+
  #xlim(0,0.01)+ylim(0,0.018)+
  #geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=interaction(chr,sex), color=sex), alpha=0.2)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean, linetype=sex))+guides(fill=FALSE)
ggplot(to_plot,aes(x=maf,y=this_mean,color=chr))+
  #xlim(0,0.01)+ylim(0,0.24)+
  #xlim(0,0.01)+ylim(0,0.018)+
  xlim(0,0.01)+ylim(0,0.0005)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean, linetype=sex))+guides(fill=FALSE)



#plot p-valuex x vartypes
ggplot(x_sigdif,aes(x=maf,y=-log(padj),group=vartype,color=vartype))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
 # xlim(0,0.01)+
  scale_color_manual(values=c("#CFAAFE", "#A474CF", "#5E2A77"))+
  ggtitle("Significance of Difference in Males vs. Females of Number of Rare Variants on the X")+xlab("MAF")+
  geom_line()
#plot p-valuex x vartypes
ggplot(aut_sigdif,aes(x=maf,y=-log(padj),group=interaction(vartype,chr),color=chr))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
  # xlim(0,0.01)+
  #scale_color_manual(values=c("#CFAAFE", "#A474CF", "#5E2A77"))+
  ggtitle("Significance of Difference in Males vs. Females of Number of Rare Variants on the Autosomes")+xlab("MAF")+
  geom_line(aes(linetype=vartype))

#plot p-valuex x vartypes
ggplot(aut_sigdif_indels,aes(x=maf,y=-log(padj),group=chr,color=chr))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
  #xlim(0,0.01)+
  #scale_color_manual(values=c("#CFAAFE", "#A474CF", "#5E2A77"))+
  ggtitle("Significance of Difference in Males vs. Females of Number of Rare Variants")+xlab("MAF")+
  geom_line()





to_plot<-x_numrv_subtypes_sv
vartype<-"SVs"
chrtype<-"X"

ggplot(to_plot,aes(x=maf,y=this_mean,color=par,group=interaction(sex,par)))+
 # geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=par), alpha=0.2, colour = NA)+
  # scale_color_manual(values=c("#704714","#2848E0","#51B1B8", "#C0456E","#C46023","#C42A23","#349078"))+
  scale_color_manual(values=c("#704714","#2848E0", "#C0456E","#C46023","#C42A23","#349078"))+
  #xlim(0,0.01)+ylim(0,.00112)+
  #xlim(0,0.01)+ylim(0,0.5)+
  #xlim(0,0.01)+ylim(0,0.043)+
  ggtitle(paste0("Average Number of RVs on the ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+ylab("Mean Number of RVs")+
  geom_line(aes(y=this_mean,linetype=sex))

#plot p-valuex x vartypes
ggplot(to_plot,aes(x=maf,y=-log(padj),color=par_region,group=interaction(vartype,par_region)))+
  geom_hline(yintercept=-log(0.043), linetype="dashed", color="gray")+
  scale_color_manual(values=c("#704714","#2848E0", "#C0456E","#C46023","#C42A23","#349078"))+
  xlim(0,0.01)+ylim(0,25)+
  ggtitle(paste0("M vs. F ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+
  #ggtitle("Significance of Difference in Males vs. Females of Number of RVs (indels)")+
geom_line(aes(linetype=vartype))

#



########OLD
# 
# 
# infile<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/all_rvs_inds_types.txt.gz"
# euro_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_v8_euro_VCFids.txt"
# sex_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# wrong_self_reported_ancestry<-c("GTEX-11TT1", "GTEX-12ZZX", "GTEX-131XF", "GTEX-147F3","GTEX-15SZO","GTEX-16XZY","GTEX-17EUY","GTEX-17HHE", 
#                                 "GTEX-18D9U", "GTEX-1C4CL", "GTEX-1IDJC", "GTEX-1RMOY", "GTEX-R53T", "GTEX-R55D", "GTEX-WHPG", "GTEX-XMD3","GTEX-YB5E")
# 
# #Read everything in
# exp_data = fread(infile,data.table=F)
# colnames(exp_data)<-c("chr","start","end","maf","gtex_sample","vartype")
# euro<-fread(euro_file)
# euro_vec_unchecked<-as.character(data.frame(euro)[,1])
# euro_vec<-euro_vec_unchecked[!(euro_vec_unchecked %in% wrong_self_reported_ancestry )]
# exp_data_euro<-exp_data[exp_data$gtex_sample %in% euro_vec,]
# 
# sex_df<-fread(sex_file,data.table=F)
# sex_hash<-sex_df$SEX
# names(sex_hash)<-sex_df$SUBJID
# #temp<-rbind(head(exp_data,5000),head(tail(head(exp_data,100000),5000)),tail(head(exp_data,35000),5000))
# #this uses data table (efficient for large dataset)
# #uses mltools empirical_cdf function to get a CDF that also contains num cumulative
# # such that it outputs the cumulative number of RV at a certain MAF
# #such that that min_y=0 and max_y=(total #RV)
# 
# #takes in cumulative MAF matrix and then across individuals for the same vartype/chr
# #calculatates the mean, median, and sd at a given MAF and returns this df
# #which has cols vartype,maf,chr,this_mean,this_median,this_sd
# get_summary<-function(this_cum_maf){
#   
#   cum_maf_sd<-as.data.table(this_cum_maf)[,
#                                           Reduce(c,lapply(.SD,stats::sd)), 
#                                           by=.(vartype,UpperBound,chr),.SDcols="N.cum"]
#   cum_maf_summary<-as.data.table(this_cum_maf)[,
#                                                Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
#                                                by=.(vartype,UpperBound,chr),.SDcols="N.cum"]
#   
#   cum_maf_df<-merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr"), all=T)
#   colnames(cum_maf_df)<-c("vartype","maf","chr","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
#   upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
#   lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
#   cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
#   return(as.data.frame(cum_maf_df_CI))
# }
# 
# #same as get_summary but includes PAR as consideration
# get_summary_par<-function(this_cum_maf, this_sdcol){
#   
#   cum_maf_sd<-as.data.table(this_cum_maf)[,
#                                           Reduce(c,lapply(.SD,stats::sd)), 
#                                           by=.(vartype,UpperBound,chr,par_region),.SDcols=this_sdcol]
#   cum_maf_summary<-as.data.table(this_cum_maf)[,
#                                                Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
#                                                by=.(vartype,UpperBound,chr,par_region),.SDcols=this_sdcol]
#   
#   cum_maf_df<-merge(x = cum_maf_summary, y = cum_maf_sd, by = c("vartype","UpperBound","chr","par_region"), all=T)
#   colnames(cum_maf_df)<-c("vartype","maf","chr","par","this_min","IQR1","this_median","this_mean","IQR3","this_max","this_sd")
#   upper=cum_maf_df$this_mean+2*cum_maf_df$this_sd
#   lower=cum_maf_df$this_mean-2*cum_maf_df$this_sd
#   cum_maf_df_CI<-cbind(cum_maf_df,upper,lower)
#   return(as.data.frame(cum_maf_df_CI))
# }
# #mann-whitney-wilcoxon test across MAF on m vs. F for a given vartype
# get_mww_test_m_f<-function(maf_m,maf_f,vartype,mult_test_type){
#   mw_test_p<-unlist(lapply(sort(unique(maf_m$UpperBound)),
#                            function(x)
#                              wilcox.test(as.numeric(as.data.frame(maf_m)[maf_m$UpperBound==x & maf_m$vartype==vartype,"N.cum"]),
#                                          as.numeric(as.data.frame(maf_f)[maf_f$UpperBound==x & maf_f$vartype==vartype ,"N.cum"]))$p.value)
#   )
#   mw.adj.p<-p.adjust(mw_test_p,method =mult_test_type)
#   mw_df<-data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p,"padj"=mw.adj.p, "vartype"=vartype)
#   return(mw_df)
#   
# }
# #again, par considered
# get_mww_test_m_f_par<-function(maf_m,maf_f,this_var,mult_test_type="BH"){
#   print(this_var)
#   mw_test_p_par1<-unlist(lapply(sort(unique(maf_m$UpperBound)),
#                                 function(x)
#                                   wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR1" & UpperBound==x) %>% pull(N.cum)),
#                                               as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="PAR1" & UpperBound==x) %>% pull(N.cum)))$p.value
#   ))
#   mw_test_p_nonpar<-unlist(lapply(sort(unique(maf_m$UpperBound)),
#                                   function(x)
#                                     wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="NONPAR" & UpperBound==x) %>% pull(N.cum)),
#                                                 as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="NONPAR" & UpperBound==x) %>% pull(N.cum)))$p.value
#   ))
#   if(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR2") %>% nrow() ==0){
#     mw_df<-rbind(
#       data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par1,"padj"=p.adjust(mw_test_p_par1,method =mult_test_type), "vartype"=this_var, "par_region"="PAR1"),
#       data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_nonpar,"padj"=p.adjust(mw_test_p_nonpar,method =mult_test_type), "vartype"=this_var, "par_region"="NONPAR")
#     )
#   }else{
#     mw_test_p_par2<-unlist(lapply(sort(unique(maf_m$UpperBound)),
#                                   function(x)
#                                     wilcox.test(as.numeric(maf_m %>% dplyr::filter(vartype==this_var & par_region=="PAR2" & UpperBound==x) %>% pull(N.cum)),
#                                                 as.numeric(maf_f %>% dplyr::filter(vartype==this_var & par_region=="PAR2" & UpperBound==x) %>% pull(N.cum)))$p.value
#     ))
#     mw_df<-rbind(
#       data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par1,"padj"=p.adjust(mw_test_p_par1,method =mult_test_type), "vartype"=this_var, "par_region"="PAR1"),
#       data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_nonpar,"padj"=p.adjust(mw_test_p_nonpar,method =mult_test_type), "vartype"=this_var, "par_region"="NONPAR"),
#       data.frame("maf"=sort(unique(maf_m$UpperBound)),"pval"=mw_test_p_par2,"padj"=p.adjust(mw_test_p_par2,method =mult_test_type), "vartype"=this_var, "par_region"="PAR2")
#     )
#   }
#   
#   return(mw_df)
#   
# }
# 
# #uses mltools package to calculate the cumulative number of RVs at each MAF for
# #each individual x vartype x chr x maf level, going from 0 to .5 in 0.001 intervals
# cum_maf<-as.data.table(exp_data_euro)[,
#                                       Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
#                                       by=.(gtex_sample,vartype,chr),.SDcols="maf"]
# 
# #subset to m/f/both
# this_cum_maf_both<-cum_maf
# this_cum_maf_m<-cum_maf[sex_hash[cum_maf$gtex_sample]==1,]
# this_cum_maf_f<-cum_maf[sex_hash[cum_maf$gtex_sample]==2,]
# this_cum_maf_f_half<-this_cum_maf_f
# this_cum_maf_f_half$N.cum<-this_cum_maf_f_half$N.cum/2
# 
# 
# #get mean, median, and sd at a given MAF summary for indiviudals
# #at same vartype/maf/chr
# cum_maf_summ_all<-get_summary(this_cum_maf_both)
# cum_maf_summ_m<-get_summary(this_cum_maf_m)
# cum_maf_summ_f<-get_summary(this_cum_maf_f)
# cum_maf_summ_f_half<-get_summary(this_cum_maf_f_half)
# #combine m+f for plotting

# 
# 
# 
# #PAR REGIONS
# par1_s=10001; par1_e=2781479; par2_s=155701383; par2_e=156030895
# xar_s=2731479;xar_e=58555579
# xcr1_s=62462543; xcr1_e=89140830;
# xtr_s=89140830;xtr_e=93428068;
# xcr2_s=93428068;xcr2_e=155701383
# len_par1=par1_e-par1_s;len_nonpar=par2_s-par1_e; len_par2=par2_e-par2_s
# len_xcr1=xcr1_e-xcr1_s; len_xcr2=xcr2_e-xcr2_s
# len_xtr=xtr_e-xtr_s; len_xar=xar_e-xar_s; 
# 
# par_hash<-c(len_par1, len_nonpar,len_par2,len_xcr1,len_xcr2,len_xtr,len_xar)
# names(par_hash)<-c("PAR1","NONPAR","PAR2", "XCR1","XCR2","XTR","XAR")
# 
# par_region=mapply(function(s,e){
#   if(s>=par1_s & e<=par1_e){
#     "PAR1"
#   }else if(s>=par2_s & e<=par2_e){
#     "PAR2"
#   } else if(){
#   }
#   else{
#     "NONPAR"
#   }
# },exp_data_euro$start,exp_data_euro$end)
# par_binary=mapply(function(s,e){
#   if(s>=par1_s & e<=par1_e){
#     "PAR"
#   }else if(s>=par2_s & e<=par2_e){
#     "PAR"
#   }else{
#     "NONPAR"
#   }
# },exp_data_euro$start,exp_data_euro$end)
# exp_data_euro_wpars<-cbind(exp_data_euro,par_region,par_binary)
# cum_maf_parreg<-as.data.table(exp_data_euro_wpars)[,
#                                                    Reduce(c,lapply(.SD, function(x) as.list(empirical_cdf(x,ubounds=seq(0, .25, by=0.001))))),
#                                                    by=.(gtex_sample,vartype,chr,par_region),.SDcols="maf"]
# 
# #adjust number per 1000 bp
# cum_maf_parreg$N.cum.adj<-cum_maf_parreg$N.cum/par_hash[cum_maf_parreg$par_region]*10000
# #subset to m/f/both
# this_cum_maf_parreg_m<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==1,]
# this_cum_maf_parreg_f<-cum_maf_parreg[sex_hash[cum_maf_parreg$gtex_sample]==2,]
# this_cum_maf_parreg_f_half<-this_cum_maf_parreg_f
# this_cum_maf_parreg_f_half$N.cum<-this_cum_maf_parreg_f$N.cum/2
# this_cum_maf_parreg_f_half$N.cum.adj<-this_cum_maf_parreg_f$N.cum.adj/2
# 
# 
# #get mean, median, and sd at a given MAF summary for indiviudals
# #at same vartype/maf/chr
# cum_maf_summ_parreg_m<-get_summary_par(this_cum_maf_parreg_m, "N.cum")
# cum_maf_summ_parreg_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum")
# cum_maf_summ_parreg_f_half_par<-get_summary_par(this_cum_maf_parreg_f_half,"N.cum")
# 
# 
# #get mean, median, and sd at a given MAF summary for indiviudals
# #at same vartype/maf/chr
# cum_maf_summ_parreg_adj_m<-get_summary_par(this_cum_maf_parreg_m,"N.cum.adj")
# cum_maf_summ_parreg_adj_f<-get_summary_par(this_cum_maf_parreg_f,"N.cum.adj")
# cum_maf_summ_parreg_f_half_adj_par<-get_summary_par(this_cum_maf_parreg_f_half,"N.cum.adj")
# 
# 
# #combine m+f for plotting
# cum_maf_summ_plot_m_and_f_parreg<-rbind(cbind(cum_maf_summ_parreg_m, "sex"="male"),
#                                         cbind(cum_maf_summ_parreg_f_half_par, "sex"="female"))
# cum_maf_summ_plot_m_and_f_parreg_adj<-rbind(cbind(cum_maf_summ_parreg_adj_m, "sex"="male"),
#                                             cbind(cum_maf_summ_parreg_f_half_adj_par, "sex"="female"))
