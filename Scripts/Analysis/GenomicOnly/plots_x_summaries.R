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
x_maf_0.01_file<-"Output/analysis_v8/genomic_only/x_maf_0.01.txt.gz"

x_numrv<-fread(infile_x_numrv,data.table=F, header=T)
x_sigdif<-fread(infile_x_sigdif,data.table=F)
aut_numrv<-fread(infile_aut_numrv,data.table=F)
aut_numrv$chr<-factor(aut_numrv$chr, levels=c(paste0("chr",1:22)))
aut_sigdif<-fread(infile_aut_sigdif,data.table=F)
x_subtypes_numrv<-fread(infile_x_subtypes_numrv,data.table=F)
x_subtypes_sigdif<-fread(infile_x_subtypes_sigdif,data.table=F)
x_maf_0.01<-fread(x_maf_0.01_file,data.table=F)

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

aut_comb_numrv<-as.data.table(aut_numrv[,c("vartype","maf","chr","this_median","sex")])[,
                                                       Reduce(c,lapply(.SD,function(x) as.list(summary(x)))), 
                                                       by=.(vartype,maf,sex),.SDcols=c("this_median")]
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
to_plot<-x_numrv_snps
vartype<-"SNPs"
chrtype<-"X"

ggplot(to_plot,aes(x=maf,y=this_median,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  geom_ribbon(aes(ymin=lower,ymax=upper, fill=sex), alpha=0.2)+
  ggtitle(paste0("Median Number of RVs/10kb on the",chrtype,  ": ", vartype, " [95% CI]"))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median))

# IQR
ggplot(to_plot,aes(x=maf,y=this_median,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  #xlim(0,0.01)+ylim(0,0.2)+
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  ggtitle(paste0("Median Number of RVs/10kb on the ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median))
# IQR zoom
ggplot(to_plot,aes(x=maf,y=this_median,color=sex, group=sex))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  xlim(0,0.01) + ylim(0,0.2)+
  #xlim(0,0.01) + ylim(0,0.012)+
 # xlim(0,0.01) + ylim(0,0.0004)+
  
  geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=sex), alpha=0.2, colour = NA)+
  ggtitle(paste0("Median Number of RVs/10kb on the ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median))

## do the chromosomes
ggplot(to_plot,aes(x=maf,y=this_median,color=chr))+
  # scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  #scale_fill_manual(values=rep("gray",length(unique(all_numrv_snps$chr))*2))+
  #xlim(0,0.01)+ylim(0,0.018)+
  #geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=interaction(chr,sex), color=sex), alpha=0.2)+
  ggtitle(paste0("Median Number of RVs/10kb on the ",chrtype,  ": ", vartype))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median, linetype=sex))+guides(fill=FALSE)
ggplot(to_plot,aes(x=maf,y=this_median,color=chr))+
  #xlim(0,0.01)+ylim(0,0.22)+
  #xlim(0,0.01)+ylim(0,0.018)+
  xlim(0,0.01)+ylim(0,0.0004)+
  ggtitle(paste0("Median Number of RVs/10kb on the ",chrtype,  ": ", vartype))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median, linetype=sex))+guides(fill=FALSE)



#plot p-valuex x vartypes
ggplot(x_sigdif,aes(x=maf,y=-log(padj),group=vartype,color=vartype))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
 # xlim(0,0.01)+
  scale_color_manual(values=c("#CFAAFE", "#A474CF", "#5E2A77"))+
  ggtitle("M vs F #RVs on the X")+xlab("MAF")+
  geom_line()
#plot p-valuex x vartypes
ggplot(aut_sigdif,aes(x=maf,y=-log(padj),group=interaction(vartype,chr),color=chr))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
   xlim(0,0.01)+
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





to_plot<-x_subtypes_sigdif_snps
vartype<-"SNPs"
chrtype<-"X"

ggplot(to_plot,aes(x=maf,y=this_median,color=par,group=interaction(sex,par)))+
 # geom_ribbon(aes(ymin=IQR1,ymax=IQR3, fill=par), alpha=0.2, colour = NA)+
  scale_color_manual(values=c("#704714","#2848E0","#51B1B8", "#C0456E","#C46023","#C42A23","#349078"))+
  #scale_color_manual(values=c("#704714","#2848E0", "#C0456E","#C46023","#C42A23","#349078"))+
  #xlim(0,0.01)+ylim(0,.00038)+
 xlim(0,0.01)+ylim(0,0.6)+
  #xlim(0,0.01)+ylim(0,0.038)+
  ggtitle(paste0("Median Number of RVs/10kb on the ",chrtype,  ": ", vartype))+xlab("MAF")+ylab("Median Number of RVs/10kb")+
  geom_line(aes(y=this_median,linetype=sex))

#plot p-valuex x vartypes
ggplot(to_plot,aes(x=maf,y=-log(padj),color=par_region,group=interaction(vartype,par_region)))+
  geom_hline(yintercept=-log(0.05), linetype="dashed", color="gray")+
  scale_color_manual(values=c("#704714","#2848E0","#51B1B8", "#C0456E","#C46023","#C42A23","#349078"))+
  #scale_color_manual(values=c("#704714","#2848E0", "#C0456E","#C46023","#C42A23","#349078"))+
  # xlim(0,0.01)+ #ylim(0,22)+
  ggtitle(paste0("M vs. F ",chrtype,  ": ", vartype,  " [IQR]"))+xlab("MAF")+
  #ggtitle("Significance of Difference in Males vs. Females of Number of RVs (indels)")+
geom_line(aes(linetype=vartype))

#MAF PLOTTING
x_maf_0.01_indels<-x_maf_0.01 %>% dplyr::filter(vartype=="indels")
x_maf_0.01_snps<-x_maf_0.01 %>% dplyr::filter(vartype=="SNPs")
x_maf_0.01_sv<-x_maf_0.01 %>% dplyr::filter(vartype=="SV")
vartype<-"SNPs"

to_plot<-x_maf_0.01_snps
set.seed(4791234)
to_plot_sample_males<-sample(which(to_plot[,"sex"]=="male"),length(which(to_plot[,"sex"]=="female")))
to_plot_females<-which(to_plot[,"sex"]=="female")
to_plot_red<-to_plot[c(to_plot_females,to_plot_sample_males ),]
ggplot(to_plot_red,aes(x=N.cum.adj,color=sex,fill=sex, group=sex))+
 # xlim(0,0.01)+ylim(0,0.24)+
  #xlim(0,0.01)+ylim(0,0.018)+
  #xlim(0,0.01)+ylim(0,0.0005)+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  geom_histogram(bins=50,alpha=0.5)+
  ggtitle(paste0("Number of RVs per 10kb on the ",chrtype,  ": ", vartype))+xlab("RVs per 10kb")+ylab("Freq")
  #geom_line(aes(y=N.cum.adj, linetype=sex))+guides(fill=FALSE)




