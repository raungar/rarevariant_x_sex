library("ggplot2")
library("dplyr")
library("mltools") #ecdf
library("data.table")
#x_subtypes_types_numrv_bins_linc_prot.txt.gz
infile_x_numrv<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/x_CADDtypesGQ5BlacklistRemovedALL_numrv_bins_linc_prot.txt.gz"
infile_x_sigdif<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/x_CADDtypesGQ5BlacklistRemovedALL_sigdif_bins_linc_prot.txt.gz"
infile_aut_numrv<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/chr7_CADDtypesSeenTwice_numrv_bins_linc_prot.txt.gz"
infile_aut_sigdif<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/chr7_CADDtypesSeenTwice_sigdif_bins_linc_prot.txt.gz"
infile_x_subtypes_numrv<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/x_subtypes_CADDtypesGQ5BlacklistRemovedALL_numrv_bins_linc_prot.txt.gz"
infile_x_subtypes_sigdif<-"/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8redo/genomic_only/x_subtypes_CADDtypesGQ5BlacklistRemovedALL_sigdif_bins_linc_prot.txt.gz"

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
to_plot<-x_numrv_snps
vartype<-"SNPs"
chrtype<-"chrX"
# aes(x=as.character(range_max), y=this_median,
#     lower = IQR1, upper = IQR3, 
#     ymin=this_min, ymax=this_max))+
ggplot(to_plot,aes(x=paste0(range_min, "-",range_max), y=this_median,fill=sex))+
  geom_bar(stat="identity", color="black",  position=position_dodge()) +
  geom_errorbar(aes(ymin=IQR1,ymax=IQR3),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c("#F0CE22","#7CCBEA")) + # ylim(c(0,0.3))+
  ggtitle(paste0("Number of RVs per 10kb on the ",chrtype,  ": ", vartype))+xlab("MAF Bin")+ylab("RVs per 10kb")
ggplot(to_plot,aes(x=paste0(range_min, "-",range_max), y=this_sd,color=sex))+
  geom_point(size=5,shape=18)+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+ylim(c(0,0.16))+
  ggtitle(paste0("SD on the ",chrtype,  ": ", vartype))+xlab("RVs per 10kb")+ylab("SD")


###SUBTYPES
x_numrv_subtypes_snps_wzeros<-rbind(x_numrv_subtypes_snps,
                                    c("SNPs",0.001,"chrX","NONPAR",0,0,0,0,0,0,0,0,0,"male",0,0.001),
                                    c("SNPs",0.001,"chrX","XCR1",0,0,0,0,0,0,0,0,0,"male",0,0.001),
                                    c("SNPs",0.001,"chrX","XCR2",0,0,0,0,0,0,0,0,0,"male",0,0.001),
                                    c("SNPs",0.001,"chrX","XTR",0,0,0,0,0,0,0,0,0,"male",0,0.001))

to_plot<-x_numrv_subtypes_snps %>% dplyr::filter(par=="PAR1" | par=="PAR2" | par=="NONPAR" ) %>%
  dplyr::filter(bin_max==0.001 | bin_max==0.25)
vartype<-"SNPs"
chrtype<-"X"
ggplot(to_plot,aes(x=par, y=as.numeric(this_median),
                   fill=sex,group=interaction(sex,par)))+
  geom_bar(stat="identity", color="black",  position=position_dodge()) +
  geom_errorbar(aes(ymin=as.numeric(IQR1),ymax=as.numeric(IQR3)),width=.2,position=position_dodge(.9))+
  # geom_boxplot(aes(x=par, y=this_median, group=interaction(sex,par),
  #                  #lower = IQR1, upper = IQR3, 
  #                  ymin=this_min, ymax=this_max),position="dodge")+
  facet_wrap(~paste0(bin_min, "-",bin_max),scales = "free_y")+
  #facet_wrap(~bin_max,scales = "free_y")+
  
  scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  ggtitle(paste0(" Median Number of RVs per 10kb on the ",chrtype,  ": ", vartype))+
  xlab("Subregion")+ylab("Median Number of RVs per 10kb")

ggplot(to_plot,aes(x=par, y=as.numeric(this_sd),color=sex, group=interaction(sex,par)))+
  geom_point(size=5,shape=18)+
  facet_wrap(~paste0(bin_min, "-",bin_max))+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  ggtitle(paste0("SD on the ",chrtype,  ": ", vartype))+xlab("RVs per 10kb")+ylab("SD")






# chrs<-c(1:22) #x
chrs<-c("7","x")
mydir<-c("/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/genomic_only")
suffix<-c("_all_numrv_bins_linc_prot.txt.gz")
all_chrs<-data.frame() #data.frame(matrix(nrow = 0,ncol=15))
for (this_chr in chrs){
  print(this_chr)
  this_file=paste0(mydir,"/",this_chr,suffix)
  
  varname<-paste0("bin_chr",this_chr)
  this_binned_chr=fread(this_file) %>% dplyr::filter(chr==paste0("chr",toupper(this_chr)))
  # load(file)
  assign(varname,this_binned_chr)
  all_chrs<-rbind(all_chrs,this_binned_chr,fill=T)
 # break
}
all_chrs$chr<-factor(all_chrs$chr, levels=c(paste0("chr",1:22))) #,"chrX"))
all_chrs<-all_chrs %>% dplyr::mutate(chr = fct_relevel(chr, 
                         c(paste0("chr",1:22)))) #,"chrX")))
all_chrs_m<-dplyr::filter(all_chrs,sex=="male")
all_chrs_f<-dplyr::filter(all_chrs,sex=="female")

ggplot(all_chrs_f%>%dplyr::filter(range_max !=0),aes(x=(chr), y=as.numeric(this_median),
                   fill=sex,group=interaction(chr)))+
  geom_bar(stat="identity", color="black",  position=position_dodge()) +
  geom_errorbar(aes(ymin=as.numeric(IQR1),ymax=as.numeric(IQR3)),width=.2,position=position_dodge(.9))+
  xlab("Subregion")+ylab("Median Number of RVs per 10kb")+ 

  # geom_boxplot(aes(x=par, y=this_median, group=interaction(sex,par),
  #                  #lower = IQR1, upper = IQR3, 
  #                  ymin=this_min, ymax=this_max),position="dodge")+
  facet_wrap(~paste0(range_min, "-",range_max),scales = "free_y")+
  #scale_fill_manual(values=c("#F0CE22","#7CCBEA"))+
  scale_fill_manual(values=c("#F0CE22"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste0("Median Number of RVs per 10kb across the Autosomes (in Females)"))+
  scale_x_discrete(c(paste0("chr",1:22))) #+
 #xlab("Subregion")+ylab("Median Number of RVs per 10kb")

ggplot(all_chrs,aes(x=(chr), y=as.numeric(this_sd),color=sex,alpha=0.7,
                    group=interaction(chr)))+
  facet_wrap(~maf,scales = "free_y")+
  geom_point(size=5,shape=18)+
  scale_color_manual(values=c("#F0CE22","#7CCBEA"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste0("SD Across Chromosomes by Sex "))+
  scale_x_discrete(c(paste0("chr",1:22)))+
  xlab("Subregion")+ylab("Median Number of RVs per 10kb")
