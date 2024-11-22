library("ggplot2")
library("dplyr")
library("epitools")
library("readr")
library("data.table")

#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8"
groups=c("m","f","both_half_regress") #, "both_half","both_half.sex","both_half.regress
chr_types=c("aut","x")
nphen=c(2,3,4,5)
z=c(2,2.5,3)
maf_min=0.01

risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(),chr=character(), sex=character(),zmin=numeric(), nphen_min=numeric(),num_outliers=numeric())


#read in all combination of CSV files, and assign to the appropriately named variables
for (this_sex in groups){
  for(this_chr in chr_types){
    for(this_nphen in nphen){
     for(this_z in z){
        
        file=paste0(mydir,"/outlier_noglobal_medz_varAnnot_zthresh",this_z,"_nphen",this_nphen,"_",
                    this_chr,"_",this_sex,"_linc_prot.txt.gz")
       varname<-paste0(this_chr,"_",this_sex,"_z",this_z,"_nphen",this_nphen)
       print(file)
        combined_genes_vars=fread(file,sep="\t")
        print("READ")
	    combined_genes_vars<-combined_genes_vars %>% dplyr::filter(vartype=="SNP" | is.na(vartype)) ###limit to only a SNP analysis
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
        
	    rm(varname)
    	print(risks)
      }
    }
  }
}

#filter to outliers and controls

risks$CI<-risks$Upper-risks$Lower


risks = risks %>% arrange(by=Risk) 
#write.csv(risks,file=paste0(mydir,"/all_rr.csv"))
#risks<-read_csv(paste0(mydir,"/all_rr.csv"))
risks_x<-dplyr::filter(risks,chr=="x")
risks_aut<-dplyr::filter(risks,chr=="aut")

library(RColorBrewer)
library(viridis)
to_plot<-risks_aut%>%as.data.frame() #%>%dplyr::filter(zmin==2)
to_plot[to_plot[,"sex"]=="both_half_regress","sex"]<-"both"
#, labeller =c("Z>2","Z>2.5","Z>3")
to_plot[to_plot[,"zmin"]=="2","zmin"]<-"Z>2"
to_plot[to_plot[,"zmin"]=="2.5","zmin"]<-"Z>2.5"
to_plot[to_plot[,"zmin"]=="3","zmin"]<-"Z>3"
# to_plot$CI=to_plot$Upper-to_plot$Lower

ggplot(to_plot,aes(x=nphen_min,y=num_outliers,group=sex,fill=sex))+
  scale_fill_manual(values=c("#B1EAA2","#FCF5A9","#97D6F2"))+
  geom_text(aes(label =round(num_outliers,2)), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(~zmin,scales = "free")+ylab("Number of Outliers")+xlab("Minimum number of tissues in which gene is outlier")+
  geom_bar(stat="identity",position = "dodge")

ggplot(to_plot,aes(x=nphen,y=sex,fill=CI))+
 geom_tile()+
  ggtitle("Z=2")+
  facet_grid(~Cat,scales = "free")+
 # scale_fill_gradient(low="#D1EDF5", high="#296D65")+
  scale_fill_gradient(low="white", high="orange")+
  geom_text(aes(label =round(CI,2))) +xlab("Number of Tissues")+ylab("Sex")
  #scale_fill_viridis(option = "magma",alpha = 0.8)+
    
  # scale_fill_brewer(palette = "BrBG")

########playing
both<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_both_zthresh3_nphen5_noglobal_medz_aut_both_half_regress.txt.gz")
both_expression<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_aut_both_half_regress.txt.gz")
female<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_both_zthresh3_nphen5_noglobal_medz_aut_f.txt.gz")
female_expression<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/gtex_2017-06-05_normalized_expression_subsetted_aut_f.txt.gz")
male<-fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8/outliers_both_zthresh3_nphen5_noglobal_medz_aut_m.txt.gz")
both_inds<-unique(both$Ind)
female<-female%>%dplyr::filter(Ind %in% both_inds)
male<-male%>%dplyr::filter(Ind %in% both_inds)
f_both<-merge(female, both,by=c("Ind","Gene"))
m_both<-merge(male, both,by=c("Ind","Gene"))
f_both_outlier<-f_both%>%dplyr::filter(Y.x=="outlier"|Y.y=="outlier")
f_both_outlier$diff<-abs(f_both_outlier$MedZ.x-f_both_outlier$MedZ.y)
m_both_outlier<-m_both%>%dplyr::filter(Y.x=="outlier"|Y.y=="outlier")
m_both_outlier$diff<-abs(m_both_outlier$MedZ.x-m_both_outlier$MedZ.y)


both_gene<-both_expression %>%dplyr::filter(Gene=="ENSG00000157326.18")
f_gene<-female_expression %>%dplyr::filter(Gene=="ENSG00000157326.18")

ENSG00000100591_tpm<-fread("/oak/stanford/groups/smontgom/raungar/Sex/tpm.tmp")
ENSG00000170606_both=fread("/oak/stanford/groups/smontgom/raungar/Sex/ENSG00000170606.both.tmp")
ENSG00000170606_m=fread("/oak/stanford/groups/smontgom/raungar/Sex/ENSG00000170606.m.tmp")
ENSG00000170606_f=fread("/oak/stanford/groups/smontgom/raungar/Sex/ENSG00000170606.f.tmp")
tiss_m<-ENSG00000170606_m%>%dplyr::filter(Tissue=="Artery_Aorta")
tiss_f<-ENSG00000170606_f%>%dplyr::filter(Tissue=="Artery_Aorta");tiss_f[,"GTEX-132AR"]
tiss_both<-ENSG00000170606_both%>%dplyr::filter(Tissue=="Artery_Aorta");tiss_both[,"GTEX-132AR"]

this_both_gene<-data.frame(both_gene[23,-c(1,2)])
colnames(this_both_gene)<-str_replace(colnames(this_both_gene),"\\.","-")
both_m<-this_both_gene[,names(this_both_gene)[!(names(this_both_gene) %in% colnames(f_gene))]]
both_f<-this_both_gene[,names(this_both_gene)[(names(this_both_gene) %in% colnames(f_gene))]]
# both_m<-this_both_gene[,!(colnames(this_both_gene)%in%colnames(f_gene))]
# both_f<-this_both_gene %>%dplyr::pull(colnames(this_both_gene) %in% colnames(f_gene))
both_separated<-melt(data.frame(m=as.numeric(both_m),f=c(as.numeric(both_f),rep(NA,7))))
m_hist=hist(as.numeric(both_m),breaks=seq(from =-4, to = 3, by = 0.25))
f_hist=hist(as.numeric(both_f),breaks=seq(from =-4, to = 3, by = 0.25))
combined_hist=data.frame(Zscore=paste0(m_hist$breaks[-length(m_hist$breaks)]," to ",m_hist$breaks[-1]),
                         female=f_hist$counts,
                         male=m_hist$counts)
combined_hist_melt<-melt(combined_hist)
colnames(combined_hist_melt)[c(2,3)]<-c("Sex","Frequency")
combined_hist_melt$Zscore<-factor(combined_hist_melt$Zscore,levels=combined_hist$Zscore)
ggplot(combined_hist_melt,aes(x=Zscore,y=Frequency,group=Sex,fill=Sex))+
  scale_fill_manual(values=c("#F0CE22","#50a4d4"))+
  geom_bar(stat="identity",position = "dodge")+
geom_vline(xintercept=as.numeric(both_gene[23,"GTEX-UJHI"]),col="red")


ggplot(both_separated, aes(x=value, fill=variable)) +
  scale_fill_manual(values=c("#50a4d4", "#F0CE22")) +
  geom_density( alpha=0.8, position = 'identity')+
  geom_vline(xintercept=as.numeric(both_gene[23,"GTEX-UJHI"]),col="red")


ggplot()+aes(na.omit(as.numeric((f_gene[23,-c(1,2)])))) +
   geom_density( fill="#F0CE22", color="#F0CE22", alpha=0.9,bins=40) +  xlim(c(-4,4))+
  # geom_density( fill="#6cb86e", color="#6cb86e", alpha=0.9,bins=40) +  xlim(c(-4,4))+
  ylab("Count")+ggtitle("DHRS4 Muscle Skeletal (Female Expression)")+xlab("Z-score")+
geom_vline(xintercept=as.numeric(f_gene[23,"GTEX-UJHI"]),col="red")
#"#F0CE22","#7CCBEA"  