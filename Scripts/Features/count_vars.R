library(data.table)

indir="/oak/stanford/groups/smontgom/raungar/Sex/"
#x_vars=fread("/oak/stanford/group/smontgom/raungar/Sex/Output/features_v8/Combined/x_all_rvs_inds_CADDtypesSeenTwice_linc_prot_window5000.txt.gz")
chrs=c("x",paste0("chr",1:22))
sex_spec_vars=data.frame()
all_rare_vars=data.frame()

for(this_chr in chrs){
  if(this_chr=="x"){
    vars=fread(paste0(indir,"Output/features_v8/Combined/",this_chr,"_all_rvs_inds_CADDtypesGQ5BlacklistRemovedALL_linc_prot_window5000.txt.gz"))
    
  }else{
    vars=fread(paste0(indir,"Output/features_v8/Combined/",this_chr,"_all_rvs_inds_CADDtypesALL_linc_prot_window5000.txt.gz"))
    
  }
  colnames(vars)<-c("chr","pos1","pos2","maf_gtex","maf_gnomad","maf_use","id","SNPs","ensg","prot_lnc","sex","loc","cadd","cadd_phred","x","gene_start","gene_end","vartype","varname")
  rare_vars=vars%>%dplyr::filter(maf_use<0.01)
  rare_vars=rare_vars%>%mutate(variant=paste0(chr,":",pos1,"-",pos2),chr=this_chr)
  all_rare_vars=rbind(all_rare_vars,rare_vars)
  sex_sp_variants=rare_vars%>%group_by(variant,sex)%>%summarise(n=n())%>%group_by(variant)%>%mutate(n=n())%>%dplyr::filter(n==1)
  print(nrow(sex_sp_variants))
  this_row=c(this_chr,table(sex_sp_variants$sex),sex_spec=nrow(sex_sp_variants),num_rare=nrow(rare_vars))
  sex_spec_vars=rbind(sex_spec_vars,this_row)
}
##RANGE??

colnames(sex_spec_vars)<-c("chr","female_only","male_only","num_only","num_rare")
sex_spec_vars_melt=melt(sex_spec_vars,id.vars="chr")
sex_spec_vars_melt$chr=factor(sex_spec_vars_melt$chr,levels=c(paste0("chr",1:22),"x"))
sex_spec_vars_melt$value=as.numeric(sex_spec_vars_melt$value)
ggplot(sex_spec_vars_melt%>%dplyr::filter(variable!="num_only"&variable!="num_rare"),
       aes(x=chr,y=value,fill=variable))+
  scale_fill_manual(values=c("#dbab3b","#5d8596"))+
  ylab("proportion of rare variants only in one sex")+xlab("chromosome")+
  geom_bar(stat="identity",position="fill")+
  theme_bw()
#fwrite(sex_spec_vars,file=paste0(indir,"Output/features_v8/sex_spec_vars.txt.gz"))
#fwrite(all_rare_vars,file=paste0(indir,"Output/features_v8/all_rare_vars.txt.gz"))

sex_spec_vars=fread(file=paste0(indir,"Output/features_v8/sex_spec_vars.txt.gz"))

chr7=fread(paste0(indir,"Output/features_v8/Combined/chr7_all_rvs_inds_CADDtypesALL_linc_prot_window5000.txt.gz"))
colnames(chr7)<-c("chr","pos1","pos2","maf_gtex","maf_gnomad","maf_use","id","SNPs","ensg","prot_lnc","sex","loc","cadd","cadd_phred","x","gene_start","gene_end","vartype","varname")
rare_chr7=chr7%>%dplyr::filter(maf_use<0.01)
rare_chr7=rare_chr7%>%mutate(variant=paste0(chr,":",pos1,"-",pos2))
rare_chr7=rare_chr7%>%mutate(var2=as.numeric(paste0(pos1,pos2)),
                             sex2=ifelse(sex=="female",1,0))%>%
  mutate(var2sex2=as.numeric(paste0(pos1,sex2)))
num_rv_x=106266

bootstrap_chr7=function(rare_chr7,num_rv_x=106266,num_boot=10000){
  chr7_boot=data.frame()
  for(i in c(1:num_boot)){
    rare_chr7_bootstrap_rows=sample(1:nrow(rare_chr7),size=num_rv_x,replace = T) #0.02
    #rare_chr7_bootstrap=rare_chr7[rare_chr7_bootstrap_rows,c("var2","sex2","var2sex2")] #.15
    rare_chr7_bootstrap=rare_chr7[rare_chr7_bootstrap_rows,"var2sex2"] #.15
    #this is stupid, but for speed
    not_duplicated=rare_chr7_bootstrap[!(duplicated(rare_chr7_bootstrap)|duplicated(rare_chr7_bootstrap, fromLast=TRUE))] 
    not_duplicated_modulo=(not_duplicated+1)%%2 #if 0, that is female; otherwise male
    # rare_chr7_bootstrap_table=table(rare_chr7_bootstrap) #.5
    # rare_chr7_bootstrap_table_onesex=rare_chr7_bootstrap_table[rare_chr7_bootstrap_table[,"female"]==0|rare_chr7_bootstrap_table[,"male"]==0,]
    chr7_bootsrap_nfemale=length(not_duplicated_modulo[not_duplicated_modulo==0])
    chr7_bootsrap_nmale=length(not_duplicated_modulo[not_duplicated_modulo==1])
    chr7_boot=rbind(chr7_boot,c(chr7_bootsrap_nfemale,chr7_bootsrap_nmale))
  }
  return(chr7_boot)
}
my_boot_sex_rvs=bootstrap_chr7(rare_chr7,num_rv_x=106266,num_boot=10000)
#fwrite(my_boot_sex_rvs,file=paste0(indir,"Output/features_v8/chr7_bootstrap.txt.gz"))
my_boot_sex_rvs_mean=mean(my_boot_sex_rvs[,1])
my_boot_sex_rvs_sd=sd(my_boot_sex_rvs[,1])
(sex_spec_vars[1,2]-my_boot_sex_rvs_mean)/my_boot_sex_rvs_sd





chrx_rv_file="/oak/stanford/groups/smontgom/raungar/Sex"
chrx_rv_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/chrXResult_Merged.RDS"
aut_rv_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/AutosomeResult_Merged.RDS"
model_f_aut_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Female_AllAutosomes_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_m_aut_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Male_AllAutosomes_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_f_x_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Female_X_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_m_x_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Male_X_Eval_AutosomeAnnotations_SameModel.txt.gz"

aut_rv=readRDS(aut_rv_file)%>%mutate(gene=gsub("\\+.*","",geneRV))
chrx_rv=readRDS(chrx_rv_file)%>%mutate(gene=gsub("\\+.*","",geneRV))%>%mutate(within_female_diff=abs(FemaleMean-FemaleMax),within_male_diff=abs(MaleMean-MaleMax))
chr7_rv=(aut_rv%>%dplyr::filter(grepl("chr7",geneRV)))%>%mutate(gene=gsub("\\+.*","",geneRV))

functional_threshold=0
aut_rv_sexbiased=aut_rv%>%dplyr::filter(abs(diff)>.2& (abs(MaleMean)>functional_threshold|abs(FemaleMean)>functional_threshold))%>%mutate(mean_to_use=ifelse(abs(FemaleMean)>abs(MaleMean),FemaleMean,MaleMean),other_mean=ifelse(abs(FemaleMean)>abs(MaleMean),MaleMean,FemaleMean))
aut_not_sexbiased=aut_rv%>%dplyr::filter(abs(diff)<0.1 & (abs(MaleMean)>functional_threshold|abs(FemaleMean)>functional_threshold)) %>%mutate(mean_to_use=ifelse(abs(FemaleMean)>abs(MaleMean),FemaleMean,MaleMean),other_mean=ifelse(abs(FemaleMean)>abs(MaleMean),MaleMean,FemaleMean))#CHANGE 0.1, MEAN>0.2
aut_rv_male_over=aut_rv_sexbiased%>%dplyr::filter((MaleMean>FemaleMean)&sign(MaleMean)==1&MaleMean>functional_threshold)
aut_rv_male_under=aut_rv_sexbiased%>%dplyr::filter((MaleMean<FemaleMean)&sign(MaleMean)==-1&MaleMean<(-functional_threshold))
aut_rv_female_over=aut_rv_sexbiased%>%dplyr::filter((FemaleMean)>MaleMean&sign(FemaleMean)==1&FemaleMean>functional_threshold)
aut_rv_female_under=aut_rv_sexbiased%>%dplyr::filter((FemaleMean<MaleMean)&sign(FemaleMean)==-1&FemaleMean<(-functional_threshold))

x_rv_sexbiased=chrx_rv%>%dplyr::filter(abs(diff)>.2& (abs(MaleMean)>functional_threshold|abs(FemaleMean)>functional_threshold))%>%mutate(mean_to_use=ifelse(abs(FemaleMean)>abs(MaleMean),FemaleMean,MaleMean),other_mean=ifelse(abs(FemaleMean)>abs(MaleMean),MaleMean,FemaleMean))
x_not_sexbiased=chrx_rv%>%dplyr::filter(abs(diff)<0.1 & (abs(MaleMean)>functional_threshold|abs(FemaleMean)>functional_threshold))%>%mutate(mean_to_use=ifelse(abs(FemaleMean)>abs(MaleMean),FemaleMean,MaleMean),other_mean=ifelse(abs(FemaleMean)>abs(MaleMean),MaleMean,FemaleMean)) #CHANGE 0.1, MEAN>0.2
x_rv_male_over=x_rv_sexbiased%>%dplyr::filter((MaleMean>FemaleMean)&sign(MaleMean)==1&MaleMean>functional_threshold)
x_rv_male_under=x_rv_sexbiased%>%dplyr::filter((MaleMean<FemaleMean)&sign(MaleMean)==-1&MaleMean<(-functional_threshold))
x_rv_female_over=x_rv_sexbiased%>%dplyr::filter((FemaleMean)>MaleMean&sign(FemaleMean)==1&FemaleMean>functional_threshold)
x_rv_female_under=x_rv_sexbiased%>%dplyr::filter((FemaleMean<MaleMean)&sign(FemaleMean)==-1&FemaleMean<(-functional_threshold))

#bias_vs_nobias_posterior
posterior_bias_unbias=data.frame(rbind(cbind(cat="sex-biased",chr="x",posterior=x_rv_sexbiased$mean_to_use,other_posterior=x_rv_sexbiased$other_mean),
                                       cbind(cat="not-biased",chr="x",posterior=x_not_sexbiased$mean_to_use,other_posterior=x_not_sexbiased$other_mean),
                                       cbind(cat="sex-biased",chr="aut",posterior=aut_rv_sexbiased$mean_to_use,other_posterior=aut_rv_sexbiased$other_mean),
                                       cbind(cat="not-biased",chr="aut",posterior=aut_not_sexbiased$mean_to_use,other_posterior=aut_not_sexbiased$other_mean)
))
posterior_bias_unbias$posterior=as.numeric(posterior_bias_unbias$posterior)
posterior_bias_unbias$other_posterior=as.numeric(posterior_bias_unbias$other_posterior)

ggplot(posterior_bias_unbias,aes(x=chr,fill=cat))+
  geom_bar(stat="count",position="fill")+
  ylab("prop")+
  ggtitle(paste0("functional thresh ",functional_threshold))+
  theme_bw()

ggplot(posterior_bias_unbias,
       aes(y=interaction(cat,chr),x=posterior, height = ..density..))+
  geom_vline(xintercept = 0.2)+
  geom_vline(xintercept = -0.2)+
  geom_density_ridges(stat = "density", trim = TRUE)+
  #geom_density_ridges()+
  ggtitle("other posterior")+
  theme_bw()

ggplot(posterior_bias_unbias,
       aes(fill=interaction(cat,chr),alpha=0.5,x=other_posterior, height = ..density..))+
  geom_vline(xintercept = 0.2)+
  geom_vline(xintercept = -0.2)+
  # geom_density(stat = "density", trim = TRUE)+
  geom_histogram()+
  facet_wrap(~chr,scales="free_y")+
  ggtitle("<0.1")+
  scale_y_continuous(trans="log10")+
  #geom_density_ridges()+
  theme_bw()

ggplot(posterior_bias_unbias,aes(x=posterior,y=other_posterior,color=cat))+
  facet_wrap(~chr)+
  geom_point()+
  theme_bw()



reformat_model_df=function(filename,my_over,my_under,my_not_sexbiased,sex){
  mymodel=fread(filename)
  mymodel=mymodel%>%rename(beta=p_value)
  mymodel$variant=gsub("(.*)chr","chr",(mymodel$SubjectID)) #sep variant into its own col
  mymodel=mymodel%>%mutate(geneRV=paste0(GeneName,"+",variant))%>%as.data.table()
  
  #filtering, using datatable framework for speed
  setkey(mymodel,geneRV)
  mymodel_over=mymodel[.(my_over$geneRV), nomatch = 0L]%>%mutate(category=paste0(sex,"-biased (over-expression)"))
  mymodel_under=mymodel[.(my_under$geneRV), nomatch = 0L]%>%mutate(category=paste0(sex,"-biased (under-expression)"))
  mymodel_neutral=mymodel[.(my_not_sexbiased$geneRV), nomatch = 0L]%>%mutate(category=paste0("non-biased"))
  
  my_sexbiased_genes=return(rbind(mymodel_over,mymodel_under,mymodel_neutral))
  
}
functional_threshold=0.2
model_f_x=fread(model_f_x_file)
model_f_x=model_f_x%>%mutate(var=gsub()) ### make rare variant only
model_f_x_filt=model_f_x%>%dplyr::filter()
###per rare variant with beta > 0.2, get range

# model_f_aut=reformat_model_df(model_f_aut_file,aut_rv_female_over,aut_rv_female_under,aut_not_sexbiased,"female")
# model_m_aut=reformat_model_df(model_m_aut_file,aut_rv_male_over,aut_rv_male_under,aut_not_sexbiased,"male")
# model_all_aut=rbind(model_f_aut,model_m_aut)
# model_f_x=reformat_model_df(model_f_x_file,x_rv_female_over,x_rv_female_under,x_not_sexbiased,"female")
# model_m_x=reformat_model_df(model_m_x_file,x_rv_male_over,x_rv_male_under,x_not_sexbiased,"male")
# model_all_x=rbind(model_f_x,model_m_x)

#chrx_melt=melt(chrx_rv,id=c("geneRV","MaleMax","MaleMean","FemaleMax","FemaleMean","RV","gene"))
# ggplot(chrx_melt,aes(x=abs(value),color=variable))+
#   geom_density(stat = "density", trim = TRUE)+
#   scale_x_continuous(trans="log10")+
#   theme_bw()



#### XCI
xci=fread("/oak/stanford/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact.tsv")
colnames(xci)[7]="combined_xci_status"
xci$ensg=gsub("\\..*","",xci$`Gene ID`)
par=fread("/oak/stanford/groups/smontgom/raungar/Sex/Files/Tukiainen_xinact_par.tsv")
par$ensg=gsub("\\..*","",par$`Gene ID`)
xci_par=merge(xci,par,by="ensg",all=TRUE)%>%dplyr::select(ensg,`Gene name.x`,combined_xci_status,PAR_STATUS,PAR_BINARY)%>%rename(hgnc=`Gene name.x`)
chrx_rv_xci_par=merge(chrx_rv,xci_par,by.x="gene",by.y="ensg",all.x = TRUE)

ggplot(chrx_rv_xci_par%>%dplyr::filter((abs(FemaleMean)>.2|abs(MaleMean)>.2)&abs(diff)>.2),aes(x=PAR_STATUS,y=abs(diff)))+
  geom_violin()+
  geom_jitter(width = .1,height=0)+
  theme_bw()


ggplot(chrx_rv_xci_par%>%dplyr::filter((abs(FemaleMean)>.2|abs(MaleMean)>.2)&abs(diff)>.2),aes(x=combined_xci_status,y=abs(diff)))+
  geom_violin()+
  geom_jitter(width = .1,height=0)+
  theme_bw()


xci_par=xci_par%>%mutate(PAR_BINARY_NUMERIC=PAR_B)
ggplot(xci_par,aes(x=combined_xci_status))
