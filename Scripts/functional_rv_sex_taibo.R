library(data.table)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(clusterProfiler)
library(readxl)
library(VennDiagram)
library(biomaRt)


chrx_rv_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/chrXResult_Merged.RDS"
aut_rv_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/AutosomeResult_Merged.RDS"
model_f_aut_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Female_AllAutosomes_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_m_aut_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Male_AllAutosomes_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_f_x_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Female_X_Eval_AutosomeAnnotations_SameModel.txt.gz"
model_m_x_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Male_X_Eval_AutosomeAnnotations_SameModel.txt.gz"
posterior_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/Combined_SexSpecificPosterior.csv.gz"
oliva_eqtl_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/oliva-eqtl-table-s11.xlsx"
jones_vars_file="/oak/stanford/groups/smontgom/raungar/Sex/Data/jones_s7.xlsx"

sexbias_vars="/oak/stanford/groups/smontgom/raungar/Sex/Data/sexbiased_vars_list.tsv.gz"

get_posteriors<-function(posterior_file,filter_min=0){
  posteriors=fread(posterior_file)
  print(nrow(posteriors))
  posteriors=posteriors%>%mutate(permute_num_diff=as.numeric((abs(Permute_1)>.2)+(abs(Permute_2)>.2)+(abs(Permute_3)>.2)+(abs(Permute_4)>.2)+(abs(Permute_5)>.2)))
  posteriors_filt=posteriors%>%dplyr::filter(permute_num_diff>=filter_min&abs(diff)>.2)
  return(posteriors_filt)
}
#reformat_model_df=function(filename,my_over,my_under,my_not_sexbiased,sex){
reformat_model_df=function(filename,myvarlist,sex){
  my_over=myvarlist%>%dplyr::filter(over_under=="over")
  my_under=myvarlist%>%dplyr::filter(over_under=="under")
  my_not_sexbiased=myvarlist%>%dplyr::filter(bias=="not")
  
    mymodel=fread(filename)
  mymodel=mymodel%>%rename(beta=p_value)
  mymodel$variant=gsub("(.*)chr","chr",(mymodel$SubjectID)) #sep variant into its own col
  mymodel=mymodel%>%mutate(geneRV=paste0(GeneName,"+",variant))
  #filtering, using datatable framework for speed
  setkey(mymodel,geneRV)
  mymodel_over=mymodel[.(my_over$geneRV), nomatch = 0L] %>%mutate(category=paste0(sex,"-biased (over-expression)"))
  mymodel_under=mymodel[.(my_under$geneRV), nomatch = 0L] %>%mutate(category=paste0(sex,"-biased (under-expression)"))
  mymodel_neutral=mymodel[.(my_not_sexbiased$geneRV), nomatch = 0L] %>%mutate(category=paste0("non-biased"))
  
  return(rbind(mymodel_over,mymodel_under,mymodel_neutral))
  
}
get_significance_comparisons<-function(mymodel,num_mymodel,mymodel_baseline,num_mymodel_baseline){
  mymodel_df=cbind(data.frame(variant_present=mymodel),nrow_mymodel=num_mymodel)%>%
    mutate(not_present=nrow_mymodel-mymodel)%>%
    select(-nrow_mymodel)
  mymodel_df_baseline=cbind(data.frame(variant_present_baseline=mymodel_baseline),nrow_mymodelbaseline=num_mymodel_baseline)%>%
    mutate(not_present_baseline=nrow_mymodelbaseline-mymodel_baseline)%>%select(-nrow_mymodelbaseline)
  combine_nums=cbind(mymodel_df,mymodel_df_baseline)
  pvals=apply(combine_nums,1,
              function(a){
                contig_table=cbind(rbind(as.numeric(a[1]),as.numeric(a[2])),
                                   rbind(as.numeric(a[3]),as.numeric(a[4])))
                colnames(contig_table)=c("sex_biased","not_sex_biased")
                rownames(contig_table)=c("variant_present","not_present")
                fisher_test_res=fisher.test(contig_table)
                fisher_test_res_df=t(data.frame(c(fisher_test_res$p.value,fisher_test_res$estimate,fisher_test_res$conf.int)))
                colnames(fisher_test_res_df)=c("pval","point_estimate","CI_lower","CI_upper")
                fisher_test_res_df=as.data.frame(fisher_test_res_df)
                #if small, fishers exact; else, chi-square test
                # if(sum(a)<20){
                #   my_pval=fisher.test(contig_table)$p.value
                # }else{
                #   my_pval=chisq.test(contig_table)
                # }
              })
  final_df=do.call(rbind.data.frame,pvals)
  final_df$pval_adj=p.adjust(final_df$pval,method="BH")
  final_df$variant=rownames(final_df)
  return(final_df)
}


aut_rv=readRDS(aut_rv_file)%>%mutate(gene=gsub("\\+.*","",geneRV))
chrx_rv=readRDS(chrx_rv_file)%>%mutate(gene=gsub("\\+.*","",geneRV))
chr7_rv=(aut_rv%>%dplyr::filter(grepl("chr7",geneRV)))%>%mutate(gene=gsub("\\+.*","",geneRV))

#replicated at least once
confident_sexbias_vars=get_posteriors(posterior_file,filter_min=1)%>%mutate(gene=gsub("\\+.*","",geneRV))
#write_confident_sexbias_vars=confident_sexbias_vars%>%select(gene,RV,diff,Permute_1,Permute_2,Permute_3,Permute_4,Permute_5,permute_num_diff)
#fwrite(

aut_rv_filt=aut_rv%>%dplyr::filter(geneRV %in%confident_sexbias_vars$geneRV&(abs(MaleMean)>.2|abs(FemaleMean)>.2))
#fwrite(aut_rv_filt,file=sexbias_vars,quote=FALSE,sep="\t")
aut_rv_all=aut_rv%>%mutate(sexbiased_pass=ifelse((geneRV %in%confident_sexbias_vars$geneRV),1,0))
x_rv_filt=chrx_rv%>%dplyr::filter(geneRV %in%confident_sexbias_vars$geneRV)




#nothing enriched
get_go_enrichment<-function(genelist){
  ego <- enrichGO(gene          =genelist, 
                  ont ="ALL", 
                  keyType = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  return(ego)
}
# my_enrich_female=get_go_enrichment(aut_rv_filt%>%dplyr::filter(abs(FemaleMean)>.2))
# my_enrich_male=get_go_enrichment(aut_rv_filt%>%dplyr::filter(abs(MaleMean)>.2))
# my_enrich=get_go_enrichment(aut_rv_filt)

get_eqtl_overlap<-function(oliva_eqtl_file,aut_rv_filt){
  eqtls=read_xlsx(oliva_eqtl_file,sheet=2)
  eqtls$ensg=gsub("\\..*","",eqtls$ENSEMBL_gene_id)
  
  jones_vars=readxl::read_xlsx(jones_vars_file,sheet = 1,skip=2)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  jones_eqtls=getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'),
                     filters = 'hgnc_symbol', 
                     values = unique(c(jones_vars$symbol)),
                     mart = ensembl)
  
  venn_list=list(eqtls_oliva=unique(eqtls$ensg),eqtls_jones=unique(jones_eqtls$ensembl_gene_id),sex_biased=unique(aut_rv_filt$gene))
  vennplot=venn.diagram(venn_list, alpha = c(0.6, 0.6,0.6),fill=c("#cf619b","#851104","#c26c1b"),filename=NULL); grid.newpage(); grid.draw(vennplot)
  
  eqtls_oliva_len=length(unique(eqtls$ensg))
  eqtls_jones_len=length(unique(jones_eqtls$ensembl_gene_id))
  eqtls_sexbias_len=length(unique(aut_rv_filt$gene))
  #phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

  phyper(length(intersect(unique(eqtls$ensg),unique(jones_eqtls$gene))),
         eqtls_jones_len,35065-eqtls_jones_len,jones_eqtls,lower.tail = FALSE)
  
  phyper(length(intersect(unique(eqtls$ensg),unique(aut_rv_filt$gene))),
         eqtls_jones_len,35065-eqtls_jones_len,eqtls_sexbias_len,lower.tail = FALSE)
  
  phyper(length(intersect(unique(jones_eqtls$ensembl_gene_id),unique(aut_rv_filt$gene))),
         eqtls_jones_len,35065-eqtls_jones_len,eqtls_sexbias_len,lower.tail = FALSE)
  intersect(unique(jones_eqtls$ensembl_gene_id),unique(aut_rv_filt$gene),unique(eqtls$ensg))
}




get_num_variants<-function(df,chr){
  rv_sexbiased=df%>%dplyr::filter(abs(diff)>.2& (abs(MaleMean)>.2|abs(FemaleMean)>.2))
  not_sexbiased=df%>%dplyr::filter(abs(diff)<0.1 & (abs(MaleMean)>.2|abs(FemaleMean)>.2)) %>%mutate(category=paste0("non-biased"),bias="not",over_under="none")#CHANGE 0.1, MEAN>0.2
  rv_male_over=rv_sexbiased%>%dplyr::filter((MaleMean>FemaleMean)&sign(MaleMean)==1&MaleMean>0.2)%>%mutate(category=paste0("male-biased (over-expression)"),bias="male",over_under="over")
  rv_male_under=rv_sexbiased%>%dplyr::filter((MaleMean<FemaleMean)&sign(MaleMean)==-1&MaleMean<(-0.2))%>%mutate(category=paste0("male-biased (under-expression)"),bias="male",over_under="under")
  rv_female_over=rv_sexbiased%>%dplyr::filter((FemaleMean)>MaleMean&sign(FemaleMean)==1&FemaleMean>0.2)%>%mutate(category=paste0("female-biased (over-expression)"),bias="female",over_under="over")
  rv_female_under=rv_sexbiased%>%dplyr::filter((FemaleMean<MaleMean)&sign(FemaleMean)==-1&FemaleMean<(-0.2))%>%mutate(category=paste0("female-biased (under-expression)"),bias="female",over_under="under")
  print(paste0("num male over/under: ", nrow(rv_male_over),"/",nrow(rv_male_under), " (total = ",nrow(rv_male_under)+nrow(rv_male_over),")",
               " ; num female over/under: ",nrow(rv_female_over),"/",nrow(rv_female_under)," (total = ",nrow(rv_female_over)+nrow(rv_female_under),")"))
  outdf=data.frame(
    rbind(
      c("male",chr,"variants",nrow(rv_male_over),nrow(rv_male_under),nrow(rv_male_under)+nrow(rv_male_over)),
      c("female",chr,"variants",nrow(rv_female_over),nrow(rv_female_under),nrow(rv_female_under)+nrow(rv_female_over)),
      
      c("male",chr,"genes",nrow(rv_male_over%>%group_by(gene)%>%summarise(n=n())),nrow(rv_male_under%>%group_by(gene)%>%summarise(n=n())),nrow(rv_male_under%>%group_by(gene)%>%summarise(n=n()))+nrow(rv_male_over%>%group_by(gene)%>%summarise(n=n()))),
      c("female",chr,"genes",nrow(rv_female_over%>%group_by(gene)%>%summarise(n=n())),nrow(rv_female_under%>%group_by(gene)%>%summarise(n=n())),nrow(rv_female_under%>%group_by(gene)%>%summarise(n=n()))+nrow(rv_female_over%>%group_by(gene)%>%summarise(n=n())))
    ))
  colnames(outdf)=c("sex","chr","var_or_genes","over","under","both")
  
  out_all_df=data.frame(rbind(not_sexbiased,rv_female_over,rv_female_under,rv_male_over,rv_male_under))
  return(list(outdf,out_all_df))
}
aut_num_sexbias_filt=get_num_variants(aut_rv_filt,"aut")
aut_num_sexbias=get_num_variants(aut_rv_all,"aut")
chrx_num_sexbias=get_num_variants(chrx_rv,"chrX")
chr7_num_sexbias=get_num_variants(chr7_rv,"chr7")

#num_sexbias_df=rbind(aut_num_sexbias[[1]],chrx_num_sexbias[[1]],chr7_num_sexbias[[1]])
num_sexbias_df=aut_num_sexbias[[1]]
##sup fig
num_sexbias_df_melted=melt(num_sexbias_df,id.vars = c('sex','chr','var_or_genes'))%>%dplyr::filter(variable!="both")
num_sexbias_df_melted$value=as.numeric(num_sexbias_df_melted$value)
ggplot(num_sexbias_df_melted%>%dplyr::filter(var_or_genes=="genes"),aes(x=sex,y=value,fill=variable))+
  geom_bar(stat = "identity")+
  facet_wrap(~chr,scales="free_y")+
  ggtitle("number of genes with different predicted functional effect by sex")+
  theme_bw()


###Fig4a
# all_sexbiased=rbind(cbind(aut_rv_female_over,sex="female-biased",over_under="over",mean_to_use=aut_rv_female_over$FemaleMean,chr="autosomes"),
#                     cbind(aut_rv_female_under,sex="female-biased",over_under="under",mean_to_use=aut_rv_female_under$FemaleMean,chr="autosomes"),
#                     cbind(aut_rv_male_over,sex="male-biased",over_under="over",mean_to_use=aut_rv_male_over$MaleMean,chr="autosomes"),
#                     cbind(aut_rv_male_under,sex="male-biased",over_under="under",mean_to_use=aut_rv_male_under$MaleMean,chr="autosomes"),
#                     cbind(x_rv_female_over,sex="female-biased",over_under="over",mean_to_use=x_rv_female_over$FemaleMean,chr="x-chromosome"),
#                     cbind(x_rv_female_under,sex="female-biased",over_under="under",mean_to_use=x_rv_female_under$FemaleMean,chr="x-chromosome"),
#                     cbind(x_rv_male_over,sex="male-biased",over_under="over",mean_to_use=x_rv_male_over$MaleMean,chr="x-chromosome"),
#                     cbind(x_rv_male_under,sex="male-biased",over_under="under",mean_to_use=x_rv_male_under$MaleMean,chr="x-chromosome"))
ggplot(aut_num_sexbias_filt[[2]],aes(x=FemaleMean,y=MaleMean,color=bias))+
  geom_point(alpha=0.5)+
  xlim(c(-1,1))+
  xlab("mean posterior (female)")+ylab("mean posterior (male)")+
  scale_color_manual(values=c("#dbab3b","#5d8596"))+
  #facet_wrap(~chr)+
  theme_bw(base_size=12)
###Fig4b
ggplot(aut_num_sexbias_filt[[2]],aes(x=diff,group=interaction(bias,over_under),fill=bias,height=..density..))+
  geom_density(stat = "density",trim=T)+
  #facet_wrap(~chr)+
  xlab("difference in mean posterior")+
  scale_fill_manual(values=alpha(c("#dbab3b","#5d8596"),0.7),"sex bias")+
  geom_vline(xintercept=-0.2,color="#7a7878")+geom_vline(xintercept=0.2,color="#7a7878")+
  theme_bw(base_size=12)

###Fig 4c

model_f_aut=reformat_model_df(model_f_aut_file,aut_num_sexbias[[2]],"female")
model_m_aut=reformat_model_df(model_m_aut_file,aut_num_sexbias[[2]],"male")
model_all_aut=rbind(cbind(model_f_aut,sex="female"),cbind(model_m_aut,sex="male"))
model_f_x=reformat_model_df(model_f_x_file,chrx_num_sexbias[[2]],"female")
model_m_x=reformat_model_df(model_m_x_file,chrx_num_sexbias[[2]],"male")
model_all_x=rbind(model_f_x,model_m_x)


get_props_functional_bias=function(model_all,this_chr,confident_sexbias_vars,num_diff=0){
  sexbiased=confident_sexbias_vars%>%dplyr::filter(permute_num_diff>=num_diff)
  
  model_f_over_df=(model_all%>%dplyr::filter(sex=="female"&sign(beta)>0&geneRV%in%sexbiased$geneRV)); model_f_over=colSums(model_f_over_df[,c(3:28)])
  model_f_under_df=(model_all%>%dplyr::filter(sex=="female"&sign(beta)<0&geneRV%in%sexbiased$geneRV)); model_f_under=colSums(model_f_under_df[,c(3:28)])
  model_f_df=(model_all%>%dplyr::filter(sex=="female"&geneRV%in%sexbiased$geneRV)); model_f=colSums(model_f_df[,c(3:28)])
  model_m_over_df=(model_all%>%dplyr::filter(sex=="male"&sign(beta)>0&geneRV%in%sexbiased$geneRV)); model_m_over=colSums(model_m_over_df[,c(3:28)])
  model_m_under_df=(model_all%>%dplyr::filter(sex=="male"&sign(beta)<0&geneRV%in%sexbiased$geneRV)); model_m_under=colSums(model_m_under_df[,c(3:28)])
  model_m_df=(model_all%>%dplyr::filter(sex=="male"&geneRV%in%sexbiased$geneRV)); model_m=colSums(model_m_df[,c(3:28)])
  model_unbiased_df=model_all%>%dplyr::filter(abs(beta)>.2&!(geneRV%in%sexbiased$geneRV)); model_unbiased=colSums(model_unbiased_df[,c(3:28)])
  model_unbiased_over_df=(model_all%>%dplyr::filter(beta>.2&!(geneRV%in%sexbiased$geneRV))); model_unbiased_over=colSums(model_unbiased_over_df[,c(3:28)])
  model_unbiased_under_df=(model_all%>%dplyr::filter(beta<(-0.2)&!(geneRV%in%sexbiased$geneRV)));  model_unbiased_under=colSums(model_unbiased_under_df[,c(3:28)])
  model_unbiased_df_nonfunctional=model_all%>%dplyr::filter(abs(beta)<.2&!(geneRV%in%sexbiased$geneRV)); model_unbiased_nonfunctional=colSums(model_unbiased_df_nonfunctional[,c(3:28)])
  
  #functional: biased vs unbiased
  sig_biased_vs_unbiased_f=get_significance_comparisons(model_f,nrow(model_f_df),
                                                        model_unbiased,nrow(model_unbiased_df))
  sig_biased_vs_unbiased_f_over=get_significance_comparisons(model_f_over, nrow(model_f_over_df),
                                                            model_unbiased,nrow(model_unbiased_df))
  sig_biased_vs_unbiased_f_under=get_significance_comparisons(model_f_under,nrow(model_f_under_df),
                                                              model_unbiased,nrow(model_unbiased_df)) #%>%as.data.frame()%>%dplyr::filter(.<0.05)
  sig_biased_vs_unbiased_m=get_significance_comparisons(model_m,nrow(model_m_df),
                                                        model_unbiased,nrow(model_unbiased_df))
  sig_biased_vs_unbiased_m_over=get_significance_comparisons(model_m_over,nrow(model_m_over_df),
                                                             model_unbiased,nrow(model_unbiased_df))
  sig_biased_vs_unbiased_m_under=get_significance_comparisons(model_m_under,nrow(model_m_under_df),
                                                              model_unbiased,nrow(model_unbiased_df))
  #function vs nonfunction
  sig_biased_vs_nonfunctional_f=get_significance_comparisons(model_f,nrow(model_f_df),
                                                        model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  sig_biased_vs_nonfunctional_f_over=get_significance_comparisons(model_f_over, nrow(model_f_over_df),
                                                             model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  sig_biased_vs_nonfunctional_f_under=get_significance_comparisons(model_f_under,nrow(model_f_under_df),
                                                              model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  sig_biased_vs_nonfunctional_m=get_significance_comparisons(model_m,nrow(model_m_df),
                                                        model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  sig_biased_vs_nonfunctional_m_over=get_significance_comparisons(model_m_over,nrow(model_m_over_df),
                                                             model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  sig_biased_vs_nonfunctional_m_under=get_significance_comparisons(model_m_under,nrow(model_m_under_df),
                                                              model_unbiased_nonfunctional,nrow(model_unbiased_df_nonfunctional))
  
  sig_f_vs_m_over=get_significance_comparisons(model_f_over,nrow(model_f_over_df),
                                                 model_m_over,nrow(model_m_over_df))
  sig_f_vs_m_under=get_significance_comparisons(model_f_under,nrow(model_f_under_df),
                                                  model_m_under,nrow(model_m_under_df))
  sig_f_vs_m=get_significance_comparisons(model_f,nrow(model_f_df),
                                                model_m,nrow(model_m_df))
  
  sig_df=rbind(cbind(sig_biased_vs_unbiased_f_over,test="biased vs. unbiased",sex="f",chr=this_chr,over_under="over"),
                 cbind(sig_biased_vs_unbiased_f_under,test="biased vs. unbiased",sex="f",chr=this_chr,over_under="under"),
               cbind(sig_biased_vs_unbiased_m_over,test="biased vs. unbiased",sex="m",chr=this_chr,over_under="over"),
               cbind(sig_biased_vs_unbiased_m_under,test="biased vs. unbiased",sex="m",chr=this_chr,over_under="under"),
               cbind(sig_biased_vs_unbiased_m,test="biased vs. unbiased",sex="m",chr=this_chr,over_under="all"),
               cbind(sig_biased_vs_unbiased_f,test="biased vs. unbiased",sex="f",chr=this_chr,over_under="all"),
               cbind(sig_biased_vs_nonfunctional_f_over,test="biased vs. nonfunctional",sex="f",chr=this_chr,over_under="over"),
               cbind(sig_biased_vs_nonfunctional_f_under,test="biased vs. nonfunctional",sex="f",chr=this_chr,over_under="under"),
               cbind(sig_biased_vs_nonfunctional_m_over,test="biased vs. nonfunctional",sex="m",chr=this_chr,over_under="over"),
               cbind(sig_biased_vs_nonfunctional_m_under,test="biased vs. nonfunctional",sex="m",chr=this_chr,over_under="under"),
               cbind(sig_biased_vs_nonfunctional_m,test="biased vs. nonfunctional",sex="m",chr=this_chr,over_under="all"),
               cbind(sig_biased_vs_nonfunctional_f,test="biased vs. nonfunctional",sex="f",chr=this_chr,over_under="all"),
               cbind(sig_f_vs_m,test="male vs. female",sex=NA,chr=this_chr,over_under="all"),
               cbind(sig_f_vs_m_over,test="male vs. female",sex=NA,chr=this_chr,over_under="over"),
               cbind(sig_f_vs_m_under,test="male vs. female",sex=NA,chr=this_chr,over_under="under"))
  # heatmap_df=rbind(model_unbiased/nrow(model_all),
  #                    model_f_over/nrow(model_all%>%dplyr::filter(category=="female-biased (over-expression)")),
  #                    model_f_under/nrow(model_all%>%dplyr::filter(category=="female-biased (under-expression)")),
  #                    model_m_over/nrow(model_all%>%dplyr::filter(category=="male-biased (over-expression)")),
  #                    model_m_under/nrow(model_all%>%dplyr::filter(category=="male-biased (under-expression)")))
  # heatmap_df=rbind(model_unbiased/nrow(model_all%>%dplyr::filter(abs(beta)>.2&!(geneRV%in%sexbiased$geneRV))),
  #                  model_f/nrow(model_all%>%dplyr::filter(sex=="female"&geneRV%in%sexbiased$geneRV)),
  #                  model_m/nrow(model_all%>%dplyr::filter(sex=="male"&geneRV%in%sexbiased$geneRV)))
  heatmap_df=rbind(model_unbiased_nonfunctional/nrow(model_unbiased_df_nonfunctional),
                   model_unbiased_under/nrow(model_unbiased_under_df),
                   model_unbiased_over/nrow(model_unbiased_over_df),
                   model_f_under/nrow(model_f_under_df),
                   model_f_over/nrow(model_f_over_df),
                   model_m_under/nrow(model_m_under_df),
                   model_m_over/nrow(model_f_over_df))
  heatmap_df_all=cbind(model_unbiased_nonfunctional/nrow(model_unbiased_df_nonfunctional),
                  model_unbiased/nrow(model_unbiased_df),
                   model_f/nrow(model_f_df),
                   model_m/nrow(model_m_df))
  colnames(heatmap_df_all)=c(paste0("unbiasedANDnonfunctional_",this_chr,"_all"),
                             paste0("unbiased_",this_chr,"_all"),
                             paste0("f_",this_chr,"all"),
                             paste0("m_",this_chr,"all"))
  # heatmap_df_melted=heatmap_df[,colSums(heatmap_df)>0.001]%>%melt()
  heatmap_df_t=t(heatmap_df)
  colnames(heatmap_df_t)=c(paste0("unbiasedANDnonfunctional_",this_chr,"_all"),
                           paste0("unbiased_",this_chr,"_under"),paste0("nonbiased_",this_chr,"_over"),
                           paste0("f_",this_chr,"_under"),paste0("f_",this_chr,"_over"),
                           paste0("m_",this_chr,"_under"),paste0("m_",this_chr,"_over"))
  heatmap_df_subtractUnbiased=cbind(heatmap_df_t[,2]-heatmap_df_t[,1],
                                      heatmap_df_t[,3]-heatmap_df_t[,1],
                                      heatmap_df_t[,4]-heatmap_df_t[,1],
                                      heatmap_df_t[,5]-heatmap_df_t[,1],
                                    heatmap_df_t[,6]-heatmap_df_t[,1],
                                    heatmap_df_t[,7]-heatmap_df_t[,1])
  colnames(heatmap_df_subtractUnbiased)=colnames(heatmap_df_t)[-1]
  return(list(heatmap_df,heatmap_df_subtractUnbiased,sig_df,not_subtracted=heatmap_df_t,heatmap_df_all))
}
aut_get_props_functional_bias_res=get_props_functional_bias(model_all_aut,"aut",confident_sexbias_vars,1)
heatmap_df_aut_subtractUnbiased=aut_get_props_functional_bias_res[[2]]
heatmap_df_aut=aut_get_props_functional_bias_res[[1]]
aut_sig_df=aut_get_props_functional_bias_res[[3]]
heatmap_df_aut_nonsubtracted=aut_get_props_functional_bias_res[[4]]
heatmap_df_aut_nonsubtracted_whole=aut_get_props_functional_bias_res[[5]]

x_get_props_functional_bias_res=get_props_functional_bias(model_all_x,"x")
heatmap_df_x_subtractUnbiased=x_get_props_functional_bias_res[[2]]
heatmap_df_x=x_get_props_functional_bias_res[[1]]
x_sig_df=x_get_props_functional_bias_res[[3]]

# heatmap_df_subtractUnbiased_melted=rbind(t(heatmap_df_x_subtractUnbiased)[,colSums(heatmap_df_x_subtractUnbiased)>0.0001]%>%melt(),
#                                          t(heatmap_df_aut_subtractUnbiased)[,colSums(heatmap_df_aut_subtractUnbiased)>0.0001]%>%melt())%>%
#   separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)
# heatmap_df_subtractUnbiased_melted=rbind(t(heatmap_df_x_subtractUnbiased)%>%melt(),
#                                          t(heatmap_df_aut_subtractUnbiased)%>%melt())%>%
#   separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)%>%
#   dplyr::filter(!(variant%in%c("incomplete_terminal_codon_variant","coding_sequence_variant","inframe_insertion","TFBS_ablation",
#                                "stop_retained_variant","stop_lost", "inframe_deletion","frameshift_variant")))
heatmap_df_subtractUnbiased_melted= t(heatmap_df_aut_subtractUnbiased)[,colSums(heatmap_df_aut_subtractUnbiased)>0.0001]%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)
heatmap_df_subtractUnbiased_melted=t(heatmap_df_aut_subtractUnbiased)%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)%>%
  dplyr::filter(!(variant%in%c("incomplete_terminal_codon_variant","coding_sequence_variant","inframe_insertion","TFBS_ablation",
                               "stop_retained_variant","stop_lost", "inframe_deletion","frameshift_variant")))

heatmap_df_nonsubtracted_melted= t(heatmap_df_aut_nonsubtracted)[,colSums(heatmap_df_aut_nonsubtracted)>0.0001]%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)
heatmap_df_nonsubtracted_melted=t(heatmap_df_aut_nonsubtracted)%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)%>%
  mutate(over_under=ifelse(is.na(over_under),"all",over_under)) %>%
  dplyr::filter(!(variant%in%c("incomplete_terminal_codon_variant","coding_sequence_variant","inframe_insertion","TFBS_ablation",
                              "stop_retained_variant","stop_lost", "inframe_deletion","frameshift_variant")))

heatmap_df_nonsubtracted_whole_melted= t(heatmap_df_aut_nonsubtracted_whole)[,colSums(heatmap_df_aut_nonsubtracted_whole)>0.0001]%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)
heatmap_df_nonsubtracted_whole_melted=t(heatmap_df_aut_nonsubtracted_whole)%>%melt()%>%
  separate(Var1,c("sex","chr","over_under"),"_")%>%rename(variant=Var2)%>%
  mutate(over_under=ifelse(is.na(over_under),"all",over_under)) %>%
  dplyr::filter(!(variant%in%c("incomplete_terminal_codon_variant","coding_sequence_variant","inframe_insertion","TFBS_ablation",
                               "stop_retained_variant","stop_lost", "inframe_deletion","frameshift_variant")))


ggplot(heatmap_df_nonsubtracted_whole_melted,mapping=aes(x=sex,y=variant,fill=value))+
  geom_tile()+
  #facet_wrap(~over_under,scales = "free_x")+
  xlab("")+ylab("")+
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(low = "#8f2e21", high = "#c26c1b",midpoint=0) +
  ggtitle("Proportion of variants in with given variant type")+
  theme_bw()
  # theme(axis.text.x = element_text(angle = 45, hjust=1))

## significance plot
all_sig_df=rbind(aut_sig_df)%>%mutate(significance=ifelse(pval_adj<0.05,"significant","not significant")) %>%dplyr::filter(pval_adj<0.05)
#all_sig_df=rbind(aut_sig_df,x_sig_df)%>%mutate(significance=ifelse(pval_adj<0.05,"significant","not significant")) %>%dplyr::filter(pval_adj<0.05)
ggplot(all_sig_df%>%dplyr::filter(test=="biased vs. unbiased"),aes(x=reorder(variant,point_estimate),y=point_estimate,color=sex))+
  facet_wrap(~interaction(over_under),scales = "free_y",nrow = 3)+
  geom_hline(yintercept=1,color="gray")+
  #scale_color_manual(values=c("black","gray","#dbab3b","#5d8596","#dbc591","#a4cdde"))+
  scale_color_manual(values=c("#dbab3b","#5d8596"))+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2,position=position_dodge(width=0.9))+
  geom_point(aes(size=significance),position=position_dodge(width=0.9))+
  theme_bw(base_size = 8)+
  scale_y_continuous(trans='log10')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("")+ylab("fisher's test")+
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  ggtitle("biased vs nonbiased")
ggplot(all_sig_df%>%dplyr::filter(test=="biased vs. nonfunctional"),
       aes(x=reorder(variant,point_estimate),y=point_estimate,color=(sex)))+
  facet_wrap(~over_under,scales = "free_y",ncol=1)+
  geom_hline(yintercept=1,color="gray")+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values=c("#dbab3b","#5d8596"))+
  #scale_color_manual(values=c("#570644","#8c497c","black"))+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2,position=position_dodge(width=0.9))+
  geom_point(aes(size=significance),position=position_dodge(width=0.9))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("")+ylab("fisher's test")+
  ggtitle("biased vs non-functional")

# scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))


ggplot(chrx_rv_filt,aes(x=MaleMean,y=FemaleMean,alpha=0.1,size=abs(diff),color=sex_bias))+
  geom_point()+
  theme_bw()+
  xlim(c(-.85,0.85))+ylim(c(-.85,0.85))+
  scale_size(limits=c(0,max(chrx_rv_filt$diff)))+
  scale_color_manual(values=c("yellow","blue"))+
  xlab("mean beta: male")+ylab("mean beta: female")+
  ggtitle("X-Chromosome")




