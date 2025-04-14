library(ggplot2)
library(tidyverse)
library(data.table)

plot_outliers<-function(props,mytitle,fisher_test_res){
  merged_props_fisher=merge(props,fisher_test_res,by=c("vartype","chr","direction"))
  merged_props_fisher$star=ifelse(merged_props_fisher$signficant=="significant"&merged_props_fisher$outlier_status=="outlier","*","")
  g=ggplot(merged_props_fisher,aes(x=vartype,y=prop,fill=outlier_status))+
    geom_bar(stat='identity',position = 'dodge')+
    theme_bw(base_size = 14)+ 
    scale_fill_manual(values=c("#ff8fd8","#f21da8"),name="outlier status")+
    scale_y_continuous(trans='sqrt',limits = c(0,1))+
    facet_wrap(~direction,ncol=1)+
    ylab('proportion of genes')+
    xlab("variant type")+
    geom_text(aes(label=star),size=15)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(mytitle)
  return(g)
}
read_file<-function(filename){
  props=fread(filename)
  colnames(props)<-c("outlier_status","vartype","chr","num","num_ind_genes","prop","direction")
  props_to_plot=props%>%mutate(vartype=ifelse(vartype=='','no_variant',vartype))%>%mutate(outlier_status=ifelse(outlier_status=="control","non-outlier","outlier"))
  props_to_plot=props%>%mutate(vartype=ifelse(vartype=='','no_variant',vartype))%>%mutate(outlier_status=ifelse(outlier_status=="control","non-outlier","outlier"))
  return(props_to_plot)
}
mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_all/VarProp"
myfile_aut=paste0(mydir,"/aut_collapsed_outliers_zthresh2.5_nphen3_allboth_CADDtypesALL_CADD15_linc_prot_maxoutliers3_window5000_maf020.01.varprops.txt.gz")
props_aut=read_file(myfile_aut)
props_aut=props_aut%>%dplyr::group_by(outlier_status,vartype)%>%summarise(num=sum(num),num_status=sum(num_status))%>%mutate(prop=num/num_status)%>%mutate(chr="aut")
myfile_x=paste0(mydir,"/x_collapsed_outliers_zthresh2.5_nphen3_allboth_CADDtypesGQ5BlacklistRemovedALL_CADD15_linc_prot_maxoutliers3_window5000_maf020.01.varprops.txt.gz")
props_x=read_file(myfile_x)
props_all=rbind(props_aut,props_x) %>%mutate(chr=ifelse(chr=="x","x","aut"))%>%group_by(outlier_status,direction,vartype,chr)%>%summarise(num=sum(num),num_ind_genes=sum(num_ind_genes))%>%mutate(prop=num/num_ind_genes)
  #mutate(num_with=prop*num_status,num_without=num_status-(prop*num_status))

library('epitools')

get_tab_outlier_vs_nonoutlier<-function(myprops,variant_to_filt,this_chr,this_direction){
  this_var=myprops%>%dplyr::filter(vartype==variant_to_filt&chr==this_chr&direction==this_direction)%>%
    mutate(num_with=num,num_without=num_ind_genes-num_with)
  if(nrow(this_var)<=1){return("only in 1")}
  tab=as.data.frame(this_var[,c("num_with","num_without")])
  rownames(tab)<-this_var$outlier_status
  return(tab)
}


get_rr<-function(tab){
  epi=epitab(tab,method="riskratio")
  return(epi$tab[2,5:ncol(epi$tab)])
}


get_fisher_test<-function(myprop){
  fisher_df_out_vs_nonout<-data.frame()
  for (vartype in unique(myprop$vartype)){
    for(this_direction in unique(myprop$direction)){
      for(this_chr in unique(myprop$chr)){
        props_all_dir=myprop%>%dplyr::filter(direction==this_direction)
        mytab=get_tab_outlier_vs_nonoutlier(props_all_dir,vartype,this_chr,this_direction)
        
        if(is.null(dim(mytab))){
          this_line=c(test="outlier_vs_nonoutlier",vartype=vartype,chr=this_chr,pval=NA,odds_ratio=NA,ci_lower=NA,ci_upper=NA,
                      out_with=NA,out_without=NA,nonout_with=NA,nonout_without=NA,direction=this_direction)
        }else{
          sig_level=fisher.test(mytab)
          this_line=c(test="outlier_vs_nonoutlier",vartype=vartype,chr=this_chr,pval=sig_level$p.value,"odds_ratio"=as.character(sig_level$estimate),
                      ci_lower=sig_level$conf.int[1],ci_upper=sig_level$conf.int[2],
                      out_with=mytab['outlier','num_with'],out_without=mytab['outlier','num_without'],nonout_with=mytab['non-outlier','num_with'],nonout_without=mytab['non-outlier','num_without'],
                      direction=this_direction)
        }
        fisher_df_out_vs_nonout<-rbind(fisher_df_out_vs_nonout,as.data.frame(t(this_line)))
      }
    }
  }
  fisher_df_out_vs_nonout=fisher_df_out_vs_nonout%>%mutate(pval_adj=p.adjust(pval,method="BH"),
                                                           signficant=ifelse(pval_adj<0.05,"significant","not significant"))
  return(fisher_df_out_vs_nonout)
}

fisher_df_out_vs_nonout_x=get_fisher_test(props_x)
fisher_df_out_vs_nonout_aut=get_fisher_test(props_aut%>%dplyr::filter(chr=="chr7"))
fisher_df_out_vs_nonout_all=get_fisher_test(props_all)
### SUPPLEMENTAL FIGURE 2 S2
plot_outliers(props_aut%>%dplyr::filter(chr=="chr7"),"chromosome 7",fisher_df_out_vs_nonout_aut)
plot_outliers(props_x,"X-chromosome",fisher_df_out_vs_nonout_x)
plot_outliers(props_all%>%dplyr::filter(chr=="aut"),"Autosomes",fisher_df_out_vs_nonout_all)





fisher_df_x_vs_aut<-data.frame()

for (vartype in unique(props_all$vartype)){
  for(outstatus in unique(props_all$outlier_status)){
    mytab=get_tab_aut_vs_chr(props_all,vartype,outstatus)
    
    if(is.null(dim(mytab))){
      this_line=c(test="aut_vs_x",vartype=vartype,outlier_status=outstatus,pval=NA,odds_ratio=NA,ci_lower=NA,ci_upper=NA,
                  aut_with=NA,aut_without=NA,x_with=NA,x_without=NA)
    }else{
      sig_level=fisher.test(mytab)
      this_line=c(test="aut_vs_x",vartype=vartype,outlier_status=outstatus,pval=sig_level$p.value,"odds_ratio"=as.character(sig_level$estimate),
                  ci_lower=sig_level$conf.int[1],ci_upper=sig_level$conf.int[2],
                  aut_with=mytab['aut','num_with'],aut_without=mytab['aut','num_without'],x_with=mytab['x','num_with'],x_without=mytab['x','num_without'])
    }
    fisher_df_x_vs_aut<-rbind(fisher_df_x_vs_aut,as.data.frame(t(this_line)))
  }
}

fisher_df_x_vs_aut=fisher_df_x_vs_aut%>%mutate(adj.pval=ifelse(is.na(pval),NA,as.numeric(pval)*nrow(fisher_df_x_vs_aut)))%>%
  mutate(is_sig=ifelse(adj.pval<0.05,"significant","non-significant")) #%>%dplyr::filter(adj.pval<0.05)

ggplot(fisher_df_x_vs_aut%>%dplyr::filter(!is.na(pval)),
       aes(x=vartype,y=as.numeric(odds_ratio),group=outlier_status,color=outlier_status,shape=is_sig))+
  geom_hline(yintercept = 1,color='gray')+
  geom_point(position=position_dodge(width=0.9))+
  geom_errorbar(aes(ymin=as.numeric(ci_lower),ymax=as.numeric(ci_upper)),width=0.3,size=0.2,position=position_dodge(width=0.9))+
  scale_y_continuous(trans="log10")+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#9A7AA0","#702632"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("x vs. aut")+ylab("Odds Ratio")+xlab("")
