library("ggplot2")
library(tidyverse)
library("forcats")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("epitools")

maxmaf<-c("0.1", "0.05","0.01", "0.001","0.0001") #,"0.001","0.0001")
maxtomin_dic<-c("0.05","0.01","0.001","0","0"); # maxtomin_dic<-c("0","0","0","0","0","0")
names(maxtomin_dic)<-maxmaf


get_risks=function(collapsetypes="collapsed",windows=5000,max_outliers=3,filter_version="typesGQ5BlacklistRemovedALL",
                   risk_cat="relative_risk",maxtomin_dic,CADD=c(0,15),z=c(2,2.5,3),nphen=c(3),chr_types=c("x","AllAut","7sub"),my_cat="",
                   mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_all/RR",groups=c("m","f","both","allboth"),outlier_types="outliers" ){
  all_risks<-data.frame(matrix(nrow=0,ncol=10))
  colnames(all_risks)<-c("Risk","Lower","Upper" ,"Pval","Cat","Type" , "z","nphen","chr","sex")
  for(this_cat in my_cat){
    for(this_collapsetype in collapsetypes){
      for(this_outlier in outlier_types){
        for(this_maxmaf in maxmaf){
          print(this_maxmaf)
          for (this_group in groups){
            for(this_chr in chr_types){
              for(this_nphen in nphen){
                for(this_z in z){
                  for(this_filt in filter_version){
                    for(cadd_min in CADD){
                      for(this_window in windows){
                        if(this_chr=="x"){
                          myfile=paste0(mydir,"/relative_risk_x_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                      "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
                                      cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
                          this_filt="GQ5BlacklistRemovedALL"
                        }else if(grepl("sub",this_chr)){
                          myfile=paste0(mydir,"/relative_risk_",this_chr,"_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                      "_x_",this_group,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",
                                      cadd_min,"_linc_prot_maxoutliers",max_outliers,"_window",this_window,".csv")
                        } else{
                          myfile=paste0(mydir,"/relative_risk_",this_collapsetype,"_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                                      "_",this_group,"_",this_chr,"_min",maxtomin_dic[this_maxmaf],"max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
                                      "_linc_prot_maxoutliers",max_outliers,"_window",this_window,".txt")
                        }
                        #print(myfile)
                        if(!file.exists(myfile)){
                          print(paste0(myfile, " WAS NOT FOUND"))
                          next
                        }
                        risks=tryCatch({read.csv(myfile[1])},
                                       error=function(err){return(paste0("ERROR FILE NOT FOUND: ",myfile[1]))},
                                       warnings=function(war){return(paste0("ERROR FILE NOT FOUND: ",myfile[1]))})
                        
                        risks$z<-this_z
                        risks$nphen<-this_nphen
                        risks$chr<-this_chr
                        risks$sex<-this_group
                        risks$maxmaf<-this_maxmaf
                        risks$filt<-this_filt
                        risks$cadd<-cadd_min
                        risks$collapse_method<-this_collapsetype
                        risks$thiscat<-this_cat
                        risks$outlierType<-this_outlier
                        risks$gene_window<-this_window
                        risks$real_or_null<-"real"
                        risks$discrete_or_continuous<-"discrete"
                        colnames(risks)[6]<-"CATEGORY"
                        #assign(varname,risks)
                        #print(varname)
                        
                        all_risks<-rbind(all_risks,risks)
                      }
                      # break
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  all_risks$outliers_tested<-all_risks$exp_yn+all_risks$exp_yy
  all_risks$nonoutliers_tested<-all_risks$exp_nn+all_risks$exp_ny
  all_risks$prop_outliers<-all_risks$outliers_tested/all_risks$nonoutliers_tested
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/ following this
  risks_to_compare<-all_risks %>% group_by(outlierType,CATEGORY,z,nphen,cadd,chr,sex,maxmaf,real_or_null,
                                           exp_type,gene_window,var_location,veptype,collapse_method) %>% 
    dplyr::summarise(across(c(Risk,Lower,Upper),log)) %>%
    mutate(SE=abs(Lower-Upper)/(2*1.96))%>%dplyr::filter(sex=="m" | sex=="f")
  risks_to_compare_side<-pivot_wider(risks_to_compare,names_from=c(sex),
                                     id_cols=c(outlierType,CATEGORY,z,nphen,cadd,chr,maxmaf,real_or_null,exp_type,
                                               gene_window,var_location,veptype,collapse_method) ,
                                     values_from=c(Risk,SE)) #id_cols=, names_from=,
  risks_to_compare_side$RiskDiff<-risks_to_compare_side$Risk_f-risks_to_compare_side$Risk_m
  risks_to_compare_side$SEDiff<-sqrt(risks_to_compare_side$SE_f**2+risks_to_compare_side$SE_m**2)
  risks_to_compare_side$zDiff<-(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
  risks_to_compare_side$p<-pnorm(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
  #n is 2 for cadd * 3 for z * 3 for nphen * 3 for MAF * 24 for chr  (23 + allAut)=1296 OR 54* 2 for outliers
  risks_to_compare_side$padj<-p.adjust(risks_to_compare_side$p,method = "BH") #,n=108)
  all_risks_p<-merge(all_risks,risks_to_compare_side[,c("outlierType","exp_type","CATEGORY","z","nphen","cadd","chr",
                                                        "maxmaf","p","padj","gene_window","var_location","veptype","collapse_method")])
  all_risks_p$gene_window <- factor(all_risks_p$gene_window, levels = c("100","500","1000","5000","10000"))
  all_risks_p$outlierType <-str_replace(all_risks_p$outlierType,"outliers","")
  all_risks_p$Pval_adj<-p.adjust(all_risks_p$Pval,method = "BH",n = nrow(all_risks_p)) 
  all_risks_p=all_risks_p%>%mutate(significant=ifelse(Pval_adj<0.05,"significant","not significant"),
                                   maf_range=paste0("[",maxtomin_dic[maxmaf],"-",maxmaf,"]"))
    
  return(all_risks_p)
}



##VEPTYPE
all_risks=get_risks(collapsetypes="collapsed",windows=c("5000"),max_outliers=3,filter_version="typesGQ5BlacklistRemovedALL",
                    risk_cat="relative_risk",maxtomin_dic=maxtomin_dic,CADD=c(0,15),z=c(2.5),nphen=c(3),chr_types=c("7","x","AllAut"),my_cat="",
                    mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_all/RR",groups=c("m","f","both"))
all_risks=all_risks%>%mutate(chr=fct_relevel(chr,"x","7sub","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7sub"="chr7 subsetted","AllAut"="Autosomes","7"="chr7")) 
  

fig1<-function(all_risks){
  to_plot<-all_risks  %>% dplyr::filter(cadd==0&veptype=="all"&gene_window==5000&var_location=="all") %>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  all_risks_renamed<-all_risks
  chr_rename<-c("chrX","chr7","autosomes")
  names(chr_rename)<-c("x","7","AllAut")
  all_risks_renamed$chr_renamed<-chr_rename[all_risks_renamed$chr]
  
  
  
  # %>% dplyr::filter(exp_type=="under" | exp_type=="over") # %>%dplyr::filter(exp_type=="all") #
  to_plot_spread<-to_plot %>% dplyr::select(exp_type,chr,outlierType,maxmaf,outliers_tested,sex,z,nphen,cadd) %>% spread(exp_type,outliers_tested)
  to_plot_spread$over_prop=to_plot_spread$over/to_plot_spread$all
  to_plot_spread$under_prop=to_plot_spread$under/to_plot_spread$all
  to_plot_gather<-to_plot_spread %>% gather(exp_type,outliers_tested,all:under_prop)
  to_plot_gather_filt<-to_plot_gather %>% dplyr::filter(exp_type=="over_prop" | exp_type=="under_prop") %>% dplyr::filter(outlierType=="outliers")
  
  to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5&veptype=="all"&
                                                                                    gene_window==5000&var_location=="all")%>% 
    mutate(mygroup=ifelse(chr=="AllAut","b","a"))%>%
    mutate(chr=fct_relevel(chr,"x","7sub","AllAut")) %>% mutate(chr=recode(chr,x="chrX","7sub"="chr7 subsetted","AllAut"="Autosomes","7"="chr7"))  %>% 
    dplyr::filter(sex!="allboth") # %>%
    #dplyr::filter(chr=="autosomes") #dplyr::filter(chr=="autosomes") 
  ###f1a
  ggplot(to_plot%>%dplyr::filter(sex=="both"&chr!="Autosomes"&exp_type=="all"),aes(x=as.factor(chr),y=outliers_tested, group=sex,fill=sex ))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme_linedraw(base_size=15)+
    # facet_wrap(~mygroup,scales="free")+
    scale_alpha_discrete(range=c(1,0.7))+
    geom_bar(stat="identity", position=position_dodge())+
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
    ylab("Number of Outliers")+xlab("")
  ggplot(to_plot%>%dplyr::filter(sex=="both"&chr=="Autosomes"&exp_type=="all"),aes(x=as.factor(chr),y=outliers_tested, group=sex,fill=sex ))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme_linedraw(base_size=15)+
    # facet_wrap(~mygroup,scales="free")+
    scale_alpha_discrete(range=c(1,0.7))+
    geom_bar(stat="identity", position=position_dodge())+
    #scale_fill_manual(values=c("#dbab3b","#5d8596"))+  
    xlab("sex") +ylab("Number of Outliers")
  
  ###s1a
  ggplot(to_plot%>%dplyr::filter(exp_type!="all"&sex=="both"&chr=="Autosomes"),aes(x=chr,y=outliers_tested,group=sex,
                     alpha=as.factor(exp_type),
                     fill="#296818"))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"),guide="none")+
    theme_classic(base_size=10)+
    scale_alpha_discrete(range=c(0.7,0.5),name="")+
    geom_bar(stat="identity",position="stack",group="exp_type") +
    geom_text(size=5, aes(group=exp_type,label=outliers_tested),alpha=1, position = position_stack(vjust = .5))+
    # facet_wrap(~as.factor(chr),ncol=3,scales="free") +
    xlab("sex") +ylab("Number of Outliers")
  
  
  to_plot<-all_risks  %>% dplyr::filter(veptype=="all"&gene_window==5000&var_location=="all"&sex=="both") #%>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  
  ##f1, f1b, s1
  ggplot(to_plot%>%dplyr::filter(cadd==15&exp_type=="all"),
         aes(x=maf_range,y=Risk,color=chr,shape=chr,group=chr)) + #  shape=sex,
    scale_color_manual(values=c("#926fa8","#47265c","#4eb5ac"))+  
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))
  ggplot(to_plot%>%dplyr::filter(cadd==0&exp_type=="all"),aes(x=maf_range,y=Risk,color=chr,shape=chr,group=chr)) + #  shape=sex,
    scale_color_manual(values=c("#4eb5ac","#47265c","#926fa8"))+  
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+ggtitle("CADD=0")
  ggplot(to_plot%>%dplyr::filter(cadd==0&exp_type!="all"),aes(x=maf_range,y=Risk,color=chr,shape=exp_type,group=exp_type)) + #  shape=sex,
    scale_color_manual(values=c("#4eb5ac","#47265c","#926fa8"))+  
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    facet_wrap(~chr)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+ggtitle("CADD=0")
  ggplot(to_plot%>%dplyr::filter(cadd==15&exp_type!="all"),aes(x=maf_range,y=Risk,color=chr,shape=exp_type,group=exp_type)) + #  shape=sex,
    scale_color_manual(values=c("#4eb5ac","#47265c","#926fa8"))+  
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    facet_wrap(~chr)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+ggtitle("CADD=15")
  
  opentargets.6=fread("/oak/stanford/groups/smontgom/raungar/Sex/Data/opentargetscoresgt6.txt")
  colnames(opentargets.6)=c("target","ensg","score","NA")
  opentargets.6_dic=opentargets.6$score
  names(opentargets.6_dic)=opentargets.6$ensg
  gene_constraint=fread("/oak/stanford/groups/smontgom/shared/gnomadv4/gnomad.v4.0.constraint_metrics.tsv")
  gene_constraint$ensg=as.character(mapIds(org.Hs.eg.db,
                    keys=gene_constraint$gene, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first"))
  gene_constraint_missense_dic=gene_constraint%>%dplyr::filter(mane_select==TRUE)%>%group_by(ensg)%>%summarise(missense_z=max(mis.z_score))%>%pull(missense_z,ensg)
  gene_constraint_syn_dic=gene_constraint%>%dplyr::filter(mane_select==TRUE)%>%group_by(ensg)%>%summarise(syn_z=max(syn.z_score))%>%pull(syn_z,ensg)
  gene_constraint_lof_dic=gene_constraint%>%dplyr::filter(mane_select==TRUE)%>%group_by(ensg)%>%summarise(lof_z=max(lof.z_score))%>%pull(lof_z,ensg)
  
  dosage_xin_in=fread("/oak/stanford/groups/smontgom/raungar/Sex/Data/xin_dosage.csv")
  dosage_xin=dosage_xin_in%>%mutate(xin_mean=rowMeans(across(Adipose_Subcutaneous:Whole_Blood),na.rm = TRUE))%>%dplyr::select(ID,gene,xin_mean)%>%mutate(ensg=gsub("\\..*","",ID))
  dosage_anevadot=fread("/oak/stanford/groups/smontgom/raungar/Sex/Data/anevadot_gene_average.csv")
  dosage_xin_dic=dosage_xin$xin_mean
  names(dosage_xin_dic)=dosage_xin$ensg
  dosage_aneva_dic=dosage_anevadot$Avg_VG
  names(dosage_aneva_dic)=dosage_anevadot$GeneID
  
  
  my_zs_x=fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_all/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_x_allboth_maxoutliers3.txt.gz")
  my_zs_x=my_zs_x%>%mutate(ensg=gsub("\\..*","",Gene))%>%mutate(outlier_direction=ifelse(sign(MedZ)==1,"over",'under'))%>%
    mutate(open_target_score=opentargets.6_dic[ensg],gene_constraint_missense=gene_constraint_missense_dic[ensg],
           aneva_score=dosage_xin_dic[ensg],xin_score=dosage_aneva_dic[ensg])
  my_zs_aut=fread("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_all/OutliersFiltered/outliers_noglobal_medz_zthresh2.5_nphen3_aut_allboth_maxoutliers3.txt.gz")
  my_zs_aut=my_zs_aut%>%mutate(ensg=gsub("\\..*","",Gene))%>%mutate(outlier_direction=ifelse(sign(MedZ)==1,"over",'under'))%>%
    mutate(open_target_score=opentargets.6_dic[ensg],gene_constraint_missense=gene_constraint_missense_dic[ensg],
           gene_constraint_syn_dic=gene_constraint_syn_dic[ensg],gene_constraint_lof_dic=gene_constraint_lof_dic[ensg],
           aneva_score=dosage_xin_dic[ensg],xin_score=dosage_aneva_dic[ensg])
  
  run_rr=function(myzs,col_to_use,mydirection){
    # constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)>2)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)<(-2))%>%nrow()),
    #                                      c(myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)>2)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)<(-2))%>%nrow()))
    if(mydirection=="all"){
      constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&(get(col_to_use))>2)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&(get(col_to_use))<(-2))%>%nrow()),
                                           c(myzs%>%dplyr::filter(Y!="outlier"&(get(col_to_use))>2)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&(get(col_to_use))<(-2))%>%nrow()))
      
    }else{
      constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&(get(col_to_use))>2)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&(get(col_to_use))<(-2))%>%nrow()),
                                           c(myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&(get(col_to_use))>2)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&(get(col_to_use))<(-2))%>%nrow()))
      
    }

    
    colnames(constraint_vs_outlier_rr_table)=c("constrained","not_constrained"); rownames(constraint_vs_outlier_rr_table)=c("outlier","non-outlier")
    print(constraint_vs_outlier_rr_table%>%as.data.frame()%>%mutate(prop=constrained/(constrained+not_constrained)))
    constraint_vs_outlier_rr_table_risk = epitab(constraint_vs_outlier_rr_table, method = 'riskratio')
    constraint_vs_outlier_rr_res=data.frame(Risk = constraint_vs_outlier_rr_table_risk$tab[2,5],Lower = constraint_vs_outlier_rr_table_risk$tab[2,6],Upper = constraint_vs_outlier_rr_table_risk$tab[2,7],Pval = constraint_vs_outlier_rr_table_risk$tab[2,8])
    mytable_to_cols=t(data.frame(as.numeric(constraint_vs_outlier_rr_table)))
    colnames(mytable_to_cols)=c("constrained_outlier","constrained_nonoutlier","notconstrained_outlier","notconstrained_nonoutlier")
    rownames(mytable_to_cols)=""
    return(cbind(constraint_vs_outlier_rr_res,mytable_to_cols,mydirection,col_to_use))
  }
  run_rr_constraintscore=function(myzs,col_to_use,mydirection,lowscore,highscore){
    # constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)>2)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)<(-2))%>%nrow()),
    #                                      c(myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)>2)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)<(-2))%>%nrow()))
    if(mydirection=="all"){
      constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&get(col_to_use)>highscore)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&get(col_to_use)<(lowscore))%>%nrow()),
                                           c(myzs%>%dplyr::filter(Y!="outlier"&get(col_to_use)>highscore)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&get(col_to_use)<(lowscore))%>%nrow()))
      
    }else{
      constraint_vs_outlier_rr_table=rbind(c(myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)>highscore)%>%nrow(),myzs%>%dplyr::filter(Y=="outlier"&outlier_direction==mydirection&get(col_to_use)<(lowscore))%>%nrow()),
                                           c(myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)>highscore)%>%nrow(),myzs%>%dplyr::filter(Y!="outlier"&outlier_direction==mydirection&get(col_to_use)<(lowscore))%>%nrow()))
      
    }
    
    
    colnames(constraint_vs_outlier_rr_table)=c("constrained","not_constrained"); rownames(constraint_vs_outlier_rr_table)=c("outlier","non-outlier")
    print(constraint_vs_outlier_rr_table%>%as.data.frame()%>%mutate(prop=constrained/(constrained+not_constrained)))
    constraint_vs_outlier_rr_table_risk = epitab(constraint_vs_outlier_rr_table, method = 'riskratio')
    constraint_vs_outlier_rr_res=data.frame(Risk = constraint_vs_outlier_rr_table_risk$tab[2,5],Lower = constraint_vs_outlier_rr_table_risk$tab[2,6],Upper = constraint_vs_outlier_rr_table_risk$tab[2,7],Pval = constraint_vs_outlier_rr_table_risk$tab[2,8])
    
    mytable_to_cols=t(data.frame(as.numeric(constraint_vs_outlier_rr_table)))
    colnames(mytable_to_cols)=c("constrained_outlier","constrained_nonoutlier","notconstrained_outlier","notconstrained_nonoutlier")
    rownames(mytable_to_cols)=""
    return(cbind(constraint_vs_outlier_rr_res,mytable_to_cols,mydirection,col_to_use))
  }
  constraint_risks=rbind(run_rr(my_zs_aut,"gene_constraint_missense_dic","all"),
        run_rr(my_zs_aut,"gene_constraint_missense_dic","over"),
        run_rr(my_zs_aut,"gene_constraint_missense_dic","under"),
        run_rr(my_zs_aut,"gene_constraint_syn_dic","all"),
        run_rr(my_zs_aut,"gene_constraint_syn_dic","over"),
        run_rr(my_zs_aut,"gene_constraint_syn_dic","under"),
        run_rr(my_zs_aut,"gene_constraint_lof_dic","all"),
        run_rr(my_zs_aut,"gene_constraint_lof_dic","over"),
        run_rr(my_zs_aut,"gene_constraint_lof_dic","under"),
        run_rr_constraintscore(my_zs_aut,"aneva_score","all",.95,.05),
        run_rr_constraintscore(my_zs_aut,"aneva_score","under",.95,.05),
        run_rr_constraintscore(my_zs_aut,"aneva_score","over",.8,.2),
        run_rr_constraintscore(my_zs_aut,"xin_score","all",.95,.05),
        run_rr_constraintscore(my_zs_aut,"xin_score","under",.95,.05),
        run_rr_constraintscore(my_zs_aut,"xin_score","over",.85,.05))
  constraint_risks$padj=p.adjust(constraint_risks$Pval,method="BH")
  constraint_risks=constraint_risks%>%mutate(significant=ifelse(padj<0.05,"significant","not significant"))%>%
    mutate(renamed_col=ifelse(col_to_use=="gene_constraint_missense_dic","missense",
                              ifelse(col_to_use=="gene_constraint_syn_dic","synonymous",
                                     ifelse(col_to_use=="gene_constraint_lof_dic","LOF",
                                            ifelse(col_to_use=="aneva_score","ANEVA","MoD")))))
  constraint_risks$renamed_col=factor(constraint_risks$renamed_col,levels=c("synonymous","missense","LOF","ANEVA","MoD"))
  ggplot(constraint_risks,aes(x=renamed_col,y=Risk,color=mydirection,group=mydirection,shape=significant))+
    geom_hline(yintercept=1)+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.6)) +
    geom_point(position=position_dodge(width=0.6),aes(size=significant))+
    theme_bw(base_size=12)+
    scale_color_manual(values=c("#296818","#123807","#91bf84"),name="outlier direction")+xlab("constraint score")

  all_outliers=data.frame()
  for(mythresh in c(2,2.5,3)){
    for(mychr in c("x","aut")){
      this_outliers=fread(paste0("/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_all/OutliersFiltered/outliers_noglobal_medz_zthresh",mythresh,"_nphen3_",mychr,"_allboth_maxoutliers3.txt.gz"))
      this_outliers=this_outliers%>%dplyr::filter(Y=="outlier")%>%mutate(chr_type=mychr,z_thresh=mythresh)
      all_outliers=rbind(all_outliers,this_outliers)
    }
  }
  all_outliers_summ=all_outliers%>%mutate(direction=ifelse(MedZ>0,"over","under"))%>%
    mutate(chr_type=ifelse(chr=="chr7"&!is.na(chr),"7",chr_type))%>%
    group_by(chr_type,z_thresh,direction)%>%summarise(n=n())%>%
    mutate(chr_cat=ifelse(chr_type=="aut","a","b"))%>%
    mutate(my_x_axis=paste0(chr_type,"; ","z=",z_thresh))
  ggplot(all_outliers_summ,aes(x=my_x_axis,fill=as.factor(z_thresh),y=n,alpha=direction))+
    ylab("number of outliers")+xlab("outlier threshold")+
    scale_alpha_discrete(range=c(0.9,0.5),name="")+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c("#e39e66","#a86a38","#61340f"),name="z threshold")+
    facet_wrap(~chr_cat,scales="free")+
    theme_bw()
  
  all_risks=get_risks(collapsetypes="collapsed",windows=c("5000"),max_outliers=3,filter_version="typesGQ5BlacklistRemovedALL",
                      risk_cat="relative_risk",maxtomin_dic=maxtomin_dic,CADD=c(0,15),z=c(3,2.5,2),nphen=c(3),chr_types=c("x","7"),my_cat="",
                      mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_all/RR",groups=c("m","f","both"))
  to_plot<-all_risks  %>% dplyr::filter(veptype=="all"&gene_window==5000&var_location=="all"&sex=="both") #%>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  
  ##z params
  ggplot(to_plot%>%dplyr::filter(cadd==15&exp_type=="all"),
         aes(x=maf_range,y=Risk,color=as.factor(z),shape=chr,group=z)) + #  shape=sex,
    scale_color_manual(values=c("#e39e66","#a86a38","#61340f"),name="z threshold")+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    facet_wrap(~chr)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))
  
}



all_risks=get_risks(collapsetypes="collapsed",windows=c("5000"),max_outliers=3,filter_version="typesGQ5BlacklistRemovedALL",
                    risk_cat="relative_risk",maxtomin_dic=maxtomin_dic,CADD=c(0,15),z=c(2,2.5,3),nphen=c(3),chr_types=c("7","x","AllAut"),my_cat="",
                    mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR",groups=c("f","m","both"))
fig2<-function(all_risks){
  all_risks=get_risks(collapsetypes="collapsed",windows=c("5000"),max_outliers=3,filter_version="typesGQ5BlacklistRemovedALL",
                      risk_cat="relative_risk",maxtomin_dic=maxtomin_dic,CADD=c(0,15),z=c(2.5),nphen=c(3),chr_types=c("7","x","AllAut"),my_cat="",
                      mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8eqtl/RR",groups=c("f","m","both"))
  ##Plot by gene window
  
  to_plot<-all_risks  %>% dplyr::filter(cadd==0&veptype=="all"&gene_window==5000&var_location=="all") %>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  to_plot_spread<-to_plot %>% dplyr::select(exp_type,chr,outlierType,maxmaf,outliers_tested,sex,z,nphen,cadd) %>% spread(exp_type,outliers_tested)
  to_plot_spread$over_prop=to_plot_spread$over/to_plot_spread$all
  to_plot_spread$under_prop=to_plot_spread$under/to_plot_spread$all
  to_plot_gather<-to_plot_spread %>% gather(exp_type,outliers_tested,all:under_prop)
  to_plot_gather_filt<-to_plot_gather %>% dplyr::filter(exp_type=="over_prop" | exp_type=="under_prop") %>% dplyr::filter(outlierType=="outliers")
  to_plot=all_risks%>%dplyr::filter(CATEGORY=="all")%>%filter(cadd==15) %>%filter(maxmaf==0.01 & z==2.5&veptype=="all"& gene_window==5000&var_location=="all")%>% 
    mutate(mygroup=ifelse(chr=="AllAut","b","a"))%>%
    dplyr::filter(sex!="allboth") # %>%
  #f2,f2a
  ggplot(to_plot%>%dplyr::filter(exp_type!="all"),
         aes(x=sex,y=outliers_tested,group=sex,alpha=as.factor(exp_type),fill=sex))+
    scale_fill_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme_classic(base_size=10)+
    scale_alpha_discrete(range=c(0.9,0.5),name="")+
    geom_bar(stat="identity",position="stack",group="exp_type") +
    geom_text(size=5, aes(group=exp_type,label=outliers_tested),alpha=1, position = position_stack(vjust = .5))+
    facet_wrap(~as.factor(chr),ncol=3,scales="free") +
    xlab("sex") +ylab("Number of Outliers")
  
  to_plot<-all_risks  %>% dplyr::filter(cadd==15&veptype=="all"&gene_window==5000&var_location=="all"&exp_type=="all") #%>%dplyr::filter(sex!="both")#%>%dplyr::filter(CATEGORY=="all") # %>% dplyr::filter(outlierType=="outliersTOP") #$%>%dplyr::filter(CATEGORY=="all")
  
  ##f2, f2c, f2d, s2
  ggplot(to_plot%>%dplyr::filter(chr=="7"),aes(x=maf_range,y=Risk,group=sex,color=sex)) + #  shape=sex,
    scale_color_manual(values=c("#296818","#dbab3b","#5d8596"))+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1,position=position_dodge(width=0.5)) +
    geom_line(position=position_dodge(width=0.5))+
    theme_classic(base_size=15)+
    geom_point(aes(size=significant),position=position_dodge(width=0.5))+
    xlab("minor allele frequency bins")+ylab("relative risk")+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    ggtitle("Chromosome 7")
  
  
  ###gene window
  df_plot_gene_window=all_risks%>%dplyr::filter(veptype=="all"&var_location=="all"&exp_type=="all")%>%mutate(chr=ifelse(chr=="x","X-chromosome","Autosomes"))
  ggplot(df_plot_gene_window,aes(x=maf_range,y=Risk,color=gene_window,shape=significant,group=gene_window))+
    theme_bw()+
    facet_wrap(~chr*cadd,scales="free")+
    geom_line(position=position_dodge(width=.6))+
    scale_y_continuous(trans="log10")+
    scale_color_manual(values=c("#e35966","#87222c","#42040a"))+
    geom_errorbar(position=position_dodge(width=.6),aes(ymin=Lower,ymax=Upper),width=.1,size=.5)+
    geom_point(size=3,position=position_dodge(width=.6))
  
  
}



fig3<-function(mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output"){
#"outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh2_nphen3_aut_both_maxoutliers3.txt.gz
  sexes=c("f","m","both")
  chrtypes=c("x","aut")
  this_z=2.5
  all_zs=data.frame()
   for(this_sex in sexes){
     for(this_chr in chrtypes){
       myfile=paste0(mydir,"/outliers_v8eqtl/OutliersFiltered/outliers_noglobal_medz_zthresh",this_z,"_nphen3_",this_chr,"_",this_sex,"_maxoutliers3.txt.gz")
       mydf=fread(myfile)
       all_zs=rbind(all_zs,mydf%>%mutate(sex=this_sex,z=this_z,chrtype=this_chr))
     }
     
   }
  all_zs_wide=all_zs%>%pivot_wider(names_from = sex,values_from = c(MedZ,Y),id_cols=c(Ind,Gene,chr,chrtype)) #id_cols=c(Ind,Gene,chr,sex,chrtype),
  #all_zs_wide=reshape(all_zs,idvar=c("Ind","Gene","chr","chrtype"))
  all_zs_wide_filt=all_zs_wide%>%dplyr::filter((!is.na(MedZ_f)&!is.na(MedZ_both))|(!is.na(MedZ_m)&!is.na(MedZ_both)))
  to_plot_all_zs_wide_filt=all_zs_wide_filt%>% dplyr::filter(abs(MedZ_m)>2.5|abs(MedZ_f)>2.5|abs(MedZ_both)>2.5)%>%
                            mutate(between_sex_z=ifelse(is.na(MedZ_m),MedZ_f,MedZ_m),
                                   between_sex_Y=ifelse(is.na(Y_m),Y_f,Y_m),
                                   sex=ifelse(is.na(MedZ_m),"female","male"),
                                   chrtype=ifelse(chrtype=="x","X-Chromosome","Autosomes"),
                                   z_diff=between_sex_z-MedZ_both,
                                   bigger_z=ifelse(abs(between_sex_z)>abs(MedZ_both),between_sex_z,MedZ_both))
  to_plot_shuffled <- to_plot_all_zs_wide_filt[sample(nrow(to_plot_all_zs_wide_filt)),]
  ggplot(to_plot_shuffled,aes(x=between_sex_z,y=MedZ_both))+
    geom_hline(yintercept = -2.5,color="#bdbbbb")+geom_hline(yintercept=2.5,color="#bdbbbb")+geom_vline(xintercept=2.5,color="#bdbbbb")+geom_vline(xintercept=-2.5,color="#bdbbbb")+
    geom_abline(slope=1,intercept=0,color="#bdbbbb")+
    #stat_density_2d(aes(fill= ..level..), geom = 'polygon') +
    #stat_density_2d(data=to_plot_all_zs_wide_filt%>%dplyr::filter(abs(between_sex_z)<2.5&abs(MedZ_both)<2.5),aes(fill = stat(density)), geom = 'raster', contour = FALSE)+
    facet_wrap(~chrtype,scales="free")+
    scale_fill_viridis_c()+
    geom_point(aes(color=sex),alpha=.4)+
    scale_color_manual(values=c("#dbab3b","#5d8596"))+ #"#296818","#dbab3b","#5d8596"
    xlab("sex-stratified z-score")+ylab("combined sex z-score")+
    theme_bw(base_size=14)
  #supp supplemental table outlier change
  supp_table_out_change=to_plot_all_zs_wide_filt%>%dplyr::filter(between_sex_Y!=Y_both)%>%
                        dplyr::select(Gene,chrtype,chr,between_sex_z,MedZ_both,z_diff,between_sex_Y,Y_both,sex)%>%
                        rename("chromosome"="chr",both_z=MedZ_both,between_sex_outlier_status=between_sex_Y,both_outlier_status=Y_both)
  fwrite(x=supp_table_out_change,file="/oak/stanford/groups/smontgom/raungar/Sex/Files/supp_table_outliers.txt",sep="\t")
  
  ##FAM9C
  FAM9C=all_zs%>%dplyr::filter(grepl("ENSG00000187268",Gene)&(sex=="f"|sex=="both"))
  ggplot(FAM9C,aes(x=MedZ,color=sex))+
    scale_color_manual(values=c("#296818","#dbab3b"))+ #"#296818","#dbab3b","#5d8596"
    geom_density()+
    xlab("multi-tissue z-score")+ylab("density")+
    geom_point(size=4,data=FAM9C%>%dplyr::filter(Ind=="GTEX-1GF9X"),aes(y=0))+
    theme_bw(base_size = 18)+ggtitle("FAM9C")
  IL1RAPL2=all_zs%>%dplyr::filter(grepl("ENSG00000189108",Gene))%>%dplyr::filter((sex=="m"|sex=="both"))
  ggplot(IL1RAPL2,aes(x=MedZ,color=sex))+
    scale_color_manual(values=c("#296818","#5d8596"))+ #"#296818","#dbab3b","#5d8596"
    geom_density()+
    xlab("multi-tissue z-score")+ylab("density")+
    geom_point(size=4,data=IL1RAPL2%>%dplyr::filter(Ind=="GTEX-13O3Q"),aes(y=0))+
    theme_bw(base_size = 18)+ggtitle("IL1RAPL2")
  G6PC2=all_zs%>%dplyr::filter(grepl("ENSG00000152254",Gene))%>%dplyr::filter((sex=="f"|sex=="both"))
  ggplot(G6PC2,aes(x=MedZ,color=sex))+
    scale_color_manual(values=c("#296818","#dbab3b"))+ #"#296818","#dbab3b","#5d8596"
    geom_density()+
    xlab("multi-tissue z-score")+ylab("density")+
    geom_point(size=4,data=G6PC2%>%dplyr::filter(Ind=="GTEX-1GMR8"),aes(y=0))+
    theme_bw(base_size = 18)+ggtitle("G6PC2")
  
  ##supplemental table
  supp_table_outliers=to_plot_all_zs_wide_filt%>%dplyr::filter(Y_both!=between_sex_Y)
  supp_table_outliers=supp_table_outliers%>%dplyr::select(Ind,Gene,chrtype,chr,MedZ_both,Y_both,between_sex_z,between_sex_Y,sex,z_diff)%>%
    rename(within_sex_z=MedZ_both,within_sex_outlier_status=Y_both,between_sex_outlier_status=between_sex_Y)
  
}

