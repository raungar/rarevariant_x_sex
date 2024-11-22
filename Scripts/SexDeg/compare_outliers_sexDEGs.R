library(data.table) #for fread
library(ggplot2)
library(matrixStats)

sexdegfile<-c("/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/combined_sexdegs.txt.gz")
sexdegs<-fread(sexdegfile)
sexdegs_minmax<-sexdegs[,`:=`(
  MIN= rowMins(as.matrix(.SD), na.rm=T),
  MAX= rowMaxs(as.matrix(.SD), na.rm=T)
), .SDcols=colnames(sexdegs)[-c(1,2)]]


t_sexdegs_minmax<-cvi
sexdegs_minmax_abs<-cbind(sexdegs_minmax,ABS=apply(sexdegs_minmax[,-c(1,2)],1,function(x) max(abs(x))))
dic_genes_beta_f<-sexdegs_minmax_abs$MAX
names(dic_genes_beta_f)<-sexdegs_minmax_abs$ENSG
dic_genes_beta_m<-sexdegs_minmax_abs$MIN
names(dic_genes_beta_m)<-sexdegs_minmax_abs$ENSG
dic_genes_beta_both<-sexdegs_minmax_abs$ABS
names(dic_genes_beta_both)<-sexdegs_minmax_abs$ENSG

dir="/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8"

groups=c("m","f","both_half.regress")
types=c("x") ##add in aut yo

for (this_group in groups){
  for(this_type in types){
    this_file=paste0(dir,"/outliers_zthresh3_nphen5_globalOutliersRemoved_",
                     this_type,"_",this_group,".txt")
    varname=paste0(this_type,"_",this_group)
    print(this_file)
    df1<-fread(this_file)
    df2<-cbind(df1,"varname"=varname,
               "beta_max"=dic_genes_beta_f[(df1)$Gene],
               "beta_min"=dic_genes_beta_m[(df1)$Gene],
               "beta_abs"=dic_genes_beta_both[(df1)$Gene])
    assign(varname,df2)
    
  }
}

get_melted_df<-function(group1,group2){
  group1_cast<-dcast(group1[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  # group1_cast_med<-group1_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group2_cast<-dcast(group2[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  # group2_cast_med<-group2_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group1vsgroup2<-merge(group1_cast,group2_cast)
  group1vsgroup2_full<-cbind(group1vsgroup2,
        "beta_max"=dic_genes_beta_f[(group1vsgroup2)$Gene],
        "beta_min"=dic_genes_beta_m[(group1vsgroup2)$Gene],
        "beta_abs"=dic_genes_beta_both[(group1vsgroup2)$Gene])
  return(group1vsgroup2_full)
}
x_m_both<-get_melted_df(x_m,x_both_half.regress)
x_f_both<-get_melted_df(x_f,x_both_half.regress)
x_m_f<-get_melted_df(x_m,x_f)

this_plot<-x_m_both
ggplot(this_plot,
       aes_string(x=colnames(this_plot)[2], y=colnames(this_plot)[3],
                  color=colnames(this_plot)[6], alpha=0.7))+
  theme(legend.position = "none")+
  ggtitle(paste0("Median Z-score Across Tissues: "))+
  geom_abline(intercept = 0, slope = 1,color="blue",alpha=0.2)+
  geom_abline(intercept = 0, slope = 0,color="blue",alpha=0.2)+
  geom_vline(xintercept = 0,color="blue",alpha=0.2)+
  # xlim(c(-6,10))+
  # ylim(-6,9)+
   scale_color_gradient2(midpoint = 0, low="pink", high="red",mid="gray")+
  geom_point()

