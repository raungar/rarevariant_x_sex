library(data.table) #for fread
library(ggplot2)

dir="/oak/stanford/groups/smontgom/raungar/Sex/Output/outliers_v8"

groups=c("both","m","f","both.sex","both.regress")
types=c("aut","x")
for (this_group in groups){
  for(this_type in types){
    this_file=paste0(dir,"/outliers_zthresh2_nphen5_globalOutliersRemoved_",
                 this_type,"_",this_group,".txt")
    varname=paste0(this_type,"_",this_group)
    print(this_file)
    df1<-fread(this_file)
    df2<-cbind(df1,"varname"=varname)
    assign(varname,df2)

  }
}
get_melted_df<-function(group1,group2){
  group1_cast<-dcast(group1[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  group1_cast_med<-group1_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group2_cast<-dcast(group2[,c(1,2,5,7)],  Gene + Ind ~varname, value.var="MedZ")
  group2_cast_med<-group2_cast[ , lapply(.SD, median), by = Gene,.SDcols=3]
  group1vsgroup2<-merge(group1_cast_med,group2_cast_med)
  # group1vsgroup2<-rbind(group1,group2)[,c(1,2,5,7)]
  # group1vsgroup2_cast<-dcast(group1vsgroup2,  Gene + Ind ~varname, value.var="MedZ")
  # #group1vsgroup2_cast_no_na<-na.omit(group1vsgroup2_cast)
  # group1vsgroup2_cast_no_na_med<-group1vsgroup2_cast_no_na[ , lapply(.SD, median), by = Gene,.SDcols=3:4]
  return(group1vsgroup2)
}

x_m_bothregress<-get_melted_df(x_m,x_both.regress)
x_f_bothregress<-get_melted_df(x_f,x_both.regress)
x_m_f<-get_melted_df(x_m,x_f)
x_bothregress_both<-get_melted_df(x_both.regress,x_both)
x_bothsex_both<-get_melted_df(x_both.sex,x_both)
x_bothregress_bothsex<-get_melted_df(x_both.regress,x_both.sex)


aut_m_bothregress<-get_melted_df(aut_m,aut_both.regress)
aut_f_bothregress<-get_melted_df(aut_f,aut_both.regress)
aut_m_f<-get_melted_df(aut_m,aut_f)
aut_bothregress_both<-get_melted_df(aut_both.regress,aut_both)
aut_bothsex_both<-get_melted_df(aut_both.sex,aut_both)
aut_bothregress_bothsex<-get_melted_df(aut_both.regress,aut_both.sex)
  
highz<-(apply(x_f_bothregress,1,function(x){any(abs(as.numeric(x[2:3]))>.5)}))
highz_genes<-x_f_bothregress[highz,]
highz_genes

plot_list<-list("x_m_bothregress"=x_m_bothregress,"x_f_bothregress"=x_f_bothregress,"x_m_f"=x_m_f,
                "x_bothregress_both"=x_bothregress_both,"x_bothsex_both"=x_bothsex_both,"x_bothregress_bothsex"=x_bothregress_bothsex,
                "aut_m_bothregress"=aut_m_bothregress,"aut_f_bothregress"=aut_f_bothregress,"aut_m_f"=aut_m_f,
                "aut_bothregress_both"=aut_bothregress_both,"aut_bothsex_both"=aut_bothsex_both,"aut_bothregress_bothsex"=aut_bothregress_bothsex)
i<-0
for (this_plot in plot_list){
  i<-i+1
  this_name<-names(plot_list)[i]
  ggplot(this_plot,
         aes_string(x=colnames(this_plot)[2], y=colnames(this_plot)[3], alpha=0.2))+
     theme(legend.position = "none")+
    ggtitle(paste("Median Z-score Gene Across Tissues: ",this_name))+
    geom_abline(intercept = 0, slope = 1,color="blue",alpha=0.2)+
    geom_abline(intercept = 0, slope = 0,color="blue",alpha=0.2)+
    geom_vline(xintercept = 0,color="blue",alpha=0.2)+
    #xlim(c(-1.05,1.05))+ylim(-1.3,1.3)+
    geom_point()
  ggsave(paste0("/oak/stanford/groups/smontgom/raungar/Sex/Plots/outliers_v8/MedzCompare/",this_name,"_genemedz.png"),
         width=10,height=10)
}
