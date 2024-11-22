Output/expression_v8/Counts/normalized_melted-m-protect_sex-nullshuffled_v8-both-binary.txt.gz


library("ggplot2")
library("data.table")
library("tidyverse")

sex="both"
method_type=c("protect_sex") #c("protect_sex","no_sex","incl_sex")
null_or_real=c("preprocessing_v8")
dir="/Volumes/groups/smontgom/raungar/Sex/Output/expression_v8/Counts"
chrtype="aut"
binary_cont=c("binary","continuous_alpha0.5")
tissue="Breast_Mammary_Tissue"

all_data<-data.table()

for(this_tissue in tissue){
  for(this_sex in sex){
    for(this_method in method_type){
      for(this_null in null_or_real){
        for(this_chrtype in chrtype){
          for(this_bc in binary_cont){
            file=paste0(dir,"/counts_melted-",this_tissue,"-",this_sex,"-",this_method,"-",this_null,"-",this_chrtype,"-",this_bc,".txt.gz")
            data=fread(file)
            # head(data)
            this_data=data[,c("Id","Ind","Counts")]
            colnames(this_data)[3]<-as.character(data[1,"sex_cont"])
            # tmp<-head(all_data)[,c("Id","Ind","Counts")]
            # tmp_data<-head(this_data)
            # inner_join(tmp,tmp_data)
            if(nrow(head(all_data))<2){
              print("creating new ")
              all_data=this_data
              print(head(all_data))
              }
            else{
              print("joining")
              all_data<-inner_join(all_data,this_data)
              print(head(all_data))
              }
              
          }
          
        }
      }
    }
  }
}


# head(all_data)
# to_plot<-all_data %>% dplyr::filter(regression_method=="incl_sex" | regression_method=="protect_sex") %>% select(Id,Ind,Counts,regression_method)
# to_plot_cast<-dcast.data.table(data=to_plot,id.vars=c("Id","Ind","Counts"),regression_method~Counts)
# to_plot_melt<-melt.data.table(to_plot,id.vars=c("Id"),variable.name = 'Ind', value.name = 'Counts')
to_plot<-all_data %>% dplyr::filter(Ind=="GTEX-T6MN")
ggplot(to_plot,aes(x=T,F))+geom_point()

all_data$diff<-all_data$F-all_data$T

all_data_lm<-lm(F~T,all_data)
summ<-summary(all_data_lm) #Blood, r2=0.95
