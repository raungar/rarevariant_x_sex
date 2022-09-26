library("ggplot2")

file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"
load(file)

if(exists("all_coefs")){
  my_df<-all_coefs
}
if(exists("risks")){
  my_df<-risks
}

my_df_aut_both<-risks

my_df_final_aut<-rbind(cbind(my_df_aut_m,pop="aut_m"),
  cbind(my_df_aut_f,pop="aut_f"),
  cbind(my_df_aut_both,pop="aut_both"))

my_df_final_x<-rbind(cbind(my_df_x_m,pop="x_m"),
                   cbind(my_df_x_f,pop="x_f"),
                   cbind(my_df_x_both,pop="x_both"))
my_df_final<-rbind(cbind(my_df_final_x,chr="x"),
                   cbind(my_df_final_aut,chr="aut"))
#melted_df<-melt(my_df[4:6,])
ggplot(my_df_final,aes(fill=pop)) + 
  scale_fill_manual(values=c("#3CAFD2","#3C6AD2","#8B3CD2",
                       "#A8C784","#71A679","#226F3A"))+
  ylim(c(0,8))+
  #geom_bar(aes(x=Cat,y=Risk),stat="identity",position="dodge")+
  geom_crossbar(aes(x=Cat,y=Risk,ymin=Lower,ymax=Upper),position="dodge")
 