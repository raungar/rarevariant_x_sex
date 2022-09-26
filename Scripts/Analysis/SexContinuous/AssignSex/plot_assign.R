library("ggplot2")
library(tidyverse)
library(data.table)
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

mydir="/Volumes/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/Assignments"
md_file="/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
md=fread(md_file)

sample_id_file<-paste0(mydir,"/test_train.txt")
sample_id<-fread(sample_id_file,data.table=F)
colnames(sample_id)<-c("Id","sex","ntiss","group","scrambled_sex","scrambled_id")
sample_to_sex<-md$SEX
names(sample_to_sex)<-as.character(md$SUBJID)

tissue=c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala",
"Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex",
"Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
"Brain_Spinal_cord_cervical_c1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Colon_Sigmoid", #"Cells_EBV",
"Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle",
"Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Pancreas","Pituitary","Skin_Not_Sun_Exposed_Suprapubic",
"Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Thyroid","Whole_Blood")
#Artery_Tibial-preprocessing_v8-accuracy-aut.txt 
counts_group=c("preprocessing_v8") #,"nullshuffled_v8")
counts_group_dic=c("original","null")
names(counts_group_dic)<-c("preprocessing_v8") #,"nullshuffled_v8")
region=c("x","aut","both")
all_accuracy<-data.frame(matrix(nrow=0,ncol=17))
all_preds<-data.frame(matrix(nrow=0,ncol=6))
tissue=c("Cells_Cultured_fibroblasts")


for(this_tiss in tissue){
  for(this_counts_group in counts_group){
    for(this_region in region){
      accuracy_file=paste0(mydir,"/",this_tiss,"-",this_counts_group,"-accuracy-",this_region,".txt")
      this_accuracy<-fread(accuracy_file,data.table=F)
      this_accuracy$tissue<-this_tiss
      this_accuracy$counts_group<-counts_group_dic[as.character(this_counts_group)]
      this_accuracy$subgroup<-this_region
      all_accuracy<-rbind(all_accuracy,this_accuracy)
      
      preds_file=paste0(mydir,"/",this_tiss,"-",this_counts_group,"-preds-",this_region,".txt")
      this_pred<-fread(preds_file)
      colnames(this_pred)[1]<-"Ind"
      this_pred_melt<-melt.data.table(this_pred)
      colnames(this_pred_melt)<-c("Ind","alpha","sex")
      this_pred_melt$tissue<-this_tiss
      this_pred_melt$subgroup<-this_region
      this_pred_melt$counts_group<-counts_group_dic[as.character(this_counts_group)]
      
      all_preds<-rbind(all_preds,this_pred_melt)
      
      
    }
  }

}


all_preds$sex_binary<-sample_to_sex[as.character(all_preds$Ind)]

##accuracy
to_plot<-all_accuracy %>%dplyr::filter(alpha==0.5) # & subgroup=="both")
ggplot(to_plot,aes(x=tissue,y=test_accuracy,color=counts_group))+
 theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_color_manual(values=c("orange","red","#ff85c2"))+  
  ggtitle(paste0("Test Accuracy for alpha=",unique(to_plot$alpha)))+
  facet_wrap(~paste0("region=",subgroup),ncol =1,scales="free_y")+geom_point()
  geom_jitter(aes(alpha=.5,size=3),width=0.3)

##preds
to_plot<-all_preds %>% dplyr::filter(tissue=="Cells_Cultured_fibroblasts")  %>% dplyr::filter(subgroup=="aut")# %>% dplyr::filter(sex_binary==1) %>%dplyr::filter(alpha==0.5)
ggplot(to_plot,aes(x=Ind,y=sex,color=interaction(as.factor(sex_binary),counts_group)))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #scale_color_manual(values=c("orange","red","#ff85c2"))+  
  theme_bw()+
  scale_color_discrete(name=guide_legend(title="Sex and Null/Original"))+
  ggtitle(paste0("Sex score across alphas"))+
  facet_wrap(~paste0("region=",subgroup),ncol =1,scales="free_y") + geom_boxplot()
  #geom_violin()

preds_summ <- all_preds[, as.list(summary(sex)), by=c("alpha","subgroup","counts_group","Ind","sex_binary")]
preds_summ_alpha0 <- (all_preds %>% dplyr::filter(alpha==0))[, as.list(summary(sex)), by=c("alpha","subgroup","counts_group","Ind","sex_binary")]
preds_summ_alpha0.5 <- (all_preds %>% dplyr::filter(alpha==0.5))[, as.list(summary(sex)), by=c("alpha","subgroup","counts_group","Ind","sex_binary")]
preds_summ_alpha1 <- (all_preds %>% dplyr::filter(alpha==1))[, as.list(summary(sex)), by=c("alpha","subgroup","counts_group","Ind","sex_binary")]
preds_summ_med_med<-preds_summ[,as.list(summary(Median)),by=c("subgroup","counts_group","Ind","sex_binary")]
to_plot<-preds_summ_med_med # %>% dplyr::filter(subgroup=="both") #%>% dplyr::filter(alpha==1)
ggplot(to_plot,aes(x=Median,fill=as.factor(sex_binary),alpha=counts_group,group=sex_binary))+
  scale_alpha_discrete(range=c(0.3,0.7))+
  ggtitle(paste0("Median of medians score across alpha for the ",unique(to_plot$subgroup)))+
  geom_density()+
  facet_wrap(~subgroup)
  # geom_histogram(bins=50)
#good inds: GTEX-OOBJ  and GTEX-1211K 
to_plot<-all_preds %>% dplyr::filter(tissue=="Cells_Cultured_fibroblasts")  %>% dplyr::filter(Ind=="GTEX-1211K")# %>% dplyr::filter(sex_binary==1) %>%dplyr::filter(alpha==0.5)
ggplot(to_plot,aes(x=Ind,y=sex,color=interaction(as.factor(sex_binary),counts_group)))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #scale_color_manual(values=c("orange","red","#ff85c2"))+  
  theme_bw()+
  scale_color_discrete(name=guide_legend(title="Sex and Null/Original"))+
  ggtitle(paste0("Sex score across alphas"))+ylim(c(0,1))+
  facet_wrap(~paste0("region=",subgroup),ncol =1,scales="free_y") + geom_boxplot()
 #geom_violin() 

###md correlate
inds=all_preds %>%pull(Ind) %>%unique()
md_red<-md%>% dplyr::filter(SUBJID %in% inds)
reduced_covars<-function(covariates){
  
  #convert hrs/mins string to mins numeric
  covariates$TRISCH<-as.numeric(str_split_fixed(covariates$TRISCH," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRISCH," ",n=4)[,3])
  covariates$TRCHSTIN<-as.numeric(str_split_fixed(covariates$TRCHSTIN," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRCHSTIN," ",n=4)[,3])
  covariates$TRCCLMP<- as.numeric(str_split_fixed(covariates$TRCCLMP," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRCCLMP," ",n=4)[,3])
  covariates$DTHPRNINT<- as.numeric(str_split_fixed(covariates$DTHPRNINT," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$DTHPRNINT," ",n=4)[,3])
  
  #if didn't smoke, smoking years to zero not NA like why yo
  covariates[is.na(covariates$MHSMKYRS) & covariates$MHSMKSTS=="No","MHSMKYRS"]<-0
  covariates$MHSMKYRS<-as.numeric(covariates$MHSMKYRS)
  covariates[(covariates$MHDRNKNMB==99) & covariates$MHDRNKSTS=="Yes","MHDRNKYRS"]<-"NA"
  covariates[(covariates$MHDRNKSTS=="No"),"MHDRNKYRS"]<-0
  covariates$MHDRNKYRS<-as.numeric(covariates$MHDRNKYRS)
  covariates[(covariates$MHCOPD==99) & !is.na(covariates$MHCOPD),"MHCOPD"]<-"NA"
  covariates$MHCOPD<-as.numeric(covariates$MHCOPD)
  covariates[(covariates$MHBCTINF==99) & !is.na(covariates$MHBCTINF),"MHBCTINF"]<-"NA"
  covariates$MHBCTINF<-as.numeric(covariates$MHBCTINF)
  covariates[(covariates$MHNPHYS4W==99) & !is.na(covariates$MHNPHYS4W),"MHNPHYS4W"]<-"NA"
  covariates$MHNPHYS4W<-as.numeric(covariates$MHNPHYS4W)
  
  
  covariates$SEX<-as.numeric(covariates$SEX-1)
  
  
  # covariates$INCEXC<-as.factor(covariates$INCEXC)
  # covariates$DTHDTRMN<-as.factor(covariates$DTHDTRMN)
  # covariates$COHORT<-as.factor(covariates$COHORT)
  # covariates$TRCRTMPU<-as.factor(covariates$TRCRTMPU)
  # covariates$TRTPTREF<-as.factor(covariates$TRTPTREF)
  # covariates$TRVNTSR<-as.factor(covariates$TRVNTSR)
  # covariates$DTHTPTREF<-as.factor(covariates$DTHTPTREF)
  # covariates$DTHMNNR<-as.factor(covariates$DTHMNNR)
  # covariates$DTHRFGD<-as.factor(covariates$DTHRFGD)
  # covariates$DTHPLCE<-as.factor(covariates$DTHPLCE)
  # covariates$MHSRC<-as.factor(covariates$MHSRC)
  # covariates$DTHSEASON<-as.factor(covariates$DTHSEASON)
  # 
  covariates_red<-covariates #[,c(2:7,9,11:18,21,22,24:32,34,37,39,41,43,45:47,49:77,79:97,99:162,164:174,188:189)]
  
  
  dth_time<-as.numeric(str_replace_all(covariates_red$DTHTIME,":","\\."))
  dthszn<-model.matrix(~0+covariates_red$DTHSEASON)
  colnames(dthszn)<-paste0("DeathSeason_",c("NA","Fall","Spring","Summer","Winter"))
  covariates_red<-cbind(covariates_red,dthszn)
  dthplc<-model.matrix(~0+as.factor(covariates_red$DTHPLCE))
  colnames(dthplc)<-paste0("DeathPlace_",c("NA",levels(as.factor(covariates_red$DTHPLCE))[-1]))
  covariates_red<-cbind(covariates_red,dthplc)
  dthmnnr<-model.matrix(~0+as.factor(covariates_red$DTHMNNR))
  colnames(dthmnnr)<-paste0("DeathManner_",c("NA",levels(as.factor(covariates_red$DTHMNNR))[-1]))
  covariates_red<-cbind(covariates_red,dthmnnr)
  race<-model.matrix(~0+(as.factor(covariates_red$RACE)))
  colnames(race)<-paste0("Race_",c(levels(as.factor(covariates_red$RACE))))
  covariates_red<-cbind(covariates_red,race)
  cohort<-model.matrix(~0+(as.factor(covariates_red$COHORT)))
  colnames(cohort)<-paste0("Cohort_",c("NA",levels(as.factor(covariates_red$COHORT))[-1]))
  covariates_red<-cbind(covariates_red,cohort)
  hardyDeath<-model.matrix(~0+factor(as.character(covariates_red$DTHHRDY),exclude=NULL))
  colnames(hardyDeath)<-paste0("hardyDeath_",c(levels(as.factor((covariates_red$DTHHRDY))),"NA"))
  covariates_red<-cbind(covariates_red,hardyDeath)
  
  
  covariates_red$dth_time_10_14<-as.numeric(dth_time >= 10 & dth_time <14)
  covariates_red$dth_time_14_18<-as.numeric(dth_time >= 14 & dth_time <18)
  covariates_red$dth_time_18_22<-as.numeric(dth_time >= 18 & dth_time < 22)
  covariates_red$dth_time_22_2<-as.numeric(dth_time >= 22 | dth_time <2)
  covariates_red$dth_time_2_6<-as.numeric(dth_time >= 2 & dth_time <6)
  covariates_red$dth_time_6_10<-as.numeric(dth_time >= 6 & dth_time <10)
  
  # covariates_red<-covariates_red%>% select(-COHORT,-DTHTIME,-DTHPLCE-DTHMNNR,-RACE,DTHHRDY)
  
  
  return(covariates_red)
}
covariates_red<-reduced_covars(md_red)

sex_alpha0<-all_preds %>% dplyr::filter(alpha==0)
sex_alpha0_reorder=sex_alpha0[match(covariates_red$SUBJID, sex_alpha0$Ind),]
sex_alpha0.5<-all_preds %>% dplyr::filter(alpha==0.5)
sex_alpha0.5_reorder=sex_alpha0.5[match(covariates_red$SUBJID, sex_alpha0.5$Ind),]
sex_alpha1<-all_preds %>% dplyr::filter(alpha==1)
sex_alpha1_reorder=sex_alpha1[match(covariates_red$SUBJID, sex_alpha1$Ind),]


covars_numeric<-apply(t(covariates_red),1,as.numeric)
corr_sex_md<-cor(covars_numeric[,"SEX"],covars_numeric)
corr_sex_md_alpha0<-cor(sex_alpha0_reorder$sex,covars_numeric)
corr_sex_md_alpha0.5<-cor(sex_alpha0.5_reorder$sex,covars_numeric)
corr_sex_md_alpha1<-cor(sex_alpha1_reorder$sex,covars_numeric)
all_cors<-rbind("sex_binary"=corr_sex_md,
                "sex_alpha0"=corr_sex_md_alpha0,
                "sex_alpha0.5"=corr_sex_md_alpha0.5,
                "sex_alpha1"=corr_sex_md_alpha1)
rownames(all_cors)<-c("sex_binary","sex_alpha0","sex_alpha0.5","sex_alpha1")
all_cors_noNA <- all_cors[,colSums(is.na(all_cors))<nrow(all_cors)]
melt_cors<-melt(all_cors_noNA,na.rm=T)
ggplot(melt_cors %>% dplyr::filter(abs(value)>0.1),aes(Var1,Var2,fill=value))+
  geom_tile()+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                    midpoint = 0, limit = c(-1,1), space = "Lab", 
                    name="Pearson\nCorrelation") +
  geom_text(aes(Var1, Var2, label = round(value,3)), color = "black", size = 4) +
  theme_minimal()
  




###OUTLIERS

library("ggplot2")
library(tidyverse)
library("forcats")
#file<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/relative_risk_aut_both.regress.RData"

# mydir="/oak/stanford/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR"
mydir="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8/RR_AllOutliers"

#groups=c("m","f","both_half.regress") #, "both_half","both_half.sex","both_half.regress")
groups=c("m","f","both") #, "both_half","both_half.sex","both_half.regress
my_cat=c("xci","par","strata","")
#my_cat=c(".xci",".par",".strata")
my_cat=""
# chr_types=c(as.character(c(1:5,7:11,15:21)),"x") #,
# chr_types=c(as.character(1:22),"x") #,
# chr_types=c(as.character(1:22),"x") #,
chr_types=c("x","7","AllAut") #,"AllAut","7") #,
#chr_types=c("x") #,
# chr_types=c("x") #aut
nphen=c(3)
z=c(2,2.5,3)
CADD<-c(0,15)
risk_cat<-c("relative_risk") #,"absolute_risk")
maxmaf<-c("0.1","0.05","0.01","0.001","0.0001") #,"0.001","0.0001")
maxmaf<-c("0.01","0.001","0.0001") #,"0.001","0.0001")
filter_version=c("typesSeenTwice","typesALL","typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                 "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                 "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesBlacklistRemovedALL","typesBlacklistRemovedSeenTwice",
                 "typesGQ10BlacklistRemovedALL","typesGQ5BlacklistRemovedALL",
                 "typesGQ10BlacklistRemovedSeenTwice",  "typesGQ5BlacklistRemovedSeenTwice")
filter_version=c("typesGQ5BlacklistRemovedALL")
filter_version=c("typesALL")
outlier_types=c("outliersTOP","outliers")

regresstype="incl_sex"
regresstype="protect_sex"
null_or_not="preprocessing_v8"
chrgroups=c("x_or_aut")
groups=c("m","f","both_half")
bin_or_cont=c("binary","continuous")
chr_types=c("x","7","aut") #,"AllAut","7") #,
all_risks<-data.frame(matrix(nrow=0,ncol=10))
colnames(all_risks)<-c("Risk","Lower","Upper" ,"Pval","Cat","Type" , "z","nphen","chr","sex")

# #relative_risk_min0max0.0001_x_outliersTOP_cadd15_noglobal_medz-zthresh3-nphen3-m-protect_sex-preprocessing_v8-x-continuous_alpha0.5.txt.gz
for(this_isnull in null_or_not){
  for(this_outlier in outlier_types){
    for(this_maxmaf in maxmaf){
      for (this_sex in groups){
        for (this_chrgroup in chrgroups){
          for(this_chr in chr_types){
          for(this_nphen in nphen){
            for(this_z in z){
              for(this_regresstype in regresstype){
                for(cadd_min in CADD){
                  for(this_bin_cont in bin_or_cont){
                      if(this_chr=="x" & this_chrgroup=="x_or_aut"){
                        use_chrgroup="x"
                      }
                    else if(this_chr!="x" & this_chrgroup=="x_or_aut"){
                      use_chrgroup="aut"
                    }
                    else{
                      use_chrgroup="both"
                      }
                    # file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                    #             "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                    #             this_isnull,"-",use_chrgroup,"-continuous_alpha0.5.txt.gz")
                    # file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                    #             "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                    #             this_isnull,"-",use_chrgroup,"-binary.txt.gz")
                    # if(this_chr=="x"){
                      if(this_bin_cont=="continuous"){
                        # file=paste0(mydir,"/relative_risk_x_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                      #             "_x_",this_sex,"_min0max",this_maxmaf,"_CADDtypesGQ5BlacklistRemovedALL_CADD",
                      #             cadd_min,"_linc_prot.csv")
                      file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                                  "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                                  this_isnull,"-",use_chrgroup,"-continuous_alpha0.5.txt.gz")
                      risks$discrete_or_continuous<-"continuous"
                      
                                  #_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                                  #"-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                                  #this_isnull,"-",use_chrgroup,"-binary.txt.gz")
                    }else{
                      # file=paste0(mydir,"/relative_risk_",this_outlier,"_z",this_z,"_nphen",this_nphen,
                      #             "_",this_sex,"_",this_chr,"_min0max",this_maxmaf,"_CADDtypesALL_CADD",cadd_min,
                      #             "_linc_prot.txt")
  
                      file=paste0(mydir,"/relative_risk_min0max",this_maxmaf,"_",this_chr,"_",this_outlier,"_cadd",cadd_min,"_noglobal_medz-zthresh",this_z,
                                  "-nphen",this_nphen,"-",this_sex,"-",this_regresstype,"-",
                                  this_isnull,"-",use_chrgroup,"-binary.txt.gz")
                      risks$discrete_or_continuous<-"discrete"
                      
                      }
  
                     risks=read.csv(file[1])
                    #load(file)
                    risks$z<-this_z
                    risks$nphen<-this_nphen
                    risks$chr<-this_chr
                    if(this_sex=="both_half"){use_sex="both"}else{use_sex=this_sex}
                    risks$sex<-use_sex
                    risks$maxmaf<-this_maxmaf
                    # risks$filt<-this_filt
                    risks$cadd<-cadd_min
                    # risks$thiscat<-this_cat
                    risks$outlierType<-this_outlier
                    colnames(risks)[6]<-"CATEGORY"
                    risks$real_or_null<-"null"
                    #assign(varname,risks)
                    #print(varname)
                    # break
  
                    all_risks<-rbind(all_risks,risks)
}}}}}}}}}}}

all_risks$outliers_tested<-all_risks$exp_yn+all_risks$exp_yy
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/ following this
risks_to_compare<-all_risks %>% group_by(outlierType,CATEGORY,z,nphen,cadd,chr,sex,maxmaf,real_or_null,exp_type) %>% 
  dplyr::summarise(across(c(Risk,Lower,Upper),log)) %>%
  mutate(SE=abs(Lower-Upper)/(2*1.96))%>%dplyr::filter(sex=="m" | sex=="f")
risks_to_compare_side<-pivot_wider(risks_to_compare,names_from=c(sex),
                                   id_cols=c(outlierType,CATEGORY,z,nphen,cadd,chr,maxmaf,real_or_null,exp_type) ,
                                   values_from=c(Risk,SE)) #id_cols=, names_from=,
risks_to_compare_side$RiskDiff<-risks_to_compare_side$Risk_f-risks_to_compare_side$Risk_m
risks_to_compare_side$SEDiff<-sqrt(risks_to_compare_side$SE_f**2+risks_to_compare_side$SE_m**2)
risks_to_compare_side$zDiff<-(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
risks_to_compare_side$p<-pnorm(risks_to_compare_side$RiskDiff/risks_to_compare_side$SEDiff)
#n is 2 for cadd * 3 for z * 3 for nphen * 3 for MAF * 24 for chr  (23 + allAut)=1296 OR 54* 2 for outliers
risks_to_compare_side$padj<-p.adjust(risks_to_compare_side$p,method = "BH") #,n=108)
all_risks_p<-merge(all_risks,risks_to_compare_side[,c("outlierType","exp_type","CATEGORY","z","nphen","cadd","chr","maxmaf","p","padj")])
