#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)
library(optparse)


option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--out_rdata_relative"), type = 'character', default = NULL, help = "path of output file (RDATA) relative risk"),
                make_option(c("--zscore"), type = 'numeric', default = NULL, help = "min z score"),
                make_option(c("--nphen"), type = 'numeric', default = NULL, help = "min tissues"),
                make_option(c("--gtf_code_file"), type = 'character', default = NULL, help = "gtf_code_file preprocessing"),
                make_option(c("--this_chr"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--this_method"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--sex"), type = 'character', default = NULL, help = "ind sex"),
                make_option(c("--min_maf"), type = 'numeric', default = NULL, help = "min MAF"),
                make_option(c("--max_maf"), type = 'numeric', default = NULL, help = "max MAF"),
                make_option(c("--cadd_min"), type = 'numeric', default = NULL, help = "min cadd score")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
out_rdata_relative <- as.character(opt$out_rdata_relative)
my_zscore <- as.numeric(opt$zscore)
cadd_min <- as.numeric(opt$cadd_min)
nphen <- as.numeric(opt$nphen)
sex <- as.character(opt$sex)
gtf_code_file<-as.character(opt$gtf_code_file)
this_chr<-as.character(opt$this_chr)
this_method<-as.character(opt$this_method)
min_maf<-as.numeric(opt$min_maf)
max_maf<-as.numeric(opt$max_maf)
# #
# infile<-"/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Outliers/outliers_noglobal_varAnnot_medz_zthresh2_nphen3_x_m_beta0.111_cadd15_x.txt.gz"
# infile="/Volumes/groups/smontgom/raungar/Sex/Output/sexdeg_v8/Outliers/outliersTOP_noglobal_varAnnot_medz_zthresh2.5_nphen3_x_f_beta0.04954_cadd15_x.txt.gz"
# # infile="/Volumes/groups/smontgom/raungar/Sex/Output/enrichments_v8/OutliersAndRVs/all_outliers_noglobal_medz_varAnnot_zthresh3_nphen5_f_CADDtypesALL_linc_prot.txt.gz"
# gtf_code_file<-"/Volumes/groups/smontgom/raungar/Sex/Output/preprocessing_v8/both_proteincoding_lncrna.gtf"
# my_zscore<-2.5
# min_maf<-0
# nphen=3
 # max_maf<-0.001
#this_method="outliers"
#sex="f"
 # cadd_min=15

# if(file.exists(out_rdata_relative)){stop("outfile - relative exists")}
#if(file.exists(out_rdata_continuous)){stop("outfile - continuous exists")}

### expression
print('Reading expression')
print(paste0("reading in: ",infile))
exp_data = fread(infile,data.table=F)
#exp_data = exp_data %>% select(indiv_id,gene_id,MedZ,Y,tier2,af_gtex,categoryOutlier,genetype,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
#exp_data = exp_data %>% dplyr::select(ind,ensg,MedZ,Y,af,isOutlier,genetype,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
#exp_data = dplyr::filter(exp_data,sv_v7==1)
#new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$genetype[x] == 'splice', exp_data$tier2[x], exp_data$genetype[x]))
#exp_data$genetype = ind	ensg	N	Df	MedZ	Y	chr	pos	vartype	sex	ref	alt	varswitch	gtex_maf	gnomad_maf_both	gnomad_maf_m	gnomad_maf_f	genetype	gnomad_maf_diff	num_rvs
colnames(exp_data)<-c("ind","ensg","N","Df","MedZ","Y","AbsBetaSex","MaxOutlierBetaSex","MinOutlierBetaSex","AbsOutlierBetaSex",
                      "chr","chrNum","start","end","vartype","sex",
                     "gtex_maf","gnomad_maf","use_maf","genetype","cadd_raw","cadd_phred","numrv")

exp_data$cadd_phred<-as.numeric(exp_data$cadd_phred)
exp_data$cadd_phred[is.na(exp_data$cadd_phred)] <- 0
# exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

print("filter for MAF")
#choose rare/common
has_variant<-apply((exp_data),1,function(x){
  this_maf<-as.numeric(x[which(colnames(exp_data)=="use_maf")])
  this_cadd_phred<-as.numeric(x[which(colnames(exp_data)=="cadd_phred")])
  #print(paste0("MAF=",this_maf,", cadd=",this_cadd_phred))
  #no variants found w/in 10kb of gene, so not rare
  if(is.na(this_maf)){"common"}
  else if (this_maf>=min_maf & this_maf<max_maf){
    if(this_cadd_phred>=cadd_min){"rare"}
    else if(this_cadd_phred<cadd_min){"common"}
    else (stop("ERROR: VARIANT NOT MAKING SENSE (at cadd level)"))
  }
  else if (this_maf>=max_maf){"common"}
  else if (this_maf<min_maf){"common"}
  else (stop("ERROR: VARIANT NOT MAKING SENSE (at maf level)"))
})
exp_data$has_variant<-has_variant



gtf_code<-fread(gtf_code_file,header=F)
genetype_dic=gtf_code$V2
names(genetype_dic)<-gtf_code$V1
exp_data$genetype<-genetype_dic[exp_data$ensg] 

print("filter for z score")

print(head(exp_data))
# exp_outliers = exp_data %>% dplyr::filter(abs(MedZ) >= my_zscore)  %>% dplyr::filter(as.numeric(Df)>=nphen)
# exp_controls = dplyr::filter(exp_data, (abs(MedZ) < my_zscore) | (abs(MedZ) >= my_zscore & as.numeric(Df)<nphen))
print("AND")
colnames(exp_data)[24]<-"some_col"
print(colnames(exp_data))
print(exp_data%>% pull(AbsBetaSex) %>% head())
exp_data_anysexdeg<-exp_data%>% dplyr::filter(as.numeric(AbsBetaSex)>0) #AbsBetaSex #was AbsOutlierBetaSex
print("0")
exp_outliers = exp_data_anysexdeg %>% dplyr::filter(Y=="outlier") #was just expdata
print("1")
exp_outliers_over = exp_data_anysexdeg %>% dplyr::filter(Y=="outlier") %>% dplyr::filter(MedZ>0)#was just expdata
print("2")
exp_outliers_under = exp_data_anysexdeg %>% dplyr::filter(Y=="outlier") %>% dplyr::filter(MedZ<0)#was just expdata
print("3")
exp_controls = exp_data_anysexdeg %>% dplyr::filter(Y=="control") #was jsut expdata

print("outliers have been split as have controls")

print("relative risk")
### Relative risk
# risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
# risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Cat = character(),
#                    exp_nn=numeric(),exp_ny=numeric(),exp_yn=numeric(),exp_yy=numeric(),
#                   sex=character(),z=numeric(),nphen=numeric(),Type=character())
risks=data.frame()
vcats =as.character(na.omit(unique(exp_data$genetype)))
sexdegs=c("MaxOutlierBetaSex","MinOutlierBetaSex","AbsOutlierBetaSex")
sexdeg_dic<-c("female","male","any")
names(sexdeg_dic)<-sexdegs
exp_types=c("all","over","under")
print(risks)
for(exp_type in exp_types){
  if(exp_type=="all"){
    this_exp_outliers<-exp_outliers
      
  }else if(exp_type=="over"){
    this_exp_outliers<-exp_outliers_over
      
  }else if(exp_type=="under"){
    this_exp_outliers<-exp_outliers_under
      
  }else{stop("ERROR: INVALID EXP TYPE")}
  for(sexdeg in sexdegs){
    print(sexdeg)
    for (vcat in (vcats)) {
      print(vcat)
      print(cadd_min)
       #any
       if(sexdeg=="AbsOutlierBetaSex"){
         outliers_passed = this_exp_outliers %>% dplyr::filter(genetype == vcat)%>% dplyr::filter(AbsOutlierBetaSex > 0)
         exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         # exp_yy = nrow(this_exp_outliers %>% dplyr::filter(genetype == vcat) %>% dplyr::filter( AbsOutlierBetaSex > 0) %>% dplyr::filter(has_variant == "rare")) # %>% dplyr::filter(cadd_phred>cadd_min))
         
      } else if(sexdeg=="MaxOutlierBetaSex"){
         ##female
        outliers_passed = this_exp_outliers %>% dplyr::filter(genetype == vcat)%>% dplyr::filter(MaxOutlierBetaSex > 0)
        exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
        
       }else if(sexdeg=="MinOutlierBetaSex"){
         outliers_passed = this_exp_outliers %>% dplyr::filter(genetype == vcat)%>% dplyr::filter(MinOutlierBetaSex < 0)
         exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg) %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
         
         
       }else(stop("ERROR, IMPROPER SEX DEG"))
  	print(paste0("the following: ",exp_nn," and ",exp_ny," and ",exp_yn," and ",exp_yy))
     #exp_nn = nrow(this_exp_outliers %>% dplyr::filter(genetype == vcat) %>% dplyr::filter( has_variant == "rare"))
     #exp_ny = nrow(this_exp_outliers %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant != "rare"))
     #exp_yn = nrow(exp_controls %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant== "rare"))
     #exp_yy = nrow(exp_controls %>% dplyr::filter(genetype == vcat) %>% dplyr::filter(has_variant != "rare"))
     exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
     print(exptable)
     riskvals=tryCatch(
      {
        err=epitab(exptable, method = 'riskratio') 
        myrisk = err$tab[2,5]
         this_lower = err$tab[2,6]
        this_upper = err$tab[2,7]
        this_p = err$tab[2,8]
        list(myrisk,this_lower,this_upper,this_p)
        },error=function(cond) {
          message(cond)
          myrisk=NA
          this_lower=NA
          this_upper=NA
          this_p=NA
          list(myrisk,this_lower,this_upper,this_p)
        }
     )
     this_risk=data.frame(Risk = unlist(riskvals)[1],
                          Lower = unlist(riskvals)[2],
                          Upper = unlist(riskvals)[3],
                          Pval = unlist(riskvals)[4],
                          Cat = vcat,
                          exp_nn=exp_nn,
                          exp_ny=exp_ny,
                          exp_yn=exp_yn,
                          exp_yy=exp_yy,max_maf=max_maf,
                          chr=this_chr,outlierType=this_method,
                          num_outliers=nrow(this_exp_outliers%>% dplyr::filter(genetype == vcat)),
                          sex=sex,z=my_zscore,nphen=nphen,cadd=cadd_min,sexdeg_type=sexdeg_dic[sexdeg],
                          outlier_direction=exp_type,Type = 'SexDEG',row.names=NULL)
     rm(myrisk);rm(this_lower);rm(this_upper); rm(this_p)
     print(this_risk)
     risks = rbind(risks,this_risk )   
      }
  }
  #DO NO FILTER VARCAT
  for(sexdeg in sexdegs){
   
     if(sexdeg=="AbsOutlierBetaSex"){
       outliers_passed = this_exp_outliers %>% dplyr::filter(AbsOutlierBetaSex > 0)
       exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)  %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)%>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       
     }
     ##femaleß
     else if(sexdeg=="MaxOutlierBetaSex"){
       outliers_passed = this_exp_outliers %>% dplyr::filter(MaxOutlierBetaSex > 0)
       exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)  %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)%>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       
     }else if(sexdeg=="MinOutlierBetaSex"){
       outliers_passed = this_exp_outliers %>% dplyr::filter(MinOutlierBetaSex < 0)
       exp_yn=nrow( outliers_passed   %>% dplyr::filter(has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_yy=nrow( outliers_passed   %>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_nn = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)  %>% dplyr::filter( has_variant != "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       exp_ny = nrow(exp_controls %>% dplyr::filter(ensg %in% outliers_passed$ensg)%>% dplyr::filter(has_variant == "rare"))  #%>% dplyr::filter(cadd_phred>cadd_min))
       
     }else(stop("ERROR, IMPROPER SEX DEG"))
    #exp_nn_all = nrow(this_exp_outliers %>% dplyr::filter( has_variant == "rare"))
    #exp_ny_all = nrow(this_exp_outliers %>% dplyr::filter(has_variant != "rare"))
    #exp_yn_all = nrow(exp_controls %>% dplyr::filter(has_variant== "rare"))
    #exp_yy_all = nrow(exp_controls  %>% dplyr::filter(has_variant != "rare"))
    exptable_all = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
    colnames(exptable_all)<-c("control","outliers")
    rownames(exptable_all)<-c("common","rare")
    riskvals=tryCatch(
      {
        err=epitab(exptable_all, method = 'riskratio') 
        myrisk = err$tab[2,5]
        this_lower = err$tab[2,6]
        this_upper = err$tab[2,7]
        this_p = err$tab[2,8]
        list(myrisk,this_lower,this_upper,this_p)
      },error=function(cond) {
        message(cond)
        myrisk=NA
        this_lower=NA
        this_upper=NA
        this_p=NA
        list(myrisk,this_lower,this_upper,this_p)
      }
    )
    this_risk=data.frame(Risk = unlist(riskvals)[1],
                         Lower = unlist(riskvals)[2],
                         Upper = unlist(riskvals)[3],
                         Pval = unlist(riskvals)[4],
                         Cat = "all",
                         exp_nn=exp_nn,
                         exp_ny=exp_ny,
                         exp_yn=exp_yn,
                         exp_yy=exp_yy,max_maf=max_maf,
                         chr=this_chr,outlierType=this_method,
                         num_outliers=nrow(this_exp_outliers),
                         sex=sex,z=my_zscore,nphen=nphen,cadd=cadd_min,sexdeg_type=sexdeg_dic[sexdeg],
                         outlier_direction=exp_type,Type = 'SexDEG',row.names=NULL)
    rm(myrisk);rm(this_lower);rm(this_upper); rm(this_p)
    
    risks = rbind(risks, this_risk)   
  }
  risks = risks %>% arrange(by=Risk) 
}
risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))
print(risks)
#save(risks,  file=out_rdata_relative)


write.csv(risks,  file=out_rdata_relative,quote=F,sep="\t")
print(paste0("COMPLETE: ",out_rdata_relative))
