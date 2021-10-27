library(optparse)
library(data.table)


option_list = list(
  make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--inclsex"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--sex_continuous"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--this_sex"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--is_null"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--chr_subgroup"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--tissue"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--alpha"), type = 'character', default = NULL, help = "path of input file"),
  make_option(c("--outfile"), type = 'character', default = NULL, help = "path of input file")
  
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
infile <- as.character(opt$infile)
regression_method<-as.character(opt$inclsex)
sex_continuous<-as.character(opt$sex_continuous)
null_or_real<-as.character(opt$is_null)
this_sex<-as.character(opt$this_sex)
chr_subgroup<-as.character(opt$chr_subgroup)
alpha<-as.character(opt$alpha)
outfile<-as.character(opt$outfile)
tissue<-as.character(opt$tissue)

# 
# infile="/Volumes/groups/smontgom/raungar/Sex/Output/expression_v8/Counts/counts-Adipose_Subcutaneous-both-incl_sex-preprocessing_v8-x-continuous_alpha0.5.txt"
# this_sex="both_half"
# sex_continuous="T"
# regression_method="incl_sex"
# chr_subgroup="both_half"
# alpha="0.5"
# null_or_real="SEX"


data=fread(infile)
data_melted<-melt.data.table(data,id.vars=c("Id"),variable.name = 'Ind', value.name = 'Counts')
full_info<-(data_melted)[,c("this_sex","sex_cont","regression_method",
                                "chr_subgroup","alpha","null_or_real","tissue"):=
                               list(this_sex,sex_continuous,regression_method,
                                    chr_subgroup,alpha,null_or_real,tissue)]

fwrite(full_info,file=outfile,quote=F,col.names=T,row.names = F,sep="\t")
#outfile is: id/gene/val/tiss/sex/regresstype/null/cont/chr_group/discrete or cont/alpha

