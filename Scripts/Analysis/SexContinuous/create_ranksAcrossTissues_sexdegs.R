library(data.table)
library(dplyr)
library("optparse") #for passing args
library("crunch") #for write.csv.gz

get_opt_parser<-function(){
        option_list = list( make_option(c("-f", "--infile"), type="character", default=NULL,help="infile (sexdeg file)"),
			make_option(c("-o", "--outfile"), type="character", default=NULL,help="outfile")
	)
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
infile=args$infile
outfile=args$outfile


infile="/oak/stanford/groups/smontgom/raungar/Sex/Output/analysis_v8/continuous/Ranks/aut_beta0.111.all.txt"
all_inds_tiss<-fread(infile)
all_inds_tiss_meds<-all_inds_tiss[,c("ind","Median","tissue","sex")]

head(all_inds_tiss_meds)

#tmp<-sexdegs[sample(1:nrow(sexdegs),100,replace=F),]
inds_summary<-all_inds_tiss_meds[, as.list(summary(Median)[1:6]), by = c("ind","sex")]

write.csv(inds_summary,file=outfile,sep="\t",quote=F,  col.names=FALSE)
