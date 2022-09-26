library(data.table)
library(dplyr)
library("optparse") #for passing args
library("crunch") #for write.csv.gz

get_opt_parser<-function(){
        option_list = list( make_option(c("-f", "--infile"), type="character", default=NULL,help="infile (sexdeg file)"),
			make_option(c("-t", "--tissue"), type="character", default=NULL,help="this tissue"),
			make_option(c("-o", "--outfile"), type="character", default=NULL,help="this tissue"),
			make_option(c("-s", "--sex"), type="character", default=NULL,help="sex")
	)
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}
opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
infile=args$infile
tissue=args$tissue
sex=args$sex
outfile=args$outfile


#infile="/oak/stanford/groups/smontgom/raungar/Sex/Output/sexdeg_v8/CombinedSingleTissue/aut_beta0.111_NERVET_m_linc_prot.txt.gz"
#tissue=""


all_genes<-fread(infile)
#get genes that are sexdegs, and get cols ind/ensg/z/beta
sexdegs<-all_genes %>% dplyr::filter(beta != 0) %>% select(ind,ensg,z,beta)
#this means that all positive z scores are "more female"
#and all negative z scores are "more male"
sexdegs$adjusted_z<-sign(sexdegs$beta)*sexdegs$z

#tmp<-sexdegs[sample(1:nrow(sexdegs),100,replace=F),]
sexdegs_summary<-sexdegs[, as.list(summary(adjusted_z)[1:6]), by = c("ind")]
sexdegs_summary$tissue<-tissue
sexdegs_summary$sex<-sex

write.csv(sexdegs_summary,file=outfile,sep="\t",quote=F,  col.names=FALSE)
