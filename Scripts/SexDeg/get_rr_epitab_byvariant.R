library(data.table)
library(epitools)
library(optparse)

option_list = list(
                make_option(c("--infile"), type = 'character', default = NULL, help = "path of input file"),
                make_option(c("--outfile"), type = 'character', default = NULL, help = "path of out rr file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

infile <- as.character(opt$infile)
outfile <- as.character(opt$outfile)
#infile="tmp"

rr_data=fread(infile,data.table=F)
# num_rvs_nonsexdegs num_nonsexdegs num_rvs_sexdegs num_sexdegs
exptable<-rbind(
	 	c(rr_data$num_rvs_sexdegs,rr_data$num_sexdegs),
		c(rr_data$num_rvs_nonsexdegs,rr_data$num_nonsexdegs)
	)
#print(exptable)

err = epitab(exptable, method = 'riskratio')
risks = data.frame(Risk = err$tab[2,5],
                                   Lower = err$tab[2,6],
                                   Upper = err$tab[2,7],
                                   Pval = err$tab[2,8],
                                   Type = 'Total expression')
#print(risks)
write.csv(risks,file=outfile)
