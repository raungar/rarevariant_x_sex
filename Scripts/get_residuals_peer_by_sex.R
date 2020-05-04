library("optparse") #for passing args


get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL)
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}


####MAIN

opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
peer_dir<-args[1]

#check for file being passed, if not quit
if(is.logical(peer_dir$help)){stop("ERROR: peer file dir not provided")}

print(peer_dir)
print("FIN")


for(file in list.files("RAREDIR/preprocessing_v8/PEER_v8/")){
	print(file)
	break
}
