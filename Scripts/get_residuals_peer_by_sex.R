library("optparse") #for passing args


get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL),
		 make_option(c("-s", "--sex_option"), type="character", default=NULL)
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}


####MAIN

opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
peer_dir<-as.character(args[1])
sex_option<-tolower(as.character(args[2]))

print(sex_option)
#check for file being passed, if not quit
if(peer_dir == "FALSE"){stop("ERROR: peer file dir not provided")}
if(sex_option != "m" & sex_option != "f" & sex_option != "b" & sex_option != "male" & sex_option != "female" & sex_option != "both" ){
	stop("ERROR: sex_option must be the following: m, male, f, female, b, both")
}

print(peer_dir)

for(this_dir in list.dirs(path=as.character(peer_dir))){
	print(this_dir)
	print("DIR ^")
	break
}

print("EOF")
