
get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL),
		 make_option(c("-s", "--sex_option"), type="character", default=NULL),
		 make_option(c("-m", "--metadata_file"), type="character", default=NULL)
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}


#checks to make sure parameters were properly passed into file
check_params<-function(peer_dir,sex_option,metadata_file){
	if(peer_dir == "FALSE"){stop("ERROR: peer file dir not provided")}
	if(sex_option != "m" & sex_option != "f" & sex_option != "b" & sex_option != "male" & sex_option != "female" & sex_option != "both"  & sex_option != "a" & sex_option != "all"){
		stop("ERROR: sex_option must be the following: m, male, f, female, b, both, a, all")
	}
	if(metadata_file == "FALSE"){stop("ERROR: metadata file not provided")}
}

####MAIN

opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
peer_dir<-as.character(args[1])
sex_option<-tolower(as.character(args[2]))
metadata_file<-as.character(args[3])

#check for parameters being passed properly, if not quit
check_params(peer_dir,sex_option,metadata_file)

print(peer_dir)

for(this_dir in list.dirs(path=as.character(peer_dir))){
	print(this_dir)
	print("DIR ^")
	break
}

print("EOF")
