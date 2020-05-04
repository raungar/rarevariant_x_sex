library("optparse") #for passing args
library("data.table") #for fread
set.seed(1234) #important: randomly subsetting individuals

#get arguments
get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL),
		 make_option(c("-s", "--sex_option"), type="character", default=NULL),
		 make_option(c("-m", "--metadata_file"), type="character", default=NULL)
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}

#get which individuals to subset via sex for particular tissue
get_individuals<-function(tissue,sex_option,md_dic,used_inds){
	#this is ugly, but just messing with the string to get the actual tissue name
	tissue_name<-gsub('.{1}$', '',sapply(strsplit(sapply(strsplit(tissue,"/"), tail,1),"Factor"),head,1))
	#read in factors
	factor_file<-paste0(tissue,"/factors.tsv")
	factors<-fread(factor_file,sep="\t",header=T)	
	inds_all<-colnames(factors)[-1] #colnames are the individuals
	sex_table<-table(md_dic[inds_all])

	#IF SEX-SPECIFIC TISSUE DON'T DO ANYTHING
	if(length(sex_table)==1){
		print(paste0(tissue_name,": no subsetting, this is a sex-specific tissue"))
		return(used_inds)	
	}

	#get min and max sex and respond appropriately
	min_sex<-names(sex_table)[sex_table %in% min(sex_table)]
	max_sex<-names(sex_table)[sex_table %in% max(sex_table)]
	if(min_sex==1){
		mix_sex_fm="male"
		inds_m<-md_dic[md_dic==min_sex]
		inds_f<-sample(md_dic[md_dic==max_sex])[1:sex_table[min_sex]]
	}else if(min_sex==2){
		min_sex_fm="female"
		inds_f<-md_dic[md_dic==min_sex]
		inds_m<-sample(md_dic[md_dic==max_sex])[1:sex_table[min_sex]]
	}else {
		min_sex_fm=stop("ERROR")
	}
	print(paste0(tissue_name,": subsetting to ", sex_table[min_sex], " individuals (",min_sex_fm," have fewest)"))
		
	inds_b<-c(inds_m,inds_f)
	used_inds[[tissue_name]]<-list(male=inds_m,female=inds_f,both=inds_b)
	used_inds_new<-unique(used_inds,inds_b)
	#inds_min<-

	#First, make a dictionary where values are individual ids, values are sex
	to_extract<-c()
	if(sex_option == "m" || sex_option == "male"){
		to_extract<-c(to_extract,"M")
	} else if (sex_option == "f" || sex_option == "female"){
		to_extract<-c(to_extract,"F")
	} else if (sex_option == "b" || sex_option == "both"){
		to_extract<-c(to_extract,"B")
	} else if (sex_option == "a" || sex_option == "all"){
		to_extract<-c(to_extract,"M","F","B")
	} else {
		#This shouldn't happen, there was an error correct option not provided
		return(-1)
	}


	return(used_inds)	
	
}

#checks to make sure parameters were properly passed into file
check_params<-function(peer_dir,sex_option,metadata_file){
	if(peer_dir == "FALSE"){stop("ERROR: peer file dir not provided")}
	if(substr(peer_dir,nchar(peer_dir),nchar(peer_dir))=="/"){stop("Remove trailing / from peer_dir for proper results")}
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

#make a dictionary where keys are individuals, vales are sex
metadata<-fread(metadata_file,sep="\t",header=T)
ind_dict<-metadata$SEX
names(ind_dict)<-metadata$SUBJID


used_inds<-list() #individuals who have been included in other tissues
for(tissue_dir in list.dirs(path=as.character(peer_dir))){
	if(tissue_dir == peer_dir){next} #for some reason list.dirs prints current dir, so skip that
	used_inds<-get_individuals(tissue_dir,sex_option,ind_dict,used_inds)
	#break
}

print("EOF")
