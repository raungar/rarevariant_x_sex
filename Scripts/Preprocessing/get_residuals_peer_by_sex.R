library("optparse") #for passing args
library("data.table") #for fread
set.seed(1234) #important: randomly subsetting individuals

#get arguments
get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL, help="PEER directory"),
		 make_option(c("-s", "--subset_do"), type="character", default="y", help="subset files: y/yes, don't do subsetting: n/no"),
		 make_option(c("-r", "--residuals_do"), type="character", default="y", help="compute residuals: y/yes, do not compute residuals: n/no"),
		 make_option(c("-m", "--metadata_file"), type="character", default=NULL, help="metadata file that contains sex")
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}

#get which individuals to subset via sex for particular tissue
get_individuals<-function(tissue_path,subset_do,residuals_do,md_dic,used_inds){
	#this is ugly, but just messing with the string to get the actual tissue name
	tissue_name<-gsub('.{1}$', '',sapply(strsplit(sapply(strsplit(tissue_path,"/"), tail,1),"Factor"),head,1))
	#read in factors
	factor_file<-paste0(tissue_path,"/factors.tsv")
	print(factor_file)
	factors<-as.data.frame(fread(factor_file,sep="\t",header=T))

	inds_all<-colnames(factors)[-1] #colnames are the individuals
	md_dic_red<-md_dic[names(md_dic) %in% inds_all] #make sure the dictionary only contains individuals that this tissue has


	####GET THE NAMES OF THE INDIVIDUALS WE ARE SUBSETTING TO
	sex_table<-table(md_dic_red[inds_all])
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
		inds_m<-md_dic_red[md_dic_red==min_sex]
		inds_f<-sample(md_dic_red[md_dic_red==max_sex])[1:sex_table[min_sex]]
	}else if(min_sex==2){
		min_sex_fm="female"
		inds_f<-md_dic_red[md_dic_red==min_sex]
		inds_m<-sample(md_dic_red[md_dic_red==max_sex])[1:sex_table[min_sex]]
	}else {
		min_sex_fm=stop("ERROR")
	}
	inds_b<-c(inds_m,inds_f)

	print(paste0(tissue_name,": subsetting to ", sex_table[min_sex], " individuals (",min_sex_fm," have fewest)"))

	#add to list where key is tissue name, value is another list with the individual ids		
	used_inds[[tissue_name]]<-list(male=names(inds_m),female=names(inds_f),both=names(inds_b))
	###used_inds_new<-unique(used_inds,inds_b)

	#substted factors
	factors_m<-factors[,c("ID",names(inds_m))]
	factors_f<-factors[,c("ID",names(inds_f))]
	factors_b<-factors[,c("ID",names(inds_b))]


	###WRITE THESE SUBSETTED FILES
	subset_dir_path<-paste0(tissue_path,"/Subset")#if does not exist, create new path for subsetted files
	if(!dir.exists(subset_dir_path)){dir.create(subset_dir_path)}
	if(subset_do == "y" || subset_do == "yes"){
		write.table(factors_m,file=paste0(subset_dir_path,"/factors_m.tsv"),sep="\t",quote=F)
		write.table(factors_f,file=paste0(subset_dir_path,"/factors_f.tsv"),sep="\t",quote=F)
		write.table(factors_b,file=paste0(subset_dir_path,"/factors_b.tsv"),sep="\t",quote=F)
	}

	###CALCULATE RESIDUALS
	resid_dir_path<-paste0(tissue_path,"/ResidualsSex")#if does not exist, create new path for subsetted files
        if(!dir.exists(resid_dir_path)){dir.create(resid_dir_path)}
	if(residuals_do == "y" || residuals_do == "yes"){
		#independent variable sex
		ind_sex<-as.matrix(as.factor(md_dic_red[colnames(factors_b)[-1]]))
		#set dependent variable
		factors_b_forlm<-as.matrix(factors_b[,-1])
		rownames(factors_b_forlm)<-factors_b[,1]

		#combine these variables into one dataframe for lm()
		lm_matrix<-cbind.data.frame(t(factors_b_forlm),ind_sex)
		#formula where it is essentially factors ~ sex, specifically cbind(Factor1, Factor2, ..., FactorN) ~ Sex
		lm_form<-as.formula(paste0("cbind(",paste0(rownames(factors_b_forlm),collapse=","),") ~ ind_sex"))
		#calculate the lm fit and residuals
		lm_fit<-lm(lm_form,data=as.data.frame(lm_matrix))
		lm_resid<-residuals(lm_fit)

		#write this to the new folder
                write.table(lm_resid,file=paste0(resid_dir_path,"/residuals_sex.tsv"),sep="\t",quote=F)

	}


	return(used_inds)	
	
}

#checks to make sure parameters were properly passed into file
check_params<-function(peer_dir,subset_do,residuals_do,metadata_file){
	if(length(peer_dir) == 0){stop("ERROR: peer file dir not provided")}
	if(substr(peer_dir,nchar(peer_dir),nchar(peer_dir))=="/"){stop("Remove trailing / from peer_dir for proper results")}
	if(subset_do != "y" & subset_do != "yes" & subset_do != "n" & subset_do != "no"){
		stop("ERROR: subset_do must be the following: y, yes, n, no")
	}
	if(residuals_do != "y" & residuals_do != "yes" & residuals_do != "n" & residuals_do != "no"){
		stop("ERROR: residuals_do must be the following: y, yes, n, no")
	}
	if(length(metadata_file) == 0){stop("ERROR: metadata file not provided")}
}

####MAIN

opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
peer_dir<-as.character(args$dir_peer)
subset_do<-tolower(as.character(args$subset_do))
residuals_do<-tolower(as.character(args$residuals_do))
metadata_file<-as.character(args$metadata_file)

#check for parameters being passed properly, if not quit
check_params(peer_dir,subset_do,residuals_do,metadata_file)

#make a dictionary where keys are individuals, vales are sex
metadata<-fread(metadata_file,sep="\t",header=T)
ind_dict<-metadata$SEX
names(ind_dict)<-metadata$SUBJID


used_inds<-list() #individuals who have been included in other tissues
for(tissue_dir in list.dirs(path=as.character(peer_dir))){
	if(tissue_dir == peer_dir){next} #for some reason list.dirs prints current dir, so skip that
	#used_inds are individuals that have been subsetted, this loops so at the end it is a master list
	#used_inds is called by used_inds$tissue_name$female where female is male/female/both
	print(tissue_dir)
	used_inds<-get_individuals(tissue_dir,subset_do,residuals_do,ind_dict,used_inds)
	break
}

print("EOF")
