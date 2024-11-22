library("optparse") #for passing args
library("data.table") #for fread
set.seed(1234) #important: randomly subsetting individuals

#get arguments
get_opt_parser<-function(){
	option_list = list(
		 make_option(c("-d", "--dir_peer"), type="character", default=NULL, help="PEER directory"),
		 make_option(c("-f", "--file_to_correct"), type="character", default=NULL, help="current file to process"),
		 #make_option(c("-s", "--subset_do"), type="character", default="y", help="subset files: y/yes, don't do subsetting: n/no"),
		 #make_option(c("-r", "--residuals_do"), type="character", default="y", help="compute residuals: y/yes, do not compute residuals: n/no"),
		 make_option(c("-m", "--metadata_file"), type="character", default=NULL, help="metadata file that contains sex")
		 ) 
	
	opt_parser = OptionParser(option_list=option_list)
	return(opt_parser)
}

#get which individuals to subset via sex for particular tissue
#get_individuals<-function(tissue_path,subset_do,residuals_do,md_dic,used_inds){
get_individuals<-function(tissue_path,md_dic,used_inds){
	#this is ugly, but just messing with the string to get the actual tissue name
	##tissue_name<-gsub('.{1}$', '',sapply(strsplit(sapply(strsplit(tissue_path,"/"), tail,1),"Factor"),head,1))
	tissue_name<-sapply(strsplit(sapply(strsplit(tissue_path,"/"), tail,1),"\\."),"[[",1)
	#read in factors
	print(tissue_name)
	path<-sapply(strsplit(tissue_path,tissue_name),"[[",1)
	factors_b<-as.data.frame(fread(tissue_path,sep="\t",header=T))

	inds_all<-colnames(factors_b)[-1] #colnames are the individuals
	md_dic_red<-md_dic[names(md_dic) %in% inds_all] #make sure the dictionary only contains individuals that this tissue has



	###CALCULATE RESIDUALS
	#resid_dir_path<-paste0(tissue_path,"/ResidualsSex")#if does not exist, create new path for subsetted files
	resid_dir_path<-tissue_path
        #if(!dir.exists(resid_dir_path)){dir.create(resid_dir_path)}
	#$if(residuals_do == "y" || residuals_do == "yes"){
		#independent variable sex
		ind_sex<-as.matrix(as.factor(md_dic_red[colnames(factors_b)[-1]]))
		ind_sex_num<-as.numeric(ind_sex)-1
		#set dependent variable
		factors_b_forlm<-t(as.matrix(factors_b[,-1]))
		colnames(factors_b_forlm)<-factors_b[,1]

		#combine these variables into one dataframe for lm()
		lm_matrix<-cbind.data.frame(factors_b_forlm,ind_sex_num)
		#print(head(lm_matrix))
		#formula where it is essentially factors ~ sex, specifically cbind(Factor1, Factor2, ..., FactorN) ~ Sex
		lm_form<-as.formula(paste0("cbind(",paste0(colnames(factors_b_forlm),collapse=","),") ~ ind_sex"))
		#calculate the lm fit and residuals
		#print(lm_form)
		lm_fit<-lm(lm_form,data=as.data.frame(lm_matrix))
		lm_fit_summary<-summary(lm_fit)
	  lm_p<-(unlist(lapply(lm_fit_summary,function(x){x$coefficients[,4][2]})))
		# lm_p<-for (item in lm_fit_summary){ print(item$coefficients[,4][2])}
		lm_resid<-residuals(lm_fit)
		lm_resid_t<-(t(lm_resid))
		lm_resid_t<-cbind("Id"=rownames(lm_resid_t),lm_resid_t)
		#write this to the new folder
                #write.table(lm_resid,file=paste0(resid_dir_path,"/residuals_sex.tsv"),sep="\t",quote=F)
		fileout=paste0(path,"/",tissue_name,".both.sex.peer.v8ciseQTLs.ztrans.txt")
		write.table(lm_resid_t,file=fileout,sep="\t",quote=F,row.names=F)
	#}


	return(used_inds)	
	
}

#checks to make sure parameters were properly passed into file
#check_params<-function(peer_dir,subset_do,residuals_do,metadata_file){
check_params<-function(peer_dir,metadata_file){
	if(length(peer_dir) == 0){stop("ERROR: peer file dir not provided")}
	if(substr(peer_dir,nchar(peer_dir),nchar(peer_dir))=="/"){stop("Remove trailing / from peer_dir for proper results")}
	#if(subset_do != "y" & subset_do != "yes" & subset_do != "n" & subset_do != "no"){
	#	stop("ERROR: subset_do must be the following: y, yes, n, no")
	#}
	#if(residuals_do != "y" & residuals_do != "yes" & residuals_do != "n" & residuals_do != "no"){
	#	stop("ERROR: residuals_do must be the following: y, yes, n, no")
	#}
	if(length(metadata_file) == 0){stop("ERROR: metadata file not provided")}
}

####MAIN

opt_parser<-get_opt_parser()
args<-parse_args(opt_parser)
peer_dir<-as.character(args$dir_peer)
#subset_do<-tolower(as.character(args$subset_do))
#residuals_do<-tolower(as.character(args$residuals_do))
metadata_file<-as.character(args$metadata_file)
file_to_correct<-as.character(args$file_to_correct)


#metadata_file<-"/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
#peer_dir<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/peer_v8"
#file_to_correct<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/Thyroid.both.sex.peer.v8ciseQTLs.ztrans.txt"

print(peer_dir)
print(metadata_file)

#check for parameters being passed properly, if not quit
#check_params(peer_dir,subset_do,residuals_do,metadata_file)
check_params(peer_dir,metadata_file)

#make a dictionary where keys are individuals, vales are sex
metadata<-fread(metadata_file,sep="\t",header=T)
ind_dict<-metadata$SEX
names(ind_dict)<-metadata$SUBJID

#peer_dirs_bothonly<-list.files("/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8",pattern="both.peer")


used_inds<-list() #individuals who have been included in other tissues
#for(file_path in peer_dirs_bothonly){
for(file_path in file_to_correct){
	print(file_path)
	if(file_path == peer_dir){next} #for some reason list.dirs prints current dir, so skip that
	#used_inds are individuals that have been subsetted, this loops so at the end it is a master list
	#used_inds is called by used_inds$tissue_name$female where female is male/female/both
	print(file_path)
	#used_inds<-get_individuals(tissue_dir,subset_do,residuals_do,ind_dict,used_inds)
	#used_inds<-get_individuals(tissue_dir,ind_dict,used_inds)
	#used_inds<-get_individuals(paste0(peer_dir,"/",file_path),ind_dict,used_inds)
	used_inds<-get_individuals(file_path,ind_dict,used_inds)
	#break
}

print("EOF")
