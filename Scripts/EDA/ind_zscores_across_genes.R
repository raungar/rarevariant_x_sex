library(data.table)



get_med_mean<-function(this_dir,this_pattern) {
  #get gtex ids
  inds<-c()
  for(file in list.files(path=this_dir,pattern=this_pattern)){
    #quick way to just get headers
    these_inds<-scan(paste0(dir,file),nlines=1,what=character()) 
    inds<-union(inds,these_inds)
  }
  
  #initialize variables
  inds_dic<-rep(NA,length(inds[-1]))
  names(inds_dic)<-inds[-1]
  exp_summary_mean<-data.frame(matrix(ncol=0,nrow=length(inds_dic)))
  rownames(exp_summary_mean)<-names(inds_dic)
  exp_summary_med<-data.frame(matrix(ncol=0,nrow=length(inds_dic)))
  rownames(exp_summary_med)<-names(inds_dic)
  
  i=1
  for(file in list.files(path=dir,pattern="*.both_half_regress.peer.ztrans.txt")){
    tissue<-sapply(strsplit(file,"\\."),"[[",1)
    #reset
    inds_dic_med<-inds_dic
    inds_dic_mean<-inds_dic
    
    #read and calculate
    exp<-fread(paste0(dir,file))
    if( (ncol(exp)-1)<50){next}
    
    exp_no_genes<-exp[,-1]
    exp_mean<-colMeans(exp_no_genes)
    exp_median<-colMedians(as.matrix(exp_no_genes))
    
    #paste into dataframe
    inds_dic_mean[names(exp_mean)]<-exp_mean
    inds_dic_med[names(exp_mean)]<-exp_median
    
    # exp_summary_mean<-rbind(exp_summary_mean,t(data.frame(inds_dic_mean)))
    # rownames(exp_summary_mean)[nrow(exp_summary_mean)]<-tissue
    exp_summary_mean<-cbind(exp_summary_mean,data.frame(inds_dic_mean))
    colnames(exp_summary_mean)[ncol(exp_summary_mean)]<-tissue
    exp_summary_med<-cbind(exp_summary_med,data.frame(inds_dic_med))
    colnames(exp_summary_med)[ncol(exp_summary_med)]<-tissue
    # if(i==5){break}; i=i+1
  }
  return(list(exp_summary_med,exp_summary_mean))
}

dir="/oak/stanford/groups/smontgom/raungar/Sex/Output/preprocessing_v8/PEER_v8/"
res_both<-get_med_mean(dir,"*.both_half_regress.peer.ztrans.txt")
exp_summary_med_both<-res_both[[1]]
exp_summary_mean_both<-res_both[[2]]


get_med_mean_summ<-function(exp_summ){
  t_exp_summary_mean<-t(exp_summ)
  t_exp_summary_meanofmean<-colMeans(t_exp_summary_mean,na.rm = T)
  t_exp_summary_medofmean<-colMedians(t_exp_summary_mean,na.rm = T)
  names(t_exp_summary_meanofmean)<-colnames(t_exp_summary_mean)
  names(t_exp_summary_medofmean)<-colnames(t_exp_summary_mean)
  exp_summary_meanofmean<-na.omit(t_exp_summary_meanofmean)
  exp_summary_medofmean<-na.omit(t_exp_summary_medofmean)
  return(list(exp_summary_medofmean,exp_summary_meanofmean))
}
exp_summary_res_both_mean<-get_med_mean_summ(exp_summary_mean_both)
both_medianofmean<-exp_summary_res_both_mean[[1]]
both_meanofmean<-exp_summary_res_both_mean[[2]]
exp_summary_res_both_med<-get_med_mean_summ(exp_summary_med_both)
both_medianofmedian<-exp_summary_res_both_med[[1]]
both_meanofmedian<-exp_summary_res_both_med[[2]]

t_exp_summary_med<-t(exp_summary_med_both)
t_exp_summary_meanofmed<-colMeans(t_exp_summary_med,na.rm = T)
t_exp_summary_medofmed<-colMedians((t_exp_summary_med),na.rm = T)
names(t_exp_summary_meanofmed)<-colnames(t_exp_summary_med)
names(t_exp_summary_medofmed)<-colnames(t_exp_summary_med)
exp_summary_meanofmed<-na.omit(t_exp_summary_meanofmed)
z_meanofmed<-scale(exp_summary_meanofmed)
exp_summary_medofmed<-na.omit(t_exp_summary_medofmed)
z_medofmed<-scale(exp_summary_medofmed)
