library(data.table)

args = commandArgs(trailingOnly = T)
outfile=args[1] #first file MUST BE OUTFILE
args=args[2:length(args)] #trailing list of files from expand
all_df<-data.frame()
#loop through all tissue files and merge
for(arg in args){
    tissue=sapply(strsplit(arg,"/"),tail,1)
    tissue=sapply(strsplit(tissue,"\\."),"[[",1)
    print(tissue)
    this_df=fread(arg)
    this_df$tissue=tissue
    if(nrow(all_df)==0){
        all_df=this_df
    }else{
        merge_by=colnames(all_df)[colnames(all_df)%in%colnames(this_df)]
        #print(head(merge_by))
        all_df<-merge(all_df,this_df,by=merge_by,all=T)
    }
    print(dim(all_df))
    #print(colnames(all_df))
}
all_df[all_df==""]<-NA
fwrite(all_df,file=outfile,sep="\t")