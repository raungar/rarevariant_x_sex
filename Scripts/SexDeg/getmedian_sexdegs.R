library(dplyr)
library(data.table)

###get summary of the betas

all_betas<-c()
for (file in list.files("Files/SexDEGs/")){

	this_f<-fread(paste0("Files/SexDEGs/",file))
	all_betas<-c(all_betas,this_f$`MASH beta`)

}

print(summary(abs(all_betas)))
num=0.04954
all_betas_gt0.04954<-all_betas[abs(all_betas)>num]
print(summary(abs(all_betas_gt0.04954)))
