library(stringr) #str_replace
library(ggplot2) #plotting
library(reshape2) #melt

#args <- commandArgs(trailingOnly = TRUE)
# covariate_file<-args[1]
# peer_file<-args[2]
# output_file<-args[3]
setwd("/oak/stanford/groups/smontgom/raungar/Sex")
covariate_file<-"/oak/stanford/groups/smontgom/raungar/Sex/Files/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
pca_file<-"Files/covariates_pcs.txt"
peer_file_blood<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/whole_bload_factors_60.tsv"
peer_file_liver<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/liver_factors_30.tsv"
peer_file_lung<-"/oak/stanford/groups/smontgom/raungar/Sex/Output/lung_factors_60.tsv"
# x_ref_file<-"Files/gencode.v26.GRCh38.genes.Xonly.bed"
# output_file<-"Plots/peer_covar_correlations_xonly_r2_20.png"


eigen_cor_pcs<-function(done_pca,md,metavars,pcs=1:10,
                        titleX = '',cexTitleX = 1.0,rotTitleX = 0, colTitleX = 'black', fontTitleX = 2,
                        titleY = '', cexTitleY = 1.0,rotTitleY = 0,  colTitleY = 'black', fontTitleY = 2,
                        cexLabX = 1.0, rotLabX = 0, colLabX = 'black', fontLabX = 2, cexLabY = 1.0,
                        rotLabY = 0, colLabY = 'black', fontLabY = 2, posLab = 'bottomleft',
                        col = c('blue4', 'blue3', 'blue2', 'blue1', 'white',
                                'red1', 'red2', 'red3', 'red4'),
                        posColKey = 'right', cexLabColKey = 1.0, cexCorval = 1.0, colCorval = 'black',
                        fontCorval = 1, scale = TRUE, main = '',  cexMain = 2, rotMain = 0, colMain = 'black',
                        fontMain = 2, corFUN = 'pearson', corUSE = 'pairwise.complete.obs', corMultipleTestCorrection = 'none',
                        signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
                        colFrame = 'white', plotRsquared = FALSE, returnPlot = TRUE)
{
  data <- done_pca
  metadata <- md
  components<-paste("PC",pcs,sep="")
  
  # issue warning if any columns to use are not numeric --- This is kind of annoying
  # --- and there should be an option to force numeric conversion as anyway I'm doing it in advance
  for (i in seq_len(length(components))) {
    if(!is.numeric(data[,components[i]])) {
      warning(components[i],
              ' is not numeric - please check the source data',
              ' as everything will be converted to a matrix')
    }
  }
  for (i in seq_len(length(metavars))) {
    if(!is.numeric(metadata[,metavars[i]])) {
      warning(metavars[i],
              ' is not numeric - please check the source data',
              ' as non-numeric variables will be coerced to numeric')
    }
  }
  xvals <- data.matrix(data[,which(colnames(data) %in% components)])
  yvals <- metadata[,which(colnames(metadata) %in% metavars)]
  chararcter_columns = unlist(lapply(yvals, is.numeric))  
  # negate it - basically if it
  chararcter_columns = !chararcter_columns
  # select only the names that are true 
  chararcter_columns = names(which(chararcter_columns))
  for (c in chararcter_columns) {
    print(c)
    yvals[, eval(quote(c))] = as.numeric(as.factor(yvals[, eval(quote(c))]))
  }
  yvals<-data.matrix(yvals)
  
  corvals <- cor(xvals, yvals, use = corUSE, method = corFUN)
  # create a new df with same dimensions as corvals and fill with P values
  # total number of tests we perform
  N <- ncol(xvals) * ncol(yvals)
  pvals <- data.frame(pval = numeric(N),
                      i = numeric(N),
                      j = numeric(N))
  k <- 0
  for (i in seq_len(ncol(xvals))) {
    for (j in seq_len(ncol(yvals))) { 
      k <- k + 1
      pvals[k,'pval'] <- cor.test(xvals[,i],
                                  yvals[,j],
                                  use = corUSE,
                                  method = corFUN)$p.value
      pvals[k,"i"] <- colnames(xvals)[i]
      pvals[k,"j"] <- colnames(yvals)[j]
      
    }
  }
  print(3)
  ### code courtesy of aleighbrown
  # -----if you want to adjust the p-values for multiple testing
  if(corMultipleTestCorrection != "none"){
    pvals$pval <- p.adjust(pvals$pval, method = corMultipleTestCorrection)
  }
  
  pvals <- reshape2::dcast(pvals, i ~ j, value.var = "pval")
  # -----make sure the pvals matchs the order of corrvals table
  rownames(pvals) <- pvals$i
  pvals$i <- NULL
  pvals <- pvals[match(rownames(corvals), rownames(pvals)), ]
  # ---make sure the columns are in the correct order
  pvals <- pvals[colnames(corvals)]
  # ------
  ### END
  
  # are we plotting R^2 values?
  if (plotRsquared==TRUE) {
    corvals <- corvals ^ 2
  }
  
  # determine max and min correlation values in order to define the range
  if (scale == FALSE && plotRsquared == TRUE) {
    iUpperRange <- 1
    iLowerRange <- 0
  } else if (scale == FALSE && plotRsquared == FALSE) {
    iUpperRange <- 1
    iLowerRange <- -1
  } else if (scale == TRUE) {
    max <- max(corvals)
    min <- min(corvals)
    if(abs(max) > abs(min)) {
      iUpperRange <- max + 0.01
      iLowerRange <- (max * (-1)) - 0.01
    } else {
      iUpperRange <- abs(min) + 0.01
      iLowerRange <- min - 0.01
    }
    if (plotRsquared==TRUE) {
      iUpperRange <- max + 0.1
      iLowerRange <- 0
    }
  }
  
  # define the colour scheme/palette
  cols <- colorRampPalette(col)
  
  # create a new df with same dimensions as corvals
  # fill with significances encoded with asterisks
  signif <- corvals
  for (i in seq_len(ncol(pvals))) {
    signif[,i] <- c(symnum(pvals[,i],
                           corr = FALSE,
                           na = FALSE,
                           cutpoints = signifCutpoints,
                           symbols = signifSymbols))
  }
  
  # create a new df with same dimensions as corvals
  # fill with r values merged with the encoded significances
  plotLabels <- corvals
  for (i in seq_len(nrow(corvals))) {
    for(j in seq_len(ncol(corvals))) {
      plotLabels[i,j] <- paste(round(corvals[i,j], 2),
                               signif[i,j],
                               sep='')
      colnames(plotLabels)[j] <- colnames(corvals)[j]
    }
    
    rownames(plotLabels)[i] <- rownames(corvals)[i]
  }
  
  # position of axis ticks
  if (posLab == 'bottomleft') {
    posLab = 1
    axisTicks = c(1,0)
  } else if (posLab == 'topright') {
    posLab = 2
    axisTicks = c(0,1)
  } else if (posLab == 'all') {
    posLab = 3
    axisTicks = c(1,1)
  } else if (posLab == 'none') {
    posLab = 0
    axisTicks = c(0,0)
  }
  
  # define a panel function for adding labels
  # labels are passed with z as a third dimension
  labels <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    ltext(x, y,
          labels = plotLabels,
          cex = cexCorval,
          col = colCorval,
          font = fontCorval)
  }
  
  # produce the levelplot
  l <- levelplot(
    data.matrix(corvals),
    xlab = list(label = titleX,
                cex = cexTitleX,
                rot = rotTitleX,
                col = colTitleX,
                font = fontTitleX),
    ylab = list(label = titleY,
                cex = cexTitleY,
                rot = rotTitleY,
                col = colTitleY,
                font = fontTitleY),
    panel = labels,
    pretty = TRUE,
    par.settings = list(panel.background = list(col = colFrame)),
    scales = list(
      x = list(cex = cexLabX,
               rot = rotLabX,
               col = colLabX,
               font = fontLabX),
      y = list(cex = cexLabY,
               rot = rotLabY,
               col = colLabY,
               font = fontLabY),
      tck = axisTicks,
      alternating = posLab),
    aspect = 'fill',
    col.regions = cols,
    cuts = 100,
    at = seq(iLowerRange, iUpperRange, 0.01),
    main = list(label = main,
                cex = cexMain,
                rot = rotMain,
                col = colMain,
                font = fontMain),
    colorkey = list(space = posColKey,
                    labels = list(cex = cexLabColKey)))
  
  # return plot?
  if (returnPlot == TRUE) {
    return(l)
  } else if (returnPlot == FALSE) {
    l
  }
}


#READ IN FILES
covariates<-read.csv(covariate_file,header=T,sep="\t")
rownames(covariates)<-covariates[,1]
peer_res_blood<-read.csv(peer_file_blood,header=T,sep="\t")
rownames(peer_res_blood)<-peer_res_blood[,1]
colnames(peer_res_blood)<-str_replace(colnames(peer_res_blood),"\\.","-") 
peer_res_liver<-read.csv(peer_file_liver,header=T,sep="\t")
rownames(peer_res_liver)<-peer_res_liver[,1]
colnames(peer_res_liver)<-str_replace(colnames(peer_res_liver),"\\.","-") 
peer_res_lung<-read.csv(peer_file_lung,header=T,sep="\t")
rownames(peer_res_lung)<-peer_res_lung[,1]
colnames(peer_res_lung)<-str_replace(colnames(peer_res_lung),"\\.","-") 
#change col names to match covariates
#x_ref<-read.csv(x_ref_file,header = F,sep = "\t")
#rownames(x_ref)<-x_ref[,4]
pca_res<-read.csv(pca_file,header=T,sep="\t")
rownames(pca_res)<-pca_res[,1]
pca_res<-pca_res[,-c(1,22)]
#fix covariate file for correlation analysis
reduced_covars<-function(covariates){
  
  #convert hrs/mins string to mins numeric
  covariates$TRISCH<-as.numeric(str_split_fixed(covariates$TRISCH," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRISCH," ",n=4)[,3])
  covariates$TRCHSTIN<-as.numeric(str_split_fixed(covariates$TRCHSTIN," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRCHSTIN," ",n=4)[,3])
  covariates$TRCCLMP<- as.numeric(str_split_fixed(covariates$TRCCLMP," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$TRCCLMP," ",n=4)[,3])
  covariates$DTHPRNINT<- as.numeric(str_split_fixed(covariates$DTHPRNINT," ",n=4)[,1])*60+as.numeric(str_split_fixed(covariates$DTHPRNINT," ",n=4)[,3])

  #if didn't smoke, smoking years to zero not NA like why yo
  covariates[is.na(covariates$MHSMKYRS) & covariates$MHSMKSTS=="No","MHSMKYRS"]<-0
  covariates$MHSMKYRS<-as.numeric(covariates$MHSMKYRS)
  covariates[(covariates$MHDRNKNMB==99) & covariates$MHDRNKSTS=="Yes","MHDRNKYRS"]<-"NA"
  covariates[(covariates$MHDRNKSTS=="No"),"MHDRNKYRS"]<-0
  covariates$MHDRNKYRS<-as.numeric(covariates$MHDRNKYRS)
  covariates[(covariates$MHCOPD==99) & !is.na(covariates$MHCOPD),"MHCOPD"]<-"NA"
  covariates$MHCOPD<-as.numeric(covariates$MHCOPD)
  covariates[(covariates$MHBCTINF==99) & !is.na(covariates$MHBCTINF),"MHBCTINF"]<-"NA"
  covariates$MHBCTINF<-as.numeric(covariates$MHBCTINF)
  covariates[(covariates$MHNPHYS4W==99) & !is.na(covariates$MHNPHYS4W),"MHNPHYS4W"]<-"NA"
  covariates$MHNPHYS4W<-as.numeric(covariates$MHNPHYS4W)
  
    
  covariates$SEX<-as.numeric(covariates$SEX-1)

  
  # covariates$INCEXC<-as.factor(covariates$INCEXC)
  # covariates$DTHDTRMN<-as.factor(covariates$DTHDTRMN)
  # covariates$COHORT<-as.factor(covariates$COHORT)
  # covariates$TRCRTMPU<-as.factor(covariates$TRCRTMPU)
  # covariates$TRTPTREF<-as.factor(covariates$TRTPTREF)
  # covariates$TRVNTSR<-as.factor(covariates$TRVNTSR)
  # covariates$DTHTPTREF<-as.factor(covariates$DTHTPTREF)
  # covariates$DTHMNNR<-as.factor(covariates$DTHMNNR)
  # covariates$DTHRFGD<-as.factor(covariates$DTHRFGD)
  # covariates$DTHPLCE<-as.factor(covariates$DTHPLCE)
  # covariates$MHSRC<-as.factor(covariates$MHSRC)
  # covariates$DTHSEASON<-as.factor(covariates$DTHSEASON)
  # 
  covariates_red<-covariates #[,c(2:7,9,11:18,21,22,24:32,34,37,39,41,43,45:47,49:77,79:97,99:162,164:174)]
  
  
  dth_time<-as.numeric(str_replace_all(covariates$DTHTIME,":","\\."))
  dthszn<-model.matrix(~0+covariates$DTHSEASON)
  colnames(dthszn)<-paste0("DeathSeason_",c("NA","Fall","Spring","Summer","Winter"))
  covariates_red<-cbind(covariates_red,dthszn)
  dthplc<-model.matrix(~0+covariates$DTHPLCE)
  colnames(dthplc)<-paste0("DeathPlace_",c("NA",levels(covariates$DTHPLCE)[-1]))
  covariates_red<-cbind(covariates_red,dthplc)
  dthmnnr<-model.matrix(~0+covariates$DTHMNNR)
  colnames(dthmnnr)<-paste0("DeathManner_",c("NA",levels(covariates$DTHMNNR)[-1]))
  covariates_red<-cbind(covariates_red,dthmnnr)
  race<-model.matrix(~0+(as.factor(covariates$RACE)))
  colnames(race)<-paste0("Race_",c(levels(as.factor(covariates$RACE))))
  covariates_red<-cbind(covariates_red,race)
  cohort<-model.matrix(~0+(as.factor(covariates$COHORT)))
  colnames(cohort)<-paste0("Cohort_",c("NA",levels(covariates$COHORT)[-1]))
  covariates_red<-cbind(covariates_red,cohort)
  hardyDeath<-model.matrix(~0+factor(as.character(covariates$DTHHRDY),exclude=NULL))
  colnames(hardyDeath)<-paste0("hardyDeath_",c(levels(as.factor((covariates$DTHHRDY))),"NA"))
  covariates_red<-cbind(covariates_red,hardyDeath)
  
  
  covariates_red$dth_time_10_14<-as.numeric(dth_time >= 10 & dth_time <14)
  covariates_red$dth_time_14_18<-as.numeric(dth_time >= 14 & dth_time <18)
  covariates_red$dth_time_18_22<-as.numeric(dth_time >= 18 & dth_time < 22)
  covariates_red$dth_time_22_2<-as.numeric(dth_time >= 22 | dth_time <2)
  covariates_red$dth_time_2_6<-as.numeric(dth_time >= 2 & dth_time <6)
  covariates_red$dth_time_6_10<-as.numeric(dth_time >= 6 & dth_time <10)
  
  
  
  return(covariates_red)
}
covariates_red<-reduced_covars(covariates)
covariates_red_pcrows<-covariates_red[rownames(pca_res),]

metavars_to_use<-c("AGE","SEX","TRISCHD") #,"DTHRFGD","DTHVNTD", "HGHT","WGHT","BMI",
                  # "MHCOPD","MHBCTINF","MHNPHYS4W","MHSMKYRS","MHDRNKYRS" ) #,
                 #  "DTHPLCE","DTHMNNR","DTHSEASON","RACE","COHORT",
                 #  "dth_time_10_14","dth_time_14_18","dth_time_18_22","dth_time_22_2","dth_time_2_6","dth_time_6_10")
metavars_added<-c(227:232, which(colnames(covariates_red) %in% metavars_to_use))
covariates_red_pcrows$MHDRNKYRS<-as.numeric(covariates_red_pcrows$MHDRNKYRS)
eigen_cor_pcs(pca_res,covariates_red_pcrows,pcs = 1:20,
              metavars=colnames(covariates_red_pcrows)[metavars_added],
              scale=F,corMultipleTestCorrection = "fdr")
#covariates w/o NAs
covariates_red_rmna<-covariates_red[ , colSums(is.na(covariates_red)) == 0]
#only get peer with x chr genes
peer_xonly<-peer_res[rownames(peer_res) %in% rownames(x_ref),-1]


#correlate to itself!! :)
cor_metavar_self<-corAndPvalue(covariates_red[-nrow(covariates_red),metavars_added],
                               method = "spearman",use="everything")
library(ggcorrplot)
library(WGCNA)
ggcorrplot(cor_metavar_self$cor,tl.cex = 6.5,lab = T)


###
eigencorplot(pca_res,
             metavars=covariates_red[1:10],
             titleX = "PCs RNA-seq", titleY = "PCs WGS",rotTitleY = 90,scale = F)

###just to check
glm_hardyDeath0<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_0,family = "binomial")
glm_hardyDeath1<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_1,family = "binomial")
glm_hardyDeath2<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_2,family = "binomial")
glm_hardyDeath3<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_3,family = "binomial")
glm_hardyDeath4<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_4,family = "binomial")
glm_hardyDeathNA<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$hardyDeath_NA,family = "binomial")
glm_DTHHRDY<-glm(covariates_red_pcrows$SEX~as.factor(covariates_red_pcrows$DTHHRDY),family = "binomial")
glm_trisch<-glm(covariates_red_pcrows$SEX~covariates_red_pcrows$TRISCH,family = "binomial")

###only for right now
covariates_red_blood<-covariates_red[rownames(covariates_red) %in% colnames(peer_res_blood),]
covariates_red_blood<-covariates_red_blood[order(match(rownames(covariates_red_blood), 
                                                       colnames(peer_res_blood))), ]
covariates_red_liver<-covariates_red[rownames(covariates_red) %in% colnames(peer_res_liver),]
covariates_red_liver<-covariates_red_liver[order(match(rownames(covariates_red_liver), 
                                                       colnames(peer_res_liver))), ]
covariates_red_lung<-covariates_red[rownames(covariates_red) %in% colnames(peer_res_lung),]
covariates_red_lung<-covariates_red_lung[order(match(rownames(covariates_red_lung), 
                                                         colnames(peer_res_lung))), ]
peer_cov_corr_blood<-cor(covariates_red_blood[,metavars_added],t(peer_res_blood[,-1]),method = "spearman",use = "everything")
peer_cov_corr_blood_pairwise<-cor(covariates_red_blood[,metavars_added],t(peer_res_blood[,-1]),method = "spearman",use = "pairwise.complete.obs")
peer_cov_corr_liver<-cor(covariates_red_liver[,metavars_added],t(peer_res_liver[,-1]),method = "spearman",use = "everything")
peer_cov_corr_liver_pairwise<-cor(covariates_red_liver[,metavars_added],t(peer_res_liver[,-1]),method = "spearman",use = "pairwise.complete.obs")
peer_cov_corr_lung<-cor(covariates_red_lung[,metavars_added],t(peer_res_lung[,-1]),method = "spearman",use = "everything")
peer_cov_corr_lung_pairwise<-cor(covariates_red_lung[,metavars_added],t(peer_res_lung[,-1]),method = "spearman",use = "pairwise.complete.obs")

corr_df<-melt(peer_cov_corr_lung)
# corr_df$r<-as.numeric(as.character(corr_df$r))

#R^2 greater than five to reduce variables
corr_df_min_10<-corr_df[abs(corr_df$value)>0.1,]
corr_df_min_20<-corr_df[abs(corr_df$value)>0.2,]

my_hmp<-ggplot(corr_df, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 2)),size=1.5) +
  scale_fill_gradient2(low="blue", high="red",midpoint=0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=6),
        axis.text.y = element_text( hjust = 1,size=6),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  labs(fill="R^2") +
  ggtitle("Lung Peer vs. Covariates")

ggsave("Plots/peer_lung_covars_cor.png",my_hmp,dpi=300)  

png(output_file)
plot(my_hmp)
dev.off()
