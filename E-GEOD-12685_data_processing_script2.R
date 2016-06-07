
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                   AD - E-GEOD-12685  - PIEPLINE V02  **ONLY MCI SAMPLES***                                             #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Afymetrix
# EXPRESSION CHIP - HG-U133A
# NUMBER OF SAMPLES - 8 control - 6 Incipient AD (MCI)
# BRAIN REGIONS - frontal cortex
# 

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

raw_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/E-GEOD-12685/Raw_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-12685/Pre-processing"

setwd(work_dir)

dir.create(paste(work_dir,"pca_plots", sep="/"))

pca_dir=paste(work_dir,"pca_plots", sep="/")

dir.create(paste(work_dir,"clean_data", sep="/"))

clean_data_dir=paste(work_dir,"clean_data", sep="/")

dir.create(paste(work_dir,"sample_network_plots", sep="/"))

sample_network_dir=paste(work_dir,"sample_network_plots", sep="/")

##### LOAD LIBRARIES ####

library(ArrayExpress)
library(affy)
library(lumi)
library(biomaRt)
library(WGCNA)
library(pamr)
library(limma)
library(sva)
library(hgu133a.db)
library(ggplot2)
library(reshape)
library(massiR)
library(stringi)

##### DOWNLOAD RAW DATA #####

setwd(raw_dir)

brain_data_raw=getAE("E-GEOD-12685", type = "full")

##### CREATE R EXPRESSION OBJECT FROM RAW DATA #####

#convert MAGE-TAB files into expresssion set

brain_data = ae2bioc(mageFiles = brain_data_raw)

brain_data

##### PLOTS OF RAW DATA ####

setwd(work_dir)

boxplot(exprs(brain_data))
pdf(file="raw_data_boxplot.pdf")
boxplot(exprs(brain_data))
dev.off()

plotDensity(exprs(brain_data), logMode=F, addLegend=F)
pdf(file="raw_data_density_plot.pdf")
plotDensity(exprs(brain_data), logMode=F, addLegend=F)
dev.off()

##### PRE-PROCESS #####

#background correct

brain_data_background_corrected<-mas5(brain_data)

brain_data_normalised_as_data_frame<-as.data.frame(brain_data_normalised)

#normalise

brain_data_normalised<-rsn(log2(exprs(brain_data_background_corrected)))

##### PLOTS OF PRE_PROCESSED DATA #####

setwd(work_dir)

boxplot(brain_data_normalised)
pdf(file="pre-processed_data_boxplot.pdf")
boxplot(brain_data_normalised)
dev.off()

plotDensity(brain_data_normalised, logMode=F, addLegend=F)
pdf(file="pre-processed_data_density_plot.pdf")
plotDensity(brain_data_normalised, logMode=F, addLegend=F)
dev.off()

##### EXTRACT PHENOTYPE INFO ####

phenotype_data<-pData(brain_data)

head(phenotype_data)

names(phenotype_data)

table(phenotype_data$Hybridization.Name)

table(phenotype_data$Comment..Sample_source_name.)

colnames(phenotype_data)

Diagnosis<-phenotype_data[17]

colnames(Diagnosis)<-"Diagnosis"

Diagnosis

Diagnosis_case<-subset(Diagnosis, grepl("AD", Diagnosis))

Diagnosis_control<-subset(Diagnosis, grepl("Control", Diagnosis))

Diagnosis_case

Diagnosis_control

Diagnosis_case[,1]<-"Incipient_AD"

Diagnosis_control[,1]<-"Control"

Diagnosis_case

Diagnosis_control

##### SAVE PROBE IDs #####

probe_ids<-rownames(brain_data_normalised)

head(probe_ids)
length(probe_ids)

setwd(work_dir)

write(file = "E-GEOD-12685_probe_IDs.txt", probe_ids, sep="\t")

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(phenotype_data))

##### GENDER CHECK ##### 
phenotype_data

# gender_information - unknown
gender_info<-phenotype_data[1]
colnames(gender_info)<-"Gender"
gender_info$Gender<-"Unknown"
table(gender_info$Gender)

# get Y choromosome genes
data(y.probes)
names(y.probes)

y_chromo_probes <- data.frame(y.probes["affy_hg_u133_plus_2"])

# extract Y chromosome genes from dataset
eset.select.out <- massi_select(brain_data_normalised_as_data_frame, y_chromo_probes, threshold=4)

massi_y_plot(eset.select.out)
massi_cluster_plot(eset.select.out)

# run gender predict
eset.results <- massi_cluster(eset.select.out)

#extract gender prediction
predicted_gender<-(eset.results$massi.results)[c(1,5)]
rownames(predicted_gender)<-predicted_gender$ID
predicted_gender$ID<-NULL
colnames(predicted_gender)<-"Predicted_Gender"

#compare to clinical Gender
#standardise

# gender_info$Gender<-as.character(gender_info$Gender)
# gender_info[gender_info$Gender=="M",]<-"male"
# gender_info[gender_info$Gender=="F",]<-"female"

#merge
gender_comparison<-merge(gender_info, predicted_gender, by="row.names")
rownames(gender_comparison)<-gender_comparison$Row.names
gender_comparison$Row.names<-NULL
colnames(gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
head(gender_comparison)

#separae male/female IDs
female_samples<-subset(gender_comparison, Predicted_Gender=="female")
male_samples<-subset(gender_comparison, Predicted_Gender=="male")

head(female_samples)
head(male_samples)

##### PROBE ID DETECTION #####

# separate case control - Factor.Value..disease. column 

Diagnosis_case

Diagnosis_control

case_ID<-rownames(Diagnosis_case)

control_ID<-rownames(Diagnosis_control)

case_ID

control_ID

case_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%case_ID]
control_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%control_ID]

head(case_exprs)
dim(case_exprs)
head(control_exprs)
dim(control_exprs)

# separate by gender

case_exprs_F<-case_exprs[colnames(case_exprs)%in%rownames(female_samples)]
case_exprs_M<-case_exprs[colnames(case_exprs)%in%rownames(male_samples)]

control_exprs_F<-control_exprs[colnames(control_exprs)%in%rownames(female_samples)]
control_exprs_M<-control_exprs[colnames(control_exprs)%in%rownames(male_samples)]


# calculate 90th percentile for each sample in each group

extract_good_probe_list<-function(dataset, probe_percentile_threshold, sample_threshold) {
  # dataset - expression dataset as dataframe
  # probe_percentile_threshold - percentile at which to use as cut-off for detected probes = 90th (compared to negative bead controls) - 0.9 = 90th percentile
  # number of samples in which probe must be expressed in - i.e 0.8 = 80% of samples
  # calculate quantile threshold for each sample
  sample_quantiles<-apply(dataset, 2, quantile, probs=c(probe_percentile_threshold))
  # count length of quantile - will be number of samples
  number_of_samples<-length(sample_quantiles)
  # convert probes values to NA in for each sample if probe expression value below sample percentile cut-off
  for (x in 1:number_of_samples) {
    is.na(dataset[x]) <- dataset[x] >= sample_quantiles[x]
  }
  # convert to dataframe
  dataset_dataframe<-as.data.frame(dataset)
  # count number of NA
  dataset_count<-as.data.frame(rowSums(is.na(dataset_dataframe)))
  colnames(dataset_count)<-"count"
  # subset good probes
  good_probes<-rownames(subset(dataset_count, dataset_count$count >= (number_of_samples*sample_threshold)))
  #print threshold used
  print(as.data.frame(sample_quantiles))
  boxplot(as.data.frame(sample_quantiles))
  # return good probes
  return(good_probes)
}

# apply function to case samples -**** percentile cut-off dropped to 0.7 based on number of probes selected + XIST gene expression in females is low

case_exprs_F_expressed_probes_list<-extract_good_probe_list(case_exprs_F, 0.7, 0.8)
length(case_exprs_F_expressed_probes_list)

case_exprs_M_expressed_probes_list<-extract_good_probe_list(case_exprs_M, 0.7, 0.8)
length(case_exprs_M_expressed_probes_list)

control_exprs_F_expressed_probes_list<-extract_good_probe_list(control_exprs_F, 0.7, 0.8)
length(control_exprs_F_expressed_probes_list)

control_exprs_M_expressed_probes_list<-extract_good_probe_list(control_exprs_M, 0.7, 0.8)
length(control_exprs_M_expressed_probes_list)

# merge list of good probes from both case + control, sort and keep unique values

good_probe_list<-unique(sort(c(case_exprs_F_expressed_probes_list,
                               case_exprs_M_expressed_probes_list,
                               control_exprs_F_expressed_probes_list,
                               control_exprs_M_expressed_probes_list)))

length(good_probe_list)

# extract good probes from dataset

brain_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]

brain_case_exprs_good_probes<-case_exprs[rownames(case_exprs)%in%good_probe_list,]

brain_control_exprs_good_probes<-control_exprs[rownames(control_exprs)%in%good_probe_list,]

head(brain_exprs_good_probes)[1:5]
dim(brain_data_normalised)
dim(brain_exprs_good_probes)
dim(brain_case_exprs_good_probes)
dim(brain_control_exprs_good_probes)

##### GENDER SPECIFIC PROBE PLOTS #####

# using dataframe before probe removal

# get gene symbol list for chip
Gene_symbols_probes <- mappedkeys(hgu133aSYMBOL)

# Convert to a list
Gene_symbols <- as.data.frame(hgu133aSYMBOL[Gene_symbols_probes])

head(Gene_symbols)
dim(Gene_symbols)

#xist gene - 
XIST_probe_ID<-subset(Gene_symbols, symbol=="XIST")
XIST_probe_ID

PRKY_probe_ID<-subset(Gene_symbols, symbol=="PRKY")
PRKY_probe_ID

NLGN4Y_probe_ID<-subset(Gene_symbols, symbol=="NLGN4Y")
NLGN4Y_probe_ID

TMSB4Y_probe_ID<-subset(Gene_symbols, symbol=="TMSB4Y")
TMSB4Y_probe_ID

USP9Y_probe_ID<-subset(Gene_symbols, symbol=="USP9Y")
USP9Y_probe_ID

UTY_probe_ID<-subset(Gene_symbols, symbol=="UTY")
UTY_probe_ID

# merge all genes 
gene_list<-rbind(XIST_probe_ID,
                 PRKY_probe_ID,
                 NLGN4Y_probe_ID,
                 TMSB4Y_probe_ID,
                 USP9Y_probe_ID,
                 UTY_probe_ID)

#create function to plot
plot_gender_specific_genes<-function(Expression_table, gender_info, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[rownames(Expression_table)%in% genes_to_extract$probe_id,]))
  #check all probes extracted
  print(c("all probes extracted:", dim(Expression_table_gene_check)[2]==dim(genes_to_extract)[1]))
  # change colnames TO GENE SYMBOL using genes to extract file
  for (x in 1:dim(Expression_table_gene_check)[2]){
    colnames(Expression_table_gene_check)[x]<-gene_list[genes_to_extract$probe_id==colnames(Expression_table_gene_check)[x],2]
  }
  # add in gender information
  Expression_table_gene_check_gender<-merge(gender_info, Expression_table_gene_check, by="row.names")
  rownames(Expression_table_gene_check_gender)<-Expression_table_gene_check_gender$Row.names
  Expression_table_gene_check_gender$Row.names<-NULL
  #melt dataframe for plot
  Expression_table_gene_check_gender_melt<-melt(Expression_table_gene_check_gender, by=Gender)
  # calculate user defined percentie threshold
  sample_quantiles<-apply(Expression_table, 2, quantile, probs=threshold)
  #print sample quantiles
  print(as.data.frame(sample_quantiles))
  # mean of used defined threshold across samples
  mean_threshold=mean(sample_quantiles)
  #plot
  qplot(variable, value, colour=get(colnames(gender_info)), data = Expression_table_gene_check_gender_melt, geom = c("boxplot", "jitter")) + 
    geom_hline(yintercept = mean_threshold) +
    ggtitle(boxplot_title) +
    labs(x="Gene",y="Expression", colour = colnames(gender_info)) 
}

# plot

setwd(work_dir)

dir.create(paste(work_dir,"Gender_specific_gene_plots", sep="/"))
Gender_plots_dir=paste(work_dir,"Gender_specific_gene_plots", sep="/")

setwd(Gender_plots_dir)

pdf("Predicted_Gender_specific_gene_plot_and_detectable_09_cut_off_threshold_used.pdf")
plot_gender_specific_genes(case_exprs, gender_comparison[2], gene_list, 0.9, "Prefrontal_cortex_case_0.9_cut-off")
plot_gender_specific_genes(control_exprs, gender_comparison[2], gene_list, 0.9, "Prefrontal_cortex_exprs_0.9_cut-off")
dev.off()

pdf("Predicted_Gender_specific_gene_plot_and_detectable_07_cut_off_threshold_used.pdf")
plot_gender_specific_genes(case_exprs, gender_comparison[2], gene_list, 0.7, "Prefrontal_cortex_case_0.7_cut-off")
plot_gender_specific_genes(control_exprs, gender_comparison[2], gene_list, 0.7, "Prefrontal_cortex_0.7_cut-off")
dev.off()

##### PCA ####

# calculate pca

pca<-prcomp(brain_exprs_good_probes)

# plot variance
plot(pca, type="l")

# sumary of pca
summary_pca<-summary(pca)

# pca importance
pca_importance_var_exp<-summary_pca$importance[2,]
pca_importance_var_exp_cum<-summary_pca$importance[3,]

par(mfrow=c(1,2))

#plot variance explained
plot(pca_importance_var_exp, ylab="PCA Proportion of Variance Explained", type="b", col="blue")

#plot variance explained cumalative
plot(pca_importance_var_exp_cum, ylab="PCA Cumulative Proportion of Variance Explained", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)

par(dev.off)
dev.off()

# pca matrix plot

#color
Diagnosis_lookup<-rbind(Diagnosis_case, Diagnosis_control)
Diagnosis_lookup


# order of samples in expression data
Diagnosis_temp<-colnames(brain_exprs_good_probes)
Diagnosis_temp

# match order
Diagnosis_pca<-Diagnosis_lookup[match(Diagnosis_temp, rownames(Diagnosis_lookup)),]
Diagnosis_pca

# assign color to group
Diagnosis_pca_color<-labels2colors(as.character(Diagnosis_pca))

# pca plot
plot(pca$rotation[,1:2], main=" PCA plot coloured by chip before_QC",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomright', unique(Diagnosis_pca), fill=unique(Diagnosis_pca_color))

setwd(work_dir)

##### SVA ####

dim(brain_case_exprs_good_probes)
dim(brain_control_exprs_good_probes)

brain_case_exprs_good_probes<-cbind(Diagnosis ="AD", brain_case_exprs_good_probes)
brain_control_exprs_good_probes<-cbind(Diagnosis = "CONTROL", brain_control_exprs_good_probes)

head(brain_case_exprs_good_probes)[1:5]
head(brain_control_exprs_good_probes)[1:5]

# check order of gene in each dataframe

any(colnames(brain_case_exprs_good_probes)==colnames(brain_control_exprs_good_probes))==F

# merge case control in to one dataset per brain region

Frontal_Cortex_good_probes<-rbind(brain_case_exprs_good_probes, brain_control_exprs_good_probes)
head(Frontal_Cortex_good_probes)[1:5]

# add Predicted gender in 
Frontal_Cortex_good_probes<-merge(gender_comparison[2], Frontal_Cortex_good_probes, by="row.names")
rownames(Frontal_Cortex_good_probes)<-Frontal_Cortex_good_probes$Row.names
Frontal_Cortex_good_probes$Row.names<-NULL
head(Frontal_Cortex_good_probes)[1:5]


# create sva function - useses predicted gender

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis+Predicted_Gender, data=dataset_pheno)
  # check number of SV in data
  num.sv(dataset_exprs, mod, method="leek")
}

# check sv

check_SV_in_data(Frontal_Cortex_good_probes)

# create function to adjust for SV - uses predicted gender

adjust_for_sva<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_sva_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_sva_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis+Predicted_Gender, data=dataset_sva_pheno)
  mod0 = model.matrix(~1, data=dataset_sva_pheno)
  # number of SV
  num.sv(dataset_sva_exprs, mod, method="leek")
  n.sv=num.sv(dataset_sva_exprs, mod, method="leek")
  # exit if n.sv=0
  if(n.sv==0){stop("No Significant Variable found, exiting....")}
  # apply sva - removed n.sv
  svobj = sva(dataset_sva_exprs, mod, mod0, n.sv=n.sv, method="two-step")
  # adjust for sva
  X = cbind(mod, svobj$sv)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(dataset_sva_exprs))
  P = ncol(mod)
  clean_data<-dataset_sva_exprs - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  # merge clean data with pheno
  clean_data_with_pheno<-merge(dataset_sva_pheno, as.data.frame(t(clean_data)), by="row.names")
  rownames(clean_data_with_pheno)<-clean_data_with_pheno$Row.names
  clean_data_with_pheno$Row.names<-NULL
  # check SVA on adjusted data
  cat("\n")
  cat("number of surrogate variables after adjustment:")
  cat("\n")
  print(num.sv(clean_data, mod, method="leek"))
  # return clean data with pheno
  return(clean_data_with_pheno)
}

# run SVA

Frontal_Cortex_good_probes_sva_adjusted<-adjust_for_sva(Frontal_Cortex_good_probes)

# separate case and control

Frontal_Cortex_good_probes_sva_adjusted_case<-Frontal_Cortex_good_probes_sva_adjusted[Frontal_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]
Frontal_Cortex_good_probes_sva_adjusted_control<-Frontal_Cortex_good_probes_sva_adjusted[Frontal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL",]

dim(Frontal_Cortex_good_probes_sva_adjusted_case)
dim(Frontal_Cortex_good_probes_sva_adjusted_control)

# remove gender
Frontal_Cortex_good_probes_sva_adjusted_case$Predicted_Gender<-NULL
Frontal_Cortex_good_probes_sva_adjusted_control$Predicted_Gender<-NULL

##### SAMPLE NETWORK PLOT #####

# sample plot function - taken from steve expression pipeline

sampleNetwork_plot <- function(dataset) {
  datExprs<-t(dataset[2:dim(dataset)[2]])
  diagnosis<-dataset[1]
  gp_col <- "group"
  cat(" setting up data for qc plots","\r","\n")
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- colnames(datExprs)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2 ## ADJACENCY MATRIX
  cat(" fundamentalNetworkConcepts","\r","\n")
  FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
  K2=FNC$ScaledConnectivity
  Z.K=(K2-mean(K2))/sd(K2)
  Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
  rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
  rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
  # set colours
  cat(" colorvec [",paste(gp_col),"]","\r","\n")
  if(gp_col=="chip") { colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode)) }
  if(gp_col=="group") { colorvec <- labels2colors(diagnosis[1]) }
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  ## samplenetwork
  local(
    {colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
      }
      n
    }
    i <- 0
    })
  cat(" begin SampleNetwork plots","\r","\n")
  group_colours<-unique(cbind(colorvec, diagnosis))
  ## Cluster for pics
  cluster1 <- hclust(as.dist(1-A.IAC),method="average")
  cluster1order <- cluster1$order
  cluster2 <- as.dendrogram(cluster1,hang=0.1)
  cluster3 <- dendrapply(cluster2,colLab,cluster1order)
  ## PLOTS
  ## cluster IAC
  par(mfrow=c(2,2))
  par(mar=c(5,6,4,2))
  plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
       main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
       xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
  mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
  ## Connectivity
  par(mar=c(5,5,4,2))
  plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
  text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
  abline(h=-2)
  abline(h=-3)
  par(mar=c(5,5,4,2))
  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec ,cex.main=1.8,cex.lab=1.4)
  abline(lm(Z.C~Z.K),col="black",lwd=2)
  mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
  abline(v=-2,lty=2,col="grey")
  abline(h=-2,lty=2,col="grey")
  ##blank plot for legend
  par(mar=c(5,5,4,2))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend(0.6, 1.4, unique(diagnosis[,1]), fill=unique(colorvec))
} #taken from steves expression pipeline

# create functio to ID outliers

names_of_outliers<-function(dataset, threshold){
  datExprs<-t(dataset[2:dim(dataset)[2]])
  IAC = cor(datExprs, method = "p", use = "p")
  diag(IAC) = 0
  A.IAC = ((1 + IAC)/2)^2  ## ADJACENCY MATRIX
  # fundamentalNetworkConcepts
  FNC = fundamentalNetworkConcepts(A.IAC)  ## WGCNA
  K2 = FNC$ScaledConnectivity
  Z.K = round((K2 - mean(K2))/sd(K2), 3)
  Z.K_outliers <- Z.K < threshold
  Z.K_outliers <- names(Z.K_outliers[Z.K_outliers == TRUE])
  n_outliers <- length(Z.K_outliers)
  return(Z.K_outliers)
}

# create function to run network analysis on each expression dataset, plot and remove bad samples 

run_sample_network_plot<-function(dataset, threshold){
  #sample network plot
  sampleNetwork_plot(dataset)
  #identify sample below Z.K -2  
  dataset_removal_1<-names_of_outliers(dataset, threshold)
  #remove samples with ZK below -2
  dataset_QC<-dataset[!(rownames(dataset)%in%dataset_removal_1),]
  #sample network plot
  sampleNetwork_plot(dataset_QC)
  #create empty count list to record samples removed
  count<-dataset_removal_1
  # reiterate above till no samples fall below threshold
  while (length(dataset_removal_1)>0) {
    # remove bad samples - 1st iteration removes none
    dataset_QC<-dataset_QC[!(rownames(dataset_QC)%in%dataset_removal_1),]
    #identify sample below Z.K -2  
    dataset_removal_1<-names_of_outliers(dataset_QC, threshold)
    #record samples removed
    count<-c(count, dataset_removal_1)
  }
  #final network plot
  sampleNetwork_plot(dataset_QC)
  # print to screen number of samples removed
  cat("\n")
  print(c("Total number of samples removed...", length(count)))
  # return clean expression set
  return(dataset_QC)
}

# run sample network on entorhinal Cortex

Frontal_Cortex_case_exprs_good_probes_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_case, -3)
Frontal_Cortex_control_exprs_good_probes_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("control_sample_network_analysis.pdf")
Frontal_Cortex_case_exprs_good_probes_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("case_sample_network_analysis.pdf")
Frontal_Cortex_control_exprs_good_probes_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

##### CREATE QC'd DATASET #####

brain_data_QCd<-Frontal_Cortex_good_probes_sva_adjusted
brain_data_QCd$Diagnosis<-NULL

dim(brain_exprs_good_probes)
dim(brain_data_QCd)

##### PCA ON CLEAN DATA ####

# calculate pca2

pca2<-prcomp(t(brain_data_QCd))

# plot variance
plot(pca2, type="l")

# sumary of pca2
summary_pca2<-summary(pca2)

# pca2 importance
pca2_importance_var_exp<-summary_pca2$importance[2,]
pca2_importance_var_exp_cum<-summary_pca2$importance[3,]

par(mfrow=c(1,2))

#plot variance explained
plot(pca2_importance_var_exp, ylab="pca Proportion of Variance Explained after QC", type="b", col="blue")

#plot variance explained cumalative
plot(pca2_importance_var_exp_cum, ylab="pca Cumulative Proportion of Variance Explained after QC", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)

par(dev.off)
dev.off()

# pca2 matrix plot

#color

# order of samples in expression data
Diagnosis_temp2<-colnames(t(brain_data_QCd))
Diagnosis_temp2

# match order
Diagnosis_pca2<-Diagnosis_lookup[match(Diagnosis_temp, rownames(Diagnosis_lookup)),]

# assign color to group
Diagnosis_pca2_color<-labels2colors(as.character(Diagnosis_pca2))

# pca2 plot
plot(pca2$rotation[,1:2], main=" pca plot coloured by chip after QC",col="black", pch=21,bg=Diagnosis_pca2_color)
legend('bottomleft', unique(Diagnosis_pca2), fill=unique(Diagnosis_pca2_color))

#plot to pdf

setwd(pca_dir)

pdf("pca_plot_before_and_after_QC.pdf")
plot(pca$rotation[,1:2], main=" PCA plot coloured by chip before_QC",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomleft', unique(Diagnosis_pca), fill=unique(Diagnosis_pca_color))
# pca2 plot
plot(pca2$rotation[,1:2], main=" pca plot coloured by Disease after QC",col="black", pch=21,bg=Diagnosis_pca2_color)
legend('bottomleft', unique(Diagnosis_pca2), fill=unique(Diagnosis_pca2_color))
dev.off()

setwd(work_dir)

##### CONVERT PROBE ID TO ENTREZ ID #####

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using hgu133a.db
mapped_probes <- mappedkeys(hgu133aENTREZID)

# Convert to a list
hgu133a.db_mapping <- as.data.frame(hgu133aENTREZID[mapped_probes])
# arrange order of column by entrezgene probe_id
hgu133a.db_mapping<-hgu133a.db_mapping[c(2,1)]
colnames(hgu133a.db_mapping)[1]<-"entrezgene"

head(hgu133a.db_mapping)
dim(hgu133a.db_mapping)

#check any duplicated probe IDs
anyDuplicated(hgu133a.db_mapping$probe_id)

#check any dupliacted entrezgene IDs
anyDuplicated(hgu133a.db_mapping$entrezgene)

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset # - removed this step
  expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}


# using hgu133a

brain_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(brain_data_QCd[2:dim(brain_data_QCd)[2]], hgu133a.db_mapping)
dim(brain_data_QCd)
dim(brain_data_QCd_entrez_id)
length(which(duplicated(colnames(brain_data_QCd_entrez_id))))

head(brain_data_QCd_entrez_id[100:110])

##### COLLAPSE MULTIPPLE ENTREZ ID BY SELECTING ONE WITH HIGHEST AVERAGE EXPRESSION ACROSS SAMPLES ######

## create function

select_duplicate_probe_by_top_expr <- function(x) {
  # transpose data frame - keep as dataframe
  x_t<-as.data.frame(t(x))
  # calculate mean expression per probe across samples - create new column - probe mean column
  x_t$probe_mean_expression<-rowMeans(x_t)
  #copy rownames (probe id) to column and truncate
  x_t$trunc_entrez_id<-trunc(as.numeric(as.character(rownames(x_t))))
  # order data frame by truncated probe id and then expression level
  x_t<-x_t[order(x_t$trunc_entrez_id, -x_t$probe_mean_expression), ]
  # remove all duplicate probe id - keep one with highest mean expression
  x_t_unique<-x_t[!duplicated(x_t$trunc_entrez_id),]
  #unique entrez column back to row name
  rownames(x_t_unique)<-x_t_unique$trunc_entrez_id
  #remove unwanted column
  x_t_unique$trunc_entrez_id<-NULL
  #remove unwanted column
  x_t_unique$probe_mean_expression<-NULL
  #transpose dataframe back
  x_unique<-as.data.frame(t(x_t_unique))
  return(x_unique)
}

## apply function to dataframes - check number of probes - check duplicates

brain_data_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(brain_data_QCd_entrez_id)
dim(brain_data_QCd_entrez_id_unique)
length(which(duplicated(colnames(brain_data_QCd_entrez_id_unique))))

head(brain_data_QCd_entrez_id_unique[1:5])

##### ATTACH DIAGNOSIS AND BRAIN REGION #####


# Diagnosis + tissue lookup

head(phenotype_data)

# keep columns - Comment..Sample_characteristics.

phenotype_to_attach<-phenotype_data[6]

head(phenotype_to_attach)

# extract MMSE score  - last 3 characters

phenotype_to_attach$MMSE<-as.numeric(stri_sub(phenotype_to_attach$Comment..Sample_characteristics., -3, -1))

phenotype_to_attach$Comment..Sample_characteristics.<-NULL

head(phenotype_to_attach)

#add Diagnosis

phenotype_to_attach<-merge(phenotype_to_attach, Diagnosis, by="row.names")
rownames(phenotype_to_attach)<- phenotype_to_attach$Row.names
phenotype_to_attach$Row.names<-NULL

# Tissue
phenotype_to_attach$Tissue<-"Prefrontal_Cortex"

# add additional columns as unknown - AGE, GENDER, BRAAK

phenotype_to_attach$BRAAK<-"Unknown"
phenotype_to_attach$Age<-"Unknown"
phenotype_to_attach$APOE<-"Unknown"

# Diagnosis_sub_cat

phenotype_to_attach$Diagnosis_sub_cat<-phenotype_to_attach$Diagnosis
table(phenotype_to_attach$Diagnosis_sub_cat)
table(phenotype_to_attach$Diagnosis)

# change diagnosis name

phenotype_to_attach[phenotype_to_attach$Diagnosis=="Incipient_AD",2]<-"MCI"

# gender

phenotype_to_attach<-merge(phenotype_to_attach, gender_comparison, by="row.names")
rownames(phenotype_to_attach)<- phenotype_to_attach$Row.names
phenotype_to_attach$Row.names<-NULL


#order
phenotype_to_attach<-phenotype_to_attach[,order(names(phenotype_to_attach), decreasing = T)]

# should be 9 colnames - Diagnosis", "Tissue", "Age", "Clinical_Gender", "Predicted_Gender", "MMSE", "BRAAK", "APOE", "Diagnosis_sub_cat"
colnames(phenotype_to_attach)

#merge

brain_data_QCd_entrez_id_unique_pheno<-merge(phenotype_to_attach, brain_data_QCd_entrez_id_unique, by="row.names")

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

rownames(brain_data_QCd_entrez_id_unique_pheno)<-brain_data_QCd_entrez_id_unique_pheno$Row.names

brain_data_QCd_entrez_id_unique_pheno$Row.names<-NULL

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

dim(brain_data_QCd_entrez_id_unique_pheno)
dim(brain_data_QCd_entrez_id_unique)

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(brain_data_QCd_entrez_id_unique_pheno, file="E-GEOD-12685_pre-processed_data2.txt", sep="\t")

##### SAVE PHENOTYPE DATAFRAME #####

setwd(clean_data_dir)

write.table(phenotype_data, file="E-GEOD-12685_phenotype_data2.txt", sep="\t")

##### WRITE LIST OF CONTROL GENES EXPRESSED #####

setwd("/media/hamel/1TB/Projects/Brain_expression/2.Expression_across_brain_regions_in_control_datasets/1.Data")

write(unique(sort(c(control_exprs_M_expressed_probes_list, control_exprs_F_expressed_probes_list))), file="E-GEOD-12685_Prefrontal_cortex.txt")

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="E-GEOD-12685_data_processing2.Rdata")
