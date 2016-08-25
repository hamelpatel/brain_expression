
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                    AMP - MayoeGWAS - PIEPLINE V02                                                      #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Illumina - DASL
# EXPRESSION CHIP 
# NUMBER OF SAMPLES - 
# BRAIN REGIONS - Cerebellum + Temporal Cortex
# 
#
# EPROCESSING DATA USING NEW PIPELINE

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

raw_dir="/media/hamel/Workspace/Projects/Brain_expression/1.Data/AD/AMP_MayoeGWAS/Raw_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MAyoeGWAS/Pre-processing"

setwd(work_dir)

dir.create(paste(work_dir,"pca_plots", sep="/"))

pca_dir=paste(work_dir,"pca_plots", sep="/")

dir.create(paste(work_dir,"clean_data", sep="/"))

clean_data_dir=paste(work_dir,"clean_data", sep="/")

dir.create(paste(work_dir,"sample_network_plots", sep="/"))

sample_network_dir=paste(work_dir,"sample_network_plots", sep="/")

##### LOAD LIBRARIES ####

library(affy)
library(lumi)
library(biomaRt)
library(WGCNA)
library(pamr)
library(limma)
library(sva)
library(illuminaHumanv3.db)
library(ggplot2)
library(reshape)
library(massiR)

##### READ IN PROCESSED DATA ######

setwd(raw_dir)

# read temporal cortex data

Temporal_cortex_exprs<-read.table("MayoEGWAS_arrayExpression_TCX.csv", sep=",", head=T, as.is=T)
Temporal_cortex_pheno<-read.table("MayoEGWAS_arrayExpression_TCX_covariates.csv", sep=",", head=T, as.is=T)

head(Temporal_cortex_exprs)[1:5]
head(Temporal_cortex_pheno)

table(Temporal_cortex_exprs$FID)
anyDuplicated(Temporal_cortex_exprs$IID)
table(Temporal_cortex_pheno$Dxn)

# move IID to rownames and remove 1st 2 columns

rownames(Temporal_cortex_exprs)<-Temporal_cortex_exprs$IID
Temporal_cortex_exprs$FID<-NULL
Temporal_cortex_exprs$IID<-NULL

head(Temporal_cortex_exprs)[1:5]

# transpose dataset for processing

Temporal_cortex_exprs_dataframe<-as.data.frame(t(Temporal_cortex_exprs))
head(Temporal_cortex_exprs_dataframe)[1:5]

# convert Dxn to Diagnosis - 0 == Control and 1 ==AD

Temporal_cortex_pheno[Temporal_cortex_pheno$Dxn=="0",3]<-"Control"
Temporal_cortex_pheno[Temporal_cortex_pheno$Dxn=="1",3]<-"AD"

# convert sex to Gender 0 == Male and 1 == Female

Temporal_cortex_pheno[Temporal_cortex_pheno$Sex=="0",4]<-"Male"
Temporal_cortex_pheno[Temporal_cortex_pheno$Sex=="1",4]<-"Female"

# change column names

colnames(Temporal_cortex_pheno)[c(3,4)]<-c("Diagnosis", "Gender")

head(Temporal_cortex_pheno)
table(Temporal_cortex_pheno$Diagnosis)

# read cerebellum data

Cerebellum_exprs<-read.table("MayoEGWAS_arrayExpression_CBE.csv", sep=",", head=T, as.is=T)
Cerebellum_pheno<-read.table("MayoEGWAS_arrayExpression_CBE_covariates.csv", sep=",", head=T, as.is=T)

head(Cerebellum_exprs)[1:5]
head(Cerebellum_pheno)

table(Cerebellum_exprs$FID)
anyDuplicated(Cerebellum_exprs$IID)
table(Cerebellum_pheno$Dxn)

# move IID to rownames and remove 1st 2 columns

rownames(Cerebellum_exprs)<-Cerebellum_exprs$IID
Cerebellum_exprs$FID<-NULL
Cerebellum_exprs$IID<-NULL

head(Cerebellum_exprs)[1:5]

# transpose dataset for processing

Cerebellum_exprs_dataframe<-as.data.frame(t(Cerebellum_exprs))
head(Cerebellum_exprs_dataframe)[1:5]

# convert Dxn to Diagnosis - 0 == Control and 1 ==AD

Cerebellum_pheno[Cerebellum_pheno$Dxn=="0",3]<-"Control"
Cerebellum_pheno[Cerebellum_pheno$Dxn=="1",3]<-"AD"

# convert sex to Gender 0 == Male and 1 == Female

Cerebellum_pheno[Cerebellum_pheno$Sex=="0",4]<-"Male"
Cerebellum_pheno[Cerebellum_pheno$Sex=="1",4]<-"Female"

# change column names

colnames(Cerebellum_pheno)[c(3,4)]<-c("Diagnosis", "Gender")

head(Cerebellum_pheno)
table(Cerebellum_pheno$Diagnosis)

##### PLOTS OF RAW DATA ####

setwd(work_dir)

pdf(file="processed_data_boxplot.pdf")
boxplot(Temporal_cortex_exprs_dataframe)
boxplot(Cerebellum_exprs_dataframe)
dev.off()

pdf(file="processed_data_density_plot.pdf")
plotDensity(as.matrix(Temporal_cortex_exprs_dataframe),logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples" )
plotDensity(as.matrix(Cerebellum_exprs_dataframe), logMode=F, addLegend=F, main="Density plot of Cerebellum Samples")
dev.off()

##### EXTRACT/STANDARDISE PHENOTYPE INFO ####

head(Cerebellum_pheno)
head(Temporal_cortex_pheno)

table(Cerebellum_pheno$Diagnosis)
table(Temporal_cortex_pheno$Diagnosis)

# move IID to rownames, and remove FID

rownames(Temporal_cortex_pheno)<-Temporal_cortex_pheno$IID
Temporal_cortex_pheno$FID<-NULL
Temporal_cortex_pheno$IID<-NULL

rownames(Cerebellum_pheno)<-Cerebellum_pheno$IID
Cerebellum_pheno$FID<-NULL
Cerebellum_pheno$IID<-NULL

head(Cerebellum_pheno)
head(Temporal_cortex_pheno)

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(Cerebellum_pheno$IID))
anyDuplicated(rownames(Temporal_cortex_pheno$IID))

anyDuplicated(colnames(Temporal_cortex_exprs_dataframe))
anyDuplicated(colnames(Cerebellum_exprs_dataframe))

##### GENDER CHECK ##### 

head(Cerebellum_pheno)
names(Cerebellum_pheno)

# gender_information - Cerebellum
Cerebellum_gender_info<-Cerebellum_pheno[2]
colnames(Cerebellum_gender_info)<-"Gender"
table(Cerebellum_gender_info$Gender)
head(Cerebellum_gender_info)

# gender_information - Temporal_cortex
Temporal_cortex_gender_info<-Temporal_cortex_pheno[2]
colnames(Temporal_cortex_gender_info)<-"Gender"
table(Temporal_cortex_gender_info$Gender)
head(Temporal_cortex_gender_info)

#standardise gender

Cerebellum_gender_info[Cerebellum_gender_info$Gender=="Female",1]<-"female"
Cerebellum_gender_info[Cerebellum_gender_info$Gender=="Male",1]<-"male"
Temporal_cortex_gender_info[Temporal_cortex_gender_info$Gender=="Female",1]<-"female"
Temporal_cortex_gender_info[Temporal_cortex_gender_info$Gender=="Male",1]<-"male"

#separae male/female IDs
Cerebellum_female_samples<-subset(Cerebellum_gender_info, Gender=="female")
Cerebellum_male_samples<-subset(Cerebellum_gender_info, Gender=="male")
Temporal_cortex_female_samples<-subset(Temporal_cortex_gender_info, Gender=="female")
Temporal_cortex_male_samples<-subset(Temporal_cortex_gender_info, Gender=="male")

head(Cerebellum_female_samples)
head(Cerebellum_male_samples)
head(Temporal_cortex_female_samples)
head(Temporal_cortex_male_samples)

# get Y choromosome genes
data(y.probes)
names(y.probes)

y_chromo_probes <- data.frame(y.probes["illumina_humanwg_6_v2"])

y_chromo_probes

# extract Y chromosome genes from dataset
Cerebellum.eset.select.out <- massi_select(Cerebellum_exprs_dataframe, y_chromo_probes)
Temporal_cortex.eset.select.out <- massi_select(Temporal_cortex_exprs_dataframe, y_chromo_probes)

massi_y_plot(Cerebellum.eset.select.out)
massi_cluster_plot(Cerebellum.eset.select.out)
massi_y_plot(Temporal_cortex.eset.select.out)
massi_cluster_plot(Temporal_cortex.eset.select.out)

# run gender predict
Cerebellum.eset.results <- massi_cluster(Cerebellum.eset.select.out)
Temporal_cortex.eset.results <- massi_cluster(Temporal_cortex.eset.select.out)

#extract gender prediction
Cerebellum_predicted_gender<-(Cerebellum.eset.results$massi.results)[c(1,5)]
rownames(Cerebellum_predicted_gender)<-Cerebellum_predicted_gender$ID
Cerebellum_predicted_gender$ID<-NULL
colnames(Cerebellum_predicted_gender)<-"Predicted_Gender"

Temporal_cortex_predicted_gender<-(Temporal_cortex.eset.results$massi.results)[c(1,5)]
rownames(Temporal_cortex_predicted_gender)<-Temporal_cortex_predicted_gender$ID
Temporal_cortex_predicted_gender$ID<-NULL
colnames(Temporal_cortex_predicted_gender)<-"Predicted_Gender"


#merge
Cerebellum_gender_comparison<-merge(Cerebellum_gender_info, Cerebellum_predicted_gender, by="row.names")
rownames(Cerebellum_gender_comparison)<-Cerebellum_gender_comparison$Row.names
Cerebellum_gender_comparison$Row.names<-NULL
colnames(Cerebellum_gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
head(Cerebellum_gender_comparison)

Temporal_cortex_gender_comparison<-merge(Temporal_cortex_gender_info, Temporal_cortex_predicted_gender, by="row.names")
rownames(Temporal_cortex_gender_comparison)<-Temporal_cortex_gender_comparison$Row.names
Temporal_cortex_gender_comparison$Row.names<-NULL
colnames(Temporal_cortex_gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
head(Temporal_cortex_gender_comparison)

#differences

Cerebellum_Gender_missmatch<-Cerebellum_gender_comparison[which(Cerebellum_gender_comparison$Clinical_Gender!=Cerebellum_gender_comparison$Predicted_Gender),]
Cerebellum_Gender_missmatch

Temporal_cortex_Gender_missmatch<-Temporal_cortex_gender_comparison[which(Temporal_cortex_gender_comparison$Clinical_Gender!=Temporal_cortex_gender_comparison$Predicted_Gender),]
Temporal_cortex_Gender_missmatch

##### PROBE ID DETECTION #####

# separate case control 

Cerebellum_case_ID<-rownames(Cerebellum_pheno[Cerebellum_pheno$Diagnosis=="AD",])
Cerebellum_control_ID<-rownames(Cerebellum_pheno[Cerebellum_pheno$Diagnosis=="Control",])
Temporal_cortex_case_ID<-rownames(Temporal_cortex_pheno[Temporal_cortex_pheno$Diagnosis=="AD",])
Temporal_cortex_control_ID<-rownames(Temporal_cortex_pheno[Temporal_cortex_pheno$Diagnosis=="Control",])

Cerebellum_case_exprs<-Cerebellum_exprs_dataframe[,colnames(Cerebellum_exprs_dataframe)%in%Cerebellum_case_ID]
Cerebellum_control_exprs<-Cerebellum_exprs_dataframe[,colnames(Cerebellum_exprs_dataframe)%in%Cerebellum_control_ID]
Temporal_cortex_case_exprs<-Temporal_cortex_exprs_dataframe[,colnames(Temporal_cortex_exprs_dataframe)%in%Temporal_cortex_case_ID]
Temporal_cortex_control_exprs<-Temporal_cortex_exprs_dataframe[,colnames(Temporal_cortex_exprs_dataframe)%in%Temporal_cortex_control_ID]

head(Cerebellum_case_exprs[1:5])
dim(Cerebellum_case_exprs)
head(Cerebellum_control_exprs[1:5])
dim(Cerebellum_control_exprs)
head(Temporal_cortex_case_exprs[1:5])
dim(Temporal_cortex_case_exprs)
head(Temporal_cortex_control_exprs[1:5])
dim(Temporal_cortex_control_exprs)


#separate by gender

Cerebellum_control_exprs_F<-Cerebellum_control_exprs[colnames(Cerebellum_control_exprs)%in%rownames(Cerebellum_female_samples)]
Temporal_cortex_control_exprs_F<-Temporal_cortex_control_exprs[colnames(Temporal_cortex_control_exprs)%in%rownames(Temporal_cortex_female_samples)]

Cerebellum_case_exprs_F<-Cerebellum_case_exprs[colnames(Cerebellum_case_exprs)%in%rownames(Cerebellum_female_samples)]
Temporal_cortex_case_exprs_F<-Temporal_cortex_case_exprs[colnames(Temporal_cortex_case_exprs)%in%rownames(Temporal_cortex_female_samples)]

Cerebellum_control_exprs_M<-Cerebellum_control_exprs[colnames(Cerebellum_control_exprs)%in%rownames(Cerebellum_male_samples)]
Temporal_cortex_control_exprs_M<-Temporal_cortex_control_exprs[colnames(Temporal_cortex_control_exprs)%in%rownames(Temporal_cortex_male_samples)]

Cerebellum_case_exprs_M<-Cerebellum_case_exprs[colnames(Cerebellum_case_exprs)%in%rownames(Cerebellum_male_samples)]
Temporal_cortex_case_exprs_M<-Temporal_cortex_case_exprs[colnames(Temporal_cortex_case_exprs)%in%rownames(Temporal_cortex_male_samples)]

# calculate 90th percentile for each sample in each group - quantile normalisedd - dataset already 25th percentile removed - adjusting to 90th percentile - using 86%

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

#apply 

Cerebellum_case_expressed_probes_list_F<-extract_good_probe_list(Cerebellum_case_exprs_F, 0.7, 0.8)
length(Cerebellum_case_expressed_probes_list_F)

Temporal_cortex_case_expressed_probes_list_F<-extract_good_probe_list(Temporal_cortex_case_exprs_F, 0.7, 0.8)
length(Temporal_cortex_case_expressed_probes_list_F)

Cerebellum_control_expressed_probes_list_F<-extract_good_probe_list(Cerebellum_control_exprs_F, 0.7, 0.8)
length(Cerebellum_control_expressed_probes_list_F)

Temporal_cortex_control_expressed_probes_list_F<-extract_good_probe_list(Temporal_cortex_control_exprs_F, 0.7, 0.8)
length(Temporal_cortex_control_expressed_probes_list_F)

Cerebellum_case_expressed_probes_list_M<-extract_good_probe_list(Cerebellum_case_exprs_M, 0.7, 0.8)
length(Cerebellum_case_expressed_probes_list_M)

Temporal_cortex_case_expressed_probes_list_M<-extract_good_probe_list(Temporal_cortex_case_exprs_M, 0.7, 0.8)
length(Temporal_cortex_case_expressed_probes_list_M)

Cerebellum_control_expressed_probes_list_M<-extract_good_probe_list(Cerebellum_control_exprs_M, 0.7, 0.8)
length(Cerebellum_control_expressed_probes_list_M)

Temporal_cortex_control_expressed_probes_list_M<-extract_good_probe_list(Temporal_cortex_control_exprs_M, 0.7, 0.8)
length(Temporal_cortex_control_expressed_probes_list_M)

# merge list of good probes from both case + control + and all brain regions, sort and keep unique values

Cerebellum_good_probe_list<-unique(sort(c(Cerebellum_control_expressed_probes_list_F,
                                          Cerebellum_case_expressed_probes_list_F,
                                          Cerebellum_control_expressed_probes_list_M,
                                          Cerebellum_case_expressed_probes_list_M)))

Temporal_cortex_good_probe_list<-unique(sort(c(Temporal_cortex_control_expressed_probes_list_F,
                                               Temporal_cortex_case_expressed_probes_list_F,
                                               Temporal_cortex_control_expressed_probes_list_M,
                                               Temporal_cortex_case_expressed_probes_list_M)))
length(Cerebellum_good_probe_list)
length(Temporal_cortex_good_probe_list)

# extract good probes from dataset

Cerebellum_case_exprs_good_probes<-Cerebellum_case_exprs[rownames(Cerebellum_case_exprs)%in%Cerebellum_good_probe_list,]

Temporal_cortex_case_exprs_good_probes<-Temporal_cortex_case_exprs[rownames(Temporal_cortex_case_exprs)%in%Temporal_cortex_good_probe_list,]

Cerebellum_control_exprs_good_probes<-Cerebellum_control_exprs[rownames(Cerebellum_control_exprs)%in%Cerebellum_good_probe_list,]

Temporal_cortex_control_exprs_good_probes<-Temporal_cortex_control_exprs[rownames(Temporal_cortex_control_exprs)%in%Temporal_cortex_good_probe_list,]

Cerebellum_exprs_good_probes<-as.data.frame(t(Cerebellum_exprs[,colnames(Cerebellum_exprs)%in%Cerebellum_good_probe_list]))
Temporal_cortex_exprs_good_probes<-as.data.frame(t(Temporal_cortex_exprs[,colnames(Temporal_cortex_exprs)%in%Temporal_cortex_good_probe_list]))

# check dataframe 

head(Cerebellum_case_exprs_good_probes)[1:5]
dim(Cerebellum_case_exprs_good_probes)

head(Temporal_cortex_case_exprs_good_probes)[1:5]
dim(Temporal_cortex_case_exprs_good_probes)

head(Cerebellum_control_exprs_good_probes)[1:5]
dim(Cerebellum_control_exprs_good_probes)

head(Temporal_cortex_control_exprs_good_probes)[1:5]
dim(Temporal_cortex_control_exprs_good_probes)

head(Cerebellum_exprs_good_probes)[1:5]
dim(Cerebellum_exprs_good_probes)

head(Temporal_cortex_exprs_good_probes)[1:5]
dim(Temporal_cortex_exprs_good_probes)

##### GENDER SPECIFIC PROBE PLOTS #####

# using dataframe before probe removal

# get gene symbol list for chip
Gene_symbols_probes <- mappedkeys(illuminaHumanv3SYMBOL)

# Convert to a list
Gene_symbols <- as.data.frame(illuminaHumanv3SYMBOL[Gene_symbols_probes])

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

# merge all genes to check

gene_list<-rbind(XIST_probe_ID,
                 PRKY_probe_ID,
                 NLGN4Y_probe_ID,
                 TMSB4Y_probe_ID,
                 USP9Y_probe_ID,
                 UTY_probe_ID)

gene_list

#create function to plot
plot_gender_specific_genes<-function(Expression_table, gender_info, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[rownames(Expression_table) %in% genes_to_extract$probe_id,]))
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

pdf("Gender_specific_gene_plot_and_detectable_cut_off_threshold_used_09.pdf")
plot_gender_specific_genes(Cerebellum_control_exprs, Cerebellum_gender_comparison[1], gene_list, 0.7, "Cerebellum_control_Clinical_Gender")
plot_gender_specific_genes(Cerebellum_control_exprs, Cerebellum_gender_comparison[2], gene_list, 0.7, "Cerebellum_control_Predicted_Gender")
plot_gender_specific_genes(Temporal_cortex_control_exprs, Temporal_cortex_gender_comparison[1], gene_list, 0.7, "Temporal_cortex_control_Clinical_Gender")
plot_gender_specific_genes(Temporal_cortex_control_exprs, Temporal_cortex_gender_comparison[2], gene_list, 0.7, "Temporal_cortex_control_Predicted_Gender")
plot_gender_specific_genes(Cerebellum_case_exprs, Cerebellum_gender_comparison[1], gene_list, 0.7, "Cerebellum_case_Clinical_Gender")
plot_gender_specific_genes(Cerebellum_case_exprs, Cerebellum_gender_comparison[2], gene_list, 0.7, "Cerebellum_case_Predicted_Gender")
plot_gender_specific_genes(Temporal_cortex_case_exprs, Temporal_cortex_gender_comparison[1], gene_list, 0.7, "Temporal_cortex_case_Clinical_Gender")
plot_gender_specific_genes(Temporal_cortex_case_exprs, Temporal_cortex_gender_comparison[2], gene_list, 0.7, "Temporal_cortex_case_Predicted_Gender")
dev.off()

##### PCA ####

# TEMPORAL CORTEX

# calculate pca

Temporal_cortex_pca<-prcomp(Temporal_cortex_exprs_good_probes)

# plot variance
plot(Temporal_cortex_pca, type="l")

# sumary of pca
Temporal_cortex_summary_pca<-summary(Temporal_cortex_pca)

# pca importance
Temporal_cortex_pca_importance_var_exp<-Temporal_cortex_summary_pca$importance[2,]
Temporal_cortex_pca_importance_var_exp_cum<-Temporal_cortex_summary_pca$importance[3,]

par(mfrow=c(1,2))

#plot variance explained
plot(Temporal_cortex_pca_importance_var_exp, ylab="PCA Proportion of Variance Explained", type="b", col="blue")

#plot variance explained cumalative
plot(Temporal_cortex_pca_importance_var_exp_cum, ylab="PCA Cumulative Proportion of Variance Explained", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)

par(dev.off)
dev.off()

# pca matrix plot

#color

head(Temporal_cortex_pheno)

Temporal_cortex_Diagnosis<-Temporal_cortex_pheno[1]
head(Temporal_cortex_Diagnosis)
table(Temporal_cortex_Diagnosis)

# order of samples in expression data
Temporal_cortex_Diagnosis_temp<-colnames(Temporal_cortex_exprs_good_probes)

# match order
Temporal_cortex_Diagnosis_pca<-Temporal_cortex_Diagnosis[match(Temporal_cortex_Diagnosis_temp, rownames(Temporal_cortex_Diagnosis)),]
head(Temporal_cortex_Diagnosis_pca)

# assign color to group
Temporal_cortex_Diagnosis_pca_color<-labels2colors(as.character(Temporal_cortex_Diagnosis_pca))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Temporal_cortex_Diagnosis_pca_color)
legend('bottomright', unique(Temporal_cortex_Diagnosis_pca), fill=unique(Temporal_cortex_Diagnosis_pca_color))

# CEREBEULLUM

# calculate pca

Cerebellum_pca<-prcomp(Cerebellum_exprs_good_probes)

# plot variance
plot(Cerebellum_pca, type="l")

# sumary of pca
Cerebellum_summary_pca<-summary(Cerebellum_pca)

# pca importance
Cerebellum_pca_importance_var_exp<-Cerebellum_summary_pca$importance[2,]
Cerebellum_pca_importance_var_exp_cum<-Cerebellum_summary_pca$importance[3,]

par(mfrow=c(1,2))

#plot variance explained
plot(Cerebellum_pca_importance_var_exp, ylab="PCA Proportion of Variance Explained", type="b", col="blue")

#plot variance explained cumalative
plot(Cerebellum_pca_importance_var_exp_cum, ylab="PCA Cumulative Proportion of Variance Explained", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)

par(dev.off)
dev.off()

# pca matrix plot

#color

head(Cerebellum_pheno)

Cerebellum_Diagnosis<-Cerebellum_pheno[1]
head(Cerebellum_Diagnosis)
table(Cerebellum_Diagnosis)

# order of samples in expression data
Cerebellum_Diagnosis_temp<-colnames(Cerebellum_exprs_good_probes)

# match order
Cerebellum_Diagnosis_pca<-Cerebellum_Diagnosis[match(Cerebellum_Diagnosis_temp, rownames(Cerebellum_Diagnosis)),]
head(Cerebellum_Diagnosis_pca)

# assign color to group
Cerebellum_Diagnosis_pca_color<-labels2colors(as.character(Cerebellum_Diagnosis_pca))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Cerebellum PCA plot coloured by Disease Status",col="black", pch=21,bg=Cerebellum_Diagnosis_pca_color)
legend('bottomright', unique(Cerebellum_Diagnosis_pca), fill=unique(Cerebellum_Diagnosis_pca_color))

#plot to pdf

setwd(pca_dir)

pdf("PCA_plot_prior_to_sample_removal.pdf")
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Temporal_cortex_Diagnosis_pca_color)
legend('bottomleft', unique(Temporal_cortex_Diagnosis_pca), fill=unique(Temporal_cortex_Diagnosis_pca_color))
plot(Cerebellum_pca$rotation[,1:2], main="Cerebellum PCA plot coloured by Disease Status",col="black", pch=21,bg=Cerebellum_Diagnosis_pca_color)
legend('bottomleft', unique(Cerebellum_Diagnosis_pca), fill=unique(Cerebellum_Diagnosis_pca_color))
dev.off()

##### PCA - BATCH EFFECT EXPLORATION ####

setwd(pca_dir)

# match order
Temporal_cortex_test<-Temporal_cortex_pheno[match(Temporal_cortex_Diagnosis_temp, rownames(Temporal_cortex_pheno)),]
head(Temporal_cortex_test)
rownames(Temporal_cortex_test)==colnames(Temporal_cortex_exprs_good_probes)

pdf("Temporal_Cortex_batch_effect_PCA_plot.pdf")

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Temporal_cortex_Diagnosis_pca_color)
legend('bottomright', unique(Temporal_cortex_Diagnosis_pca), fill=unique(Temporal_cortex_Diagnosis_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,5]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_test[,5])), fill=unique(Temporal_cortex_test_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,6]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_test[,6])), fill=unique(Temporal_cortex_test_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,7]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_test[,7])), fill=unique(Temporal_cortex_test_pca_color))

# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,8]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_test[,8])), fill=unique(Temporal_cortex_test_pca_color))

# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,2]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by gender",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_test[,2])), fill=unique(Temporal_cortex_test_pca_color))


dev.off()


#
#
# Cerebellum
#
#


# match order
Cerebellum_test<-Cerebellum_pheno[match(Cerebellum_Diagnosis_temp, rownames(Cerebellum_pheno)),]
head(Cerebellum_test)
rownames(Cerebellum_test)==colnames(Cerebellum_exprs_good_probes)

pdf("Cerebellum_batch_effect_PCA_plot.pdf")

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Cerebellum_Diagnosis_pca_color)
legend('bottomright', unique(Cerebellum_Diagnosis_pca), fill=unique(Cerebellum_Diagnosis_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,5]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate0 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,5])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,6]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate1 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,6])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,7]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,7])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,8]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,8])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,9]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,9])), fill=unique(Cerebellum_test_pca_color))

# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,2]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Gender",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomright', as.character(unique(Cerebellum_test[,2])), fill=unique(Cerebellum_test_pca_color))

dev.off()

##### SVA ####

# add diagnosis in

Cerebellum_pheno
Temporal_cortex_pheno

Cerebellum_exprs_good_probes<-merge (Cerebellum_pheno[1], t(Cerebellum_exprs_good_probes), by="row.names", all.y=T)
rownames(Cerebellum_exprs_good_probes)<-Cerebellum_exprs_good_probes$Row.names
Cerebellum_exprs_good_probes$Row.names<-NULL
head(Cerebellum_exprs_good_probes)[1:5]

Temporal_cortex_exprs_good_probes<-merge (Temporal_cortex_pheno[1], t(Temporal_cortex_exprs_good_probes), by="row.names", all.y=T)
rownames(Temporal_cortex_exprs_good_probes)<-Temporal_cortex_exprs_good_probes$Row.names
Temporal_cortex_exprs_good_probes$Row.names<-NULL
head(Temporal_cortex_exprs_good_probes)[1:5]

# add gender in 

Cerebellum_exprs_good_probes<-merge(Cerebellum_gender_comparison[1], Cerebellum_exprs_good_probes, by="row.names")
rownames(Cerebellum_exprs_good_probes)<-Cerebellum_exprs_good_probes$Row.names
Cerebellum_exprs_good_probes$Row.names<-NULL

Temporal_cortex_exprs_good_probes<-merge(Temporal_cortex_gender_comparison[1], Temporal_cortex_exprs_good_probes, by="row.names")
rownames(Temporal_cortex_exprs_good_probes)<-Temporal_cortex_exprs_good_probes$Row.names
Temporal_cortex_exprs_good_probes$Row.names<-NULL

head(Cerebellum_exprs_good_probes)[1:5]
head(Temporal_cortex_exprs_good_probes)[1:5]

# create sva function

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis+Clinical_Gender, data=dataset_pheno)
  # check number of SV in data
  num.sv(dataset_exprs, mod, method="leek")
}

# apply function
check_SV_in_data(Cerebellum_exprs_good_probes)
check_SV_in_data(Temporal_cortex_exprs_good_probes)

# create function to adjust for SV

adjust_for_sva<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_sva_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_sva_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis+Clinical_Gender, data=dataset_sva_pheno)
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

Temporal_cortex_exprs_good_probes_sva_adjusted<-adjust_for_sva(Temporal_cortex_exprs_good_probes)

# change name

Cerebellum_exprs_good_probes_sva_adjusted<-Cerebellum_exprs_good_probes

# check SV again

check_SV_in_data(Cerebellum_exprs_good_probes_sva_adjusted)
check_SV_in_data(Temporal_cortex_exprs_good_probes_sva_adjusted)

# remove gender

Temporal_cortex_exprs_good_probes_sva_adjusted$Clinical_Gender<-NULL
Cerebellum_exprs_good_probes_sva_adjusted$Clinical_Gender<-NULL

# separate case and control

Cerebellum_good_probes_sva_adjusted_case<-Cerebellum_exprs_good_probes_sva_adjusted[Cerebellum_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Cerebellum_good_probes_sva_adjusted_control<-Cerebellum_exprs_good_probes_sva_adjusted[Cerebellum_exprs_good_probes_sva_adjusted$Diagnosis=="Control",]

dim(Cerebellum_good_probes_sva_adjusted_case)
dim(Cerebellum_good_probes_sva_adjusted_control)

Temporal_cortex_good_probes_sva_adjusted_case<-Temporal_cortex_exprs_good_probes_sva_adjusted[Temporal_cortex_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Temporal_cortex_good_probes_sva_adjusted_control<-Temporal_cortex_exprs_good_probes_sva_adjusted[Temporal_cortex_exprs_good_probes_sva_adjusted$Diagnosis=="Control",]

dim(Temporal_cortex_good_probes_sva_adjusted_case)
dim(Temporal_cortex_good_probes_sva_adjusted_control)

##### PLOTS AFTER SVA #####

# plot density plot on clean data
plotDensity(t(as.matrix(Temporal_cortex_exprs_good_probes_sva_adjusted[2:dim(Temporal_cortex_exprs_good_probes_sva_adjusted)[1]])),logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )

# plot PCA plot on clean data

Temporal_cortex_SVA_pca<-prcomp(Temporal_cortex_exprs_good_probes_sva_adjusted[2:dim(Temporal_cortex_exprs_good_probes_sva_adjusted)[1]])

# plot variance
plot(Temporal_cortex_SVA_pca, type="l")

#color

head(Temporal_cortex_pheno)

Temporal_cortex_SVA_Diagnosis<-Temporal_cortex_pheno[1]
head(Temporal_cortex_SVA_Diagnosis)
table(Temporal_cortex_SVA_Diagnosis)

# order of samples in expression data
Temporal_cortex_SVA_Diagnosis_temp<-rownames(Temporal_cortex_exprs_good_probes_sva_adjusted)

# match order
Temporal_cortex_SVA_test<-Temporal_cortex_pheno[match(Temporal_cortex_SVA_Diagnosis_temp, rownames(Temporal_cortex_pheno)),]
head(Temporal_cortex_SVA_test)
rownames(Temporal_cortex_SVA_test)==rownames(Temporal_cortex_exprs_good_probes_sva_adjusted)

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,1]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,1])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,5]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,5])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,6]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,6])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,7]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,7])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,8]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,8])), fill=unique(Temporal_cortex_SVA_test_pca_color))


#plot to pdf
setwd(pca_dir)

pdf("Temporal_cortex_plots_after_SVA_adjustment.pdf")
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,1]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,1])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,5])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,6]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,6])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,7]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,7])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,8]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,8])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,2]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by gender",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomright', as.character(unique(Temporal_cortex_SVA_test[,2])), fill=unique(Temporal_cortex_SVA_test_pca_color))
dev.off()

#plot density plot
setwd(work_dir)

pdf("Temporal_Cortex_Density_plot_after_SVA.pdf")
plotDensity(t(as.matrix(Temporal_cortex_exprs_good_probes_sva_adjusted[2:dim(Temporal_cortex_exprs_good_probes_sva_adjusted)[1]])),logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )
dev.off()

###
#
# Cerebellum
#
###

# plot density plot on clean data
plotDensity(t(as.matrix(Cerebellum_exprs_good_probes_sva_adjusted[2:dim(Cerebellum_exprs_good_probes_sva_adjusted)[1]])),logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )

# plot PCA plot on clean data

Cerebellum_SVA_pca<-prcomp(Cerebellum_exprs_good_probes_sva_adjusted[2:dim(Cerebellum_exprs_good_probes_sva_adjusted)[1]])

# plot variance
plot(Cerebellum_SVA_pca, type="l")

#color

head(Cerebellum_pheno)

Cerebellum_SVA_Diagnosis<-Cerebellum_pheno[1]
head(Cerebellum_SVA_Diagnosis)
table(Cerebellum_SVA_Diagnosis)

# order of samples in expression data
Cerebellum_SVA_Diagnosis_temp<-rownames(Cerebellum_exprs_good_probes_sva_adjusted)

# match order
Cerebellum_SVA_test<-Cerebellum_pheno[match(Cerebellum_SVA_Diagnosis_temp, rownames(Cerebellum_pheno)),]
head(Cerebellum_SVA_test)
rownames(Cerebellum_SVA_test)==rownames(Cerebellum_exprs_good_probes_sva_adjusted)

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,1]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,1])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,5]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate0 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,5])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,6]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate1 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,6])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,7]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,7])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,8]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,8])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,9]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,9])), fill=unique(Cerebellum_SVA_test_pca_color))

# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,2]))

# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by gender",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,2])), fill=unique(Cerebellum_SVA_test_pca_color))


#plot to pdf
setwd(pca_dir)
pdf("Cerebellum_plots_after_SVA_adjustment.pdf")
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,1]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,1])), fill=unique(Cerebellum_SVA_test_pca_color))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate0 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,5])), fill=unique(Cerebellum_SVA_test_pca_color))
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,6]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate1 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,6])), fill=unique(Cerebellum_SVA_test_pca_color))
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,7]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,7])), fill=unique(Cerebellum_SVA_test_pca_color))
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,8]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,8])), fill=unique(Cerebellum_SVA_test_pca_color))
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,9]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,9])), fill=unique(Cerebellum_SVA_test_pca_color))
# assign color to group
Cerebellum_SVA_test_pca_color<-labels2colors(as.character(Cerebellum_SVA_test[,2]))
# pca plot - color by disease - case/control
plot(Cerebellum_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by gender",col="black", pch=21,bg=Cerebellum_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_SVA_test[,2])), fill=unique(Cerebellum_SVA_test_pca_color))
dev.off()

#plot density plot
setwd(work_dir)
pdf("Cerebellum_Density_plot_after_SVA.pdf")
Cerebellum_SVA_pca<-prcomp(Cerebellum_exprs_good_probes_sva_adjusted[2:dim(Cerebellum_exprs_good_probes_sva_adjusted)[1]])
dev.off()

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

# run sample network on Cerebellum - rho is low, however moving Z.K to -2 removes 28 samples, and rho is still low.

Cerebellum_case_exprs_good_probes_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case, -3)
Cerebellum_control_exprs_good_probes_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control, -3)

# run sample network on Temporal_cortex

Temporal_cortex_case_exprs_good_probes_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_case, -3)
Temporal_cortex_control_exprs_good_probes_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_control, -3)

#separate by gender - additional check to see if samples are influenced by gender

Cerebellum_good_probes_sva_adjusted_control_F<-Cerebellum_good_probes_sva_adjusted_control[rownames(Cerebellum_good_probes_sva_adjusted_control)%in%rownames(Cerebellum_female_samples),]
Temporal_cortex_good_probes_sva_adjusted_control_F<-Temporal_cortex_good_probes_sva_adjusted_control[rownames(Temporal_cortex_good_probes_sva_adjusted_control)%in%rownames(Temporal_cortex_female_samples),]

Cerebellum_good_probes_sva_adjusted_case_F<-Cerebellum_good_probes_sva_adjusted_case[rownames(Cerebellum_good_probes_sva_adjusted_case)%in%rownames(Cerebellum_female_samples),]
Temporal_cortex_good_probes_sva_adjusted_case_F<-Temporal_cortex_good_probes_sva_adjusted_case[rownames(Temporal_cortex_good_probes_sva_adjusted_case)%in%rownames(Temporal_cortex_female_samples),]

Cerebellum_good_probes_sva_adjusted_control_M<-Cerebellum_good_probes_sva_adjusted_control[rownames(Cerebellum_good_probes_sva_adjusted_control)%in%rownames(Cerebellum_male_samples),]
Temporal_cortex_good_probes_sva_adjusted_control_M<-Temporal_cortex_good_probes_sva_adjusted_control[rownames(Temporal_cortex_good_probes_sva_adjusted_control)%in%rownames(Temporal_cortex_male_samples),]

Cerebellum_good_probes_sva_adjusted_case_M<-Cerebellum_good_probes_sva_adjusted_case[rownames(Cerebellum_good_probes_sva_adjusted_case)%in%rownames(Cerebellum_male_samples),]
Temporal_cortex_good_probes_sva_adjusted_case_M<-Temporal_cortex_good_probes_sva_adjusted_case[rownames(Temporal_cortex_good_probes_sva_adjusted_case)%in%rownames(Temporal_cortex_male_samples),]

# run sample network analysis

Cerebellum_good_probes_sva_adjusted_control_F_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control_F, -3)
Cerebellum_good_probes_sva_adjusted_case_F_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case_F, -3)
Cerebellum_good_probes_sva_adjusted_control_M_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control_M, -3)
Cerebellum_good_probes_sva_adjusted_case_M_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case_M, -3)

Temporal_cortex_good_probes_sva_adjusted_control_F_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_control_F, -3)
Temporal_cortex_good_probes_sva_adjusted_case_F_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_case_F, -3)
Temporal_cortex_good_probes_sva_adjusted_control_M_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_control_M, -3)
Temporal_cortex_good_probes_sva_adjusted_case_M_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_case_M, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("Cerebellum_CA1_case_sample_network_analysis.pdf")
Cerebellum_case_exprs_good_probes_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Cerebellum_CA1_control_sample_network_analysis.pdf")
Cerebellum_control_exprs_good_probes_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Temporal_cortex_case_sample_network_analysis.pdf")
Temporal_cortex_case_exprs_good_probes_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Temporal_cortex_control_sample_network_analysis.pdf")
Temporal_cortex_control_exprs_good_probes_QC<-run_sample_network_plot(Temporal_cortex_good_probes_sva_adjusted_control, -3)
dev.off()

##### CREATE QC'd DATASET #####

Temporal_cortex_good_sample_ID<-c(rownames(Temporal_cortex_case_exprs_good_probes_QC),
                                  rownames(Temporal_cortex_control_exprs_good_probes_QC))

Cerebellum_good_sample_ID<-c(rownames(Cerebellum_case_exprs_good_probes_QC),
                             rownames(Cerebellum_control_exprs_good_probes_QC))

Temporal_cortex_QCd<-Temporal_cortex_exprs_good_probes_sva_adjusted[rownames(Temporal_cortex_exprs_good_probes_sva_adjusted)%in%Temporal_cortex_good_sample_ID,]
dim(Temporal_cortex_exprs_good_probes_sva_adjusted)
dim(Temporal_cortex_QCd)

Cerebellum_QCd<-Cerebellum_exprs_good_probes_sva_adjusted[rownames(Cerebellum_exprs_good_probes_sva_adjusted)%in%Cerebellum_good_sample_ID,]
dim(Cerebellum_exprs_good_probes_sva_adjusted)
dim(Cerebellum_QCd)

# remove gender and diagnosis

Temporal_cortex_QCd$Clinical_Gender<-NULL
Temporal_cortex_QCd$Diagnosis<-NULL

Cerebellum_QCd$Clinical_Gender<-NULL
Cerebellum_QCd$Diagnosis<-NULL

##### CONVERT PROBE ID TO ENTREZ ID #####

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using hgu133a.db
mapped_probes <- mappedkeys(illuminaHumanv3ENTREZID)

# Convert to a list
illuminaHumanv3.db_mapping <- as.data.frame(illuminaHumanv3ENTREZID[mapped_probes])
# arrange order of column by entrezgene probe_id
illuminaHumanv3.db_mapping<-illuminaHumanv3.db_mapping[c(2,1)]
colnames(illuminaHumanv3.db_mapping)[1]<-"entrezgene"

head(illuminaHumanv3.db_mapping)
dim(illuminaHumanv3.db_mapping)

#check any duplicated probe IDs
anyDuplicated(illuminaHumanv3.db_mapping$probe_id)

#check any dupliacted entrezgene IDs
anyDuplicated(illuminaHumanv3.db_mapping$entrezgene)

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset
  expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}

# using hgu133plus2.db

Temporal_cortex_QCd_entrez_id<-convert_probe_id_to_entrez_id(Temporal_cortex_QCd, illuminaHumanv3.db_mapping)
dim(Temporal_cortex_QCd)
dim(Temporal_cortex_QCd_entrez_id)
length(which(duplicated(colnames(Temporal_cortex_QCd_entrez_id))))

head(Temporal_cortex_QCd_entrez_id)[1:5]

# apply function to Cerebellum

Cerebellum_QCd_entrez_id<-convert_probe_id_to_entrez_id(Cerebellum_QCd, illuminaHumanv3.db_mapping)
dim(Cerebellum_QCd)
dim(Cerebellum_QCd_entrez_id)
length(which(duplicated(colnames(Cerebellum_QCd_entrez_id))))

head(Cerebellum_QCd_entrez_id[100:110])

##### COLLAPSE MULTIPLE ENTREZ ID BY SELECTING ONE WITH HIGHEST AVERAGE EXPRESSION ACROSS SAMPLES ######

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

## apply function to Temporal Cortex - check number of probes - check duplicates

Temporal_cortex_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(Temporal_cortex_QCd_entrez_id)
dim(Temporal_cortex_QCd_entrez_id_unique)
length(which(duplicated(colnames(Temporal_cortex_QCd_entrez_id_unique))))

head(Temporal_cortex_QCd_entrez_id_unique[1:5])

## apply function to Cerebellum- check number of probes - check duplicates

Cerebellum_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(Cerebellum_QCd_entrez_id)
dim(Cerebellum_QCd_entrez_id_unique)
length(which(duplicated(colnames(Cerebellum_QCd_entrez_id_unique))))

head(Cerebellum_QCd_entrez_id_unique[1:5])

##### ATTACH DIAGNOSIS AND BRAIN REGION #####

#Cerebellum_QCd_entrez_id_unique
#Temporal_cortex_QCd_entrez_id_unique

# "Tissue", "MMSE", "Gender", "Diagnosis_sub_cat", "Diagnosis", "BRAAK", "APOE", "Age"

head(Cerebellum_pheno)

Cerebellum_pheno_to_attach<-Cerebellum_pheno
colnames(Cerebellum_pheno_to_attach)<-c("Diagnosis", "Gender", "Age", "APOE", "Tissue", "MMSE", "Diagnosis_sub_cat",  "BRAAK")

Cerebellum_pheno_to_attach$Tissue<-"Cerebellum"
Cerebellum_pheno_to_attach$MMSE<-"Unknown"
Cerebellum_pheno_to_attach$Diagnosis_sub_cat<-"Unknown"
Cerebellum_pheno_to_attach$BRAAK<-"Unknown"
Cerebellum_pheno_to_attach<-Cerebellum_pheno_to_attach[1:8]

head(Cerebellum_pheno_to_attach)

# attach predicted gender

Cerebellum_pheno_to_attach<-merge(Cerebellum_pheno_to_attach, Cerebellum_gender_comparison[2], by="row.names")
rownames(Cerebellum_pheno_to_attach)<-Cerebellum_pheno_to_attach$Row.names
Cerebellum_pheno_to_attach$Row.names<-NULL

# change gender to clinical gender
colnames(Cerebellum_pheno_to_attach)[2]<-"Clinical_Gender"

# order pheno information

Cerebellum_pheno_to_attach<-Cerebellum_pheno_to_attach[,order(names(Cerebellum_pheno_to_attach), decreasing = T)]

# should be 9 colnames - Diagnosis", "Tissue", "Age", "Clinical_Gender", "Predicted_Gender", "MMSE", "BRAAK", "APOE", "Diagnosis_sub_cat"
colnames(Cerebellum_pheno_to_attach)

# same for temporal cortex

head(Temporal_cortex_pheno)

Temporal_cortex_pheno_to_attach<-Temporal_cortex_pheno
colnames(Temporal_cortex_pheno_to_attach)<-c("Diagnosis", "Gender", "Age", "APOE", "Tissue", "MMSE", "Diagnosis_sub_cat",  "BRAAK")

Temporal_cortex_pheno_to_attach$Tissue<-"Temporal_Cortex"
Temporal_cortex_pheno_to_attach$MMSE<-"Unknown"
Temporal_cortex_pheno_to_attach$Diagnosis_sub_cat<-"Unknown"
Temporal_cortex_pheno_to_attach$BRAAK<-"Unknown"
Temporal_cortex_pheno_to_attach<-Temporal_cortex_pheno_to_attach[1:8]

head(Temporal_cortex_pheno_to_attach)

# attach predicted gender

Temporal_cortex_pheno_to_attach<-merge(Temporal_cortex_pheno_to_attach, Temporal_cortex_gender_comparison[2], by="row.names")
rownames(Temporal_cortex_pheno_to_attach)<-Temporal_cortex_pheno_to_attach$Row.names
Temporal_cortex_pheno_to_attach$Row.names<-NULL

# change gender to clinical gender
colnames(Temporal_cortex_pheno_to_attach)[2]<-"Clinical_Gender"

# order pheno information

Temporal_cortex_pheno_to_attach<-Temporal_cortex_pheno_to_attach[,order(names(Temporal_cortex_pheno_to_attach), decreasing = T)]

# should be 9 colnames - Diagnosis", "Tissue", "Age", "Clinical_Gender", "Predicted_Gender", "MMSE", "BRAAK", "APOE", "Diagnosis_sub_cat"
colnames(Temporal_cortex_pheno_to_attach)

# # merge Temporal_cortex pheno info with expr

Temporal_cortex_QCd_entrez_id_unique_pheno<-merge(Temporal_cortex_pheno_to_attach, Temporal_cortex_QCd_entrez_id_unique, by="row.names")
head(Temporal_cortex_QCd_entrez_id_unique_pheno)[1:5]
rownames(Temporal_cortex_QCd_entrez_id_unique_pheno)<-Temporal_cortex_QCd_entrez_id_unique_pheno$Row.names
Temporal_cortex_QCd_entrez_id_unique_pheno$Row.names<-NULL
head(Temporal_cortex_QCd_entrez_id_unique_pheno)[1:5]
dim(Temporal_cortex_QCd_entrez_id_unique)
dim(Temporal_cortex_pheno_to_attach)
dim(Temporal_cortex_QCd_entrez_id_unique_pheno)

# merge Cerebellum pheno info with expr

Cerebellum_QCd_entrez_id_unique_pheno<-merge(Cerebellum_pheno_to_attach, Cerebellum_QCd_entrez_id_unique, by="row.names")
head(Cerebellum_QCd_entrez_id_unique_pheno)[1:5]
rownames(Cerebellum_QCd_entrez_id_unique_pheno)<-Cerebellum_QCd_entrez_id_unique_pheno$Row.names
Cerebellum_QCd_entrez_id_unique_pheno$Row.names<-NULL
head(Cerebellum_QCd_entrez_id_unique_pheno)[1:5]
dim(Cerebellum_QCd_entrez_id_unique)
dim(Cerebellum_pheno_to_attach)
dim(Cerebellum_QCd_entrez_id_unique_pheno)

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(Temporal_cortex_QCd_entrez_id_unique_pheno, file="AMP_MayoeGWAS_Temporal_Cortex_pre-processed_data2.txt", sep="\t")
write.table(Cerebellum_QCd_entrez_id_unique_pheno, file="AMP_MayoeGWAS_Cerebellum_pre-processed_data2.txt", sep="\t")

##### WRITE LIST OF CONTROL GENES EXPRESSED #####

setwd("/media/hamel/1TB/Projects/Brain_expression/2.Expression_across_brain_regions_in_control_datasets/1.Data")

write(unique(sort(c(Cerebellum_control_expressed_probes_list_F, 
                    Cerebellum_control_expressed_probes_list_M))), file="AMP_MayoeGWAS_Cerebellum.txt")

write(unique(sort(c(Temporal_cortex_control_expressed_probes_list_F, 
                    Temporal_cortex_control_expressed_probes_list_M))), file="AMP_MayoeGWAS_Temporal_Cortex.txt")

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="AMP_MayoeGWAS_data_processing2.Rdata")

