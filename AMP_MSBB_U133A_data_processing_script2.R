
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                   AMP - download/process MSBB - U133A  - PIEPLINE V02                                  #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

##### DESCRIPTION OF ANALYSIS ####
## MSBB - Mount Sinai Brain Bank
## Download raw data and process
## 
## MICROARRAY PLATFORM - Afymetrix
## EXPRESSION CHIP - HG-U133A
## 961 samples - 570 samples belonging to regions used
##
## REPROCESSING DATA USING NEW PIPELINE
##
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir_U133A="/media/hamel/Workspace/Projects/Brain_expression/1.Data/AD/AMP_MSBB/Raw_Data/AffymetrixHG-U133A"

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MSBB/Pre-Processing/AffymetrixHG-U133A"

pheno_dir="/media/hamel/Workspace/Projects/Brain_expression/1.Data/AD/AMP_MSBB/Phenotype_info"

setwd(work_dir)

dir.create(paste(work_dir,"pca_plots", sep="/"))

pca_dir=paste(work_dir,"pca_plots", sep="/")

dir.create(paste(work_dir,"clean_data", sep="/"))

clean_data_dir=paste(work_dir,"clean_data", sep="/")

dir.create(paste(work_dir,"sample_network_plots", sep="/"))

sample_network_dir=paste(work_dir,"sample_network_plots", sep="/")

##### LOAD LIBRARIES ####

library(Biobase)
library(affy)
library(lumi)
library(biomaRt)
library(WGCNA)
library(pamr)
library(limma)
library(sva)
library(hgu133a.db)
library(reshape)
library(massiR)

###### DOWNLOAD DATA #####

# data downloaded using bash script.

###### READ DATA #####

setwd(data_dir_U133A)

MyData<-ReadAffy()

##### PRE-PROCESS #####

#background correct

brain_data_background_corrected<-mas5(MyData)

#normalise

brain_data_normalised<-rsn(log2(exprs(brain_data_background_corrected)))

brain_data_normalised_as_data_frame<-as.data.frame(brain_data_normalised)

head(brain_data_normalised_as_data_frame)[1:5]

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

setwd(pheno_dir)

phenotype_data<-read.table("AMP-AD_MSBB_MSSM_AffymetrixArray_sample_key.txt", head=T)

phenotype_data2<-read.table("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv", head=T)

brain_regions<-read.table("Brain regions explained.csv", head=T)

head(phenotype_data)

head(phenotype_data2)

colnames(phenotype_data2)[1]<-"BB_ID"

head(brain_regions)

table(phenotype_data$REGION)
table(brain_regions$REGION)

# merge phenotype data and brain region

phenotype_data_brain_region<-merge(phenotype_data, brain_regions, by="REGION")

phenotype_data_brain_region<-merge(phenotype_data_brain_region, phenotype_data2, by="BB_ID")

head(phenotype_data_brain_region)

# create column with individual and brain region

phenotype_data_brain_region$individual_ID<-paste(phenotype_data_brain_region$Full_Region_Name, phenotype_data_brain_region$BB_ID, sep="_")
head(phenotype_data_brain_region)
anyDuplicated(phenotype_data_brain_region$individual_ID)

# replace expression matric colnames with individual ID (BB_ID_brain_region)

head(brain_data_normalised_as_data_frame)[1:5]

phenotype_data_brain_region[phenotype_data_brain_region$fileName==colnames(brain_data_normalised_as_data_frame)[1],][,17]

for (x in 1:dim(brain_data_normalised_as_data_frame)[2]){
  colnames(brain_data_normalised_as_data_frame)[x]<-phenotype_data_brain_region[phenotype_data_brain_region$fileName==colnames(brain_data_normalised_as_data_frame)[x],][,17]
}

head(brain_data_normalised_as_data_frame)[1:5]

colnames(brain_data_normalised_as_data_frame)

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(phenotype_data))

##### GENDER CHECK ##### 

# Tried to convert probe ID from massiR to entrez and then back to affy U133a chip - less probes hit.
# Using directly affy_hg_u133_plus_2 probe list
#
#

head(phenotype_data_brain_region)

# gender_information
gender_info<-unique(phenotype_data_brain_region[c(8,17)])
rownames(gender_info)<-gender_info$individual_ID
gender_info$individual_ID<-NULL
colnames(gender_info)<-"Gender"
table(gender_info$Gender)

female_samples<-subset(gender_info, Gender=="F")

male_samples<-subset(gender_info, Gender=="M")

head(female_samples)
head(male_samples)

head(brain_data_normalised_as_data_frame)[1:5]

# get Y choromosome genes
data(y.probes)
names(y.probes)

y_chromo_probes <- data.frame(y.probes["affy_hg_u133_plus_2"])

# extract Y chromosome genes from dataset
eset.select.out <- massi_select(brain_data_normalised_as_data_frame, y_chromo_probes)

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

gender_info$Gender<-as.character(gender_info$Gender)
gender_info[gender_info$Gender=="M",]<-"male"
gender_info[gender_info$Gender=="F",]<-"female"

#merge
gender_comparison<-merge(gender_info, predicted_gender, by="row.names")
rownames(gender_comparison)<-gender_comparison$Row.names
gender_comparison$Row.names<-NULL
colnames(gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
head(gender_comparison)

#differences

Gender_missmatch<-gender_comparison[which(gender_comparison$Clinical_Gender!=gender_comparison$Predicted_Gender),]
Gender_missmatch

# look at missmatch 

Diagnosis[rownames(Gender_missmatch),]

gender_check_results<-eset.results$massi.results

gender_check_results[gender_check_results$ID %in% rownames(Gender_missmatch),]

##### PROBE ID DETECTION #####

# separate case control - Factor.Value..disease. column 

head(phenotype_data_brain_region)

# assign diagnosis status to samples based on CDR - 0=control, 1=MCI, 2-5=AD (based on AMP website)

phenotype_data_brain_region$Diagnosis<-"Unknown"
names(phenotype_data_brain_region)

phenotype_data_brain_region[phenotype_data_brain_region$CDR=="0",18]<-"Control"
phenotype_data_brain_region[phenotype_data_brain_region$CDR=="0.5",18]<-"MCI"
phenotype_data_brain_region[phenotype_data_brain_region$CDR=="1"|
                              phenotype_data_brain_region$CDR=="2"|
                              phenotype_data_brain_region$CDR=="3"|
                              phenotype_data_brain_region$CDR=="4"|
                              phenotype_data_brain_region$CDR=="5",18]<-"AD"

table(phenotype_data_brain_region$Diagnosis)
table(phenotype_data_brain_region$CDR)

case_ID<-phenotype_data_brain_region[phenotype_data_brain_region$Diagnosis=="AD",17]

control_ID<-phenotype_data_brain_region[phenotype_data_brain_region$Diagnosis=="Control",17]

length(case_ID)

length(control_ID)

case_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%case_ID]
control_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%control_ID]

head(case_exprs)
dim(case_exprs)
head(control_exprs)
dim(control_exprs)


# separate by brain region - keep only following - ones that map to inhouse regions

#frontal pole
#precentral gyrus
#inferior frontal gyrus
#dorsolateral prefrontal cortex
#superior parietal lobule
#prefrontal cortex
#parahippocampal gyrus
#hippocampus
#inferior temporal gyrus
#middle temporal gyrus
#superior temporal gyrus
#temporal pole
#amygdala

Frontal_Pole<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Frontal_Pole",17]
Precentral_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Precentral_Gyrus",17]
Inferior_Frontal_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Inferior_Frontal_Gyrus",17]
Dorsolateral_Prefrontal_Cortex<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Dorsolateral_Prefrontal_Cortex",17]
Superior_Parietal_Lobule<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Superior_Parietal_Lobule",17]
Prefrontal_Cortex<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Prefrontal_Cortex",17]
Parahippocampal_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Parahippocampal_Gyrus",17]
Hippocampus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Hippocampus",17]
Inferior_Temporal_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Inferior_Temporal_Gyrus",17]
Middle_Temporal_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Middle_Temporal_Gyrus",17]
Superior_Temporal_Gyrus<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Superior_Temporal_Gyrus",17]
Temporal_Pole<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Temporal_Pole",17]
Amygdala<-phenotype_data_brain_region[phenotype_data_brain_region$Full_Region_Name=="Amygdala",17]

#check each list

Frontal_Pole
Precentral_Gyrus
Inferior_Frontal_Gyrus
Dorsolateral_Prefrontal_Cortex
Superior_Parietal_Lobule
Prefrontal_Cortex
Parahippocampal_Gyrus
Hippocampus
Inferior_Temporal_Gyrus
Middle_Temporal_Gyrus
Superior_Temporal_Gyrus
Temporal_Pole
Amygdala

# separate by brain region and disorder

Frontal_Pole_case_exprs<-case_exprs[,colnames(case_exprs)%in%Frontal_Pole]
Precentral_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Precentral_Gyrus]
Inferior_Frontal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Inferior_Frontal_Gyrus]
Dorsolateral_Prefrontal_Cortex_case_exprs<-case_exprs[,colnames(case_exprs)%in%Dorsolateral_Prefrontal_Cortex]
Superior_Parietal_Lobule_case_exprs<-case_exprs[,colnames(case_exprs)%in%Superior_Parietal_Lobule]
Prefrontal_Cortex_case_exprs<-case_exprs[,colnames(case_exprs)%in%Prefrontal_Cortex]
Parahippocampal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Parahippocampal_Gyrus]
Hippocampus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Hippocampus]
Inferior_Temporal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Inferior_Temporal_Gyrus]
Middle_Temporal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Middle_Temporal_Gyrus]
Superior_Temporal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Superior_Temporal_Gyrus]
Temporal_Pole_case_exprs<-case_exprs[,colnames(case_exprs)%in%Temporal_Pole]
Amygdala_case_exprs<-case_exprs[,colnames(case_exprs)%in%Amygdala]

dim(Frontal_Pole_case_exprs)
dim(Precentral_Gyrus_case_exprs)
dim(Inferior_Frontal_Gyrus_case_exprs)
dim(Dorsolateral_Prefrontal_Cortex_case_exprs)
dim(Superior_Parietal_Lobule_case_exprs)
dim(Prefrontal_Cortex_case_exprs)
dim(Parahippocampal_Gyrus_case_exprs)
dim(Hippocampus_case_exprs)
dim(Inferior_Temporal_Gyrus_case_exprs)
dim(Middle_Temporal_Gyrus_case_exprs)
dim(Superior_Temporal_Gyrus_case_exprs)
dim(Temporal_Pole_case_exprs)
dim(Amygdala_case_exprs)

Frontal_Pole_control_exprs<-control_exprs[,colnames(control_exprs)%in%Frontal_Pole]
Precentral_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Precentral_Gyrus]
Inferior_Frontal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Inferior_Frontal_Gyrus]
Dorsolateral_Prefrontal_Cortex_control_exprs<-control_exprs[,colnames(control_exprs)%in%Dorsolateral_Prefrontal_Cortex]
Superior_Parietal_Lobule_control_exprs<-control_exprs[,colnames(control_exprs)%in%Superior_Parietal_Lobule]
Prefrontal_Cortex_control_exprs<-control_exprs[,colnames(control_exprs)%in%Prefrontal_Cortex]
Parahippocampal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Parahippocampal_Gyrus]
Hippocampus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Hippocampus]
Inferior_Temporal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Inferior_Temporal_Gyrus]
Middle_Temporal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Middle_Temporal_Gyrus]
Superior_Temporal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Superior_Temporal_Gyrus]
Temporal_Pole_control_exprs<-control_exprs[,colnames(control_exprs)%in%Temporal_Pole]
Amygdala_control_exprs<-control_exprs[,colnames(control_exprs)%in%Amygdala]

dim(Frontal_Pole_control_exprs)
dim(Precentral_Gyrus_control_exprs)
dim(Inferior_Frontal_Gyrus_control_exprs)
dim(Dorsolateral_Prefrontal_Cortex_control_exprs)
dim(Superior_Parietal_Lobule_control_exprs)
dim(Prefrontal_Cortex_control_exprs)
dim(Parahippocampal_Gyrus_control_exprs)
dim(Hippocampus_control_exprs)
dim(Inferior_Temporal_Gyrus_control_exprs)
dim(Middle_Temporal_Gyrus_control_exprs)
dim(Superior_Temporal_Gyrus_control_exprs)
dim(Temporal_Pole_control_exprs)
dim(Amygdala_control_exprs)

# dropped Amygdala

# separate by gender
Frontal_Pole_case_exprs_F<-Frontal_Pole_case_exprs[colnames(Frontal_Pole_case_exprs)%in%rownames(female_samples)]
Precentral_Gyrus_case_exprs_F<-Precentral_Gyrus_case_exprs[colnames(Precentral_Gyrus_case_exprs)%in%rownames(female_samples)]
Inferior_Frontal_Gyrus_case_exprs_F<-Inferior_Frontal_Gyrus_case_exprs[colnames(Inferior_Frontal_Gyrus_case_exprs)%in%rownames(female_samples)]
Dorsolateral_Prefrontal_Cortex_case_exprs_F<-Dorsolateral_Prefrontal_Cortex_case_exprs[colnames(Dorsolateral_Prefrontal_Cortex_case_exprs)%in%rownames(female_samples)]
Superior_Parietal_Lobule_case_exprs_F<-Superior_Parietal_Lobule_case_exprs[colnames(Superior_Parietal_Lobule_case_exprs)%in%rownames(female_samples)]
Prefrontal_Cortex_case_exprs_F<-Prefrontal_Cortex_case_exprs[colnames(Prefrontal_Cortex_case_exprs)%in%rownames(female_samples)]
Parahippocampal_Gyrus_case_exprs_F<-Parahippocampal_Gyrus_case_exprs[colnames(Parahippocampal_Gyrus_case_exprs)%in%rownames(female_samples)]
Hippocampus_case_exprs_F<-Hippocampus_case_exprs[colnames(Hippocampus_case_exprs)%in%rownames(female_samples)]
Inferior_Temporal_Gyrus_case_exprs_F<-Inferior_Temporal_Gyrus_case_exprs[colnames(Inferior_Temporal_Gyrus_case_exprs)%in%rownames(female_samples)]
Middle_Temporal_Gyrus_case_exprs_F<-Middle_Temporal_Gyrus_case_exprs[colnames(Middle_Temporal_Gyrus_case_exprs)%in%rownames(female_samples)]
Superior_Temporal_Gyrus_case_exprs_F<-Superior_Temporal_Gyrus_case_exprs[colnames(Superior_Temporal_Gyrus_case_exprs)%in%rownames(female_samples)]
Temporal_Pole_case_exprs_F<-Temporal_Pole_case_exprs[colnames(Temporal_Pole_case_exprs)%in%rownames(female_samples)]
Frontal_Pole_control_exprs_F<-Frontal_Pole_control_exprs[colnames(Frontal_Pole_control_exprs)%in%rownames(female_samples)]
Precentral_Gyrus_control_exprs_F<-Precentral_Gyrus_control_exprs[colnames(Precentral_Gyrus_control_exprs)%in%rownames(female_samples)]
Inferior_Frontal_Gyrus_control_exprs_F<-Inferior_Frontal_Gyrus_control_exprs[colnames(Inferior_Frontal_Gyrus_control_exprs)%in%rownames(female_samples)]
Dorsolateral_Prefrontal_Cortex_control_exprs_F<-Dorsolateral_Prefrontal_Cortex_control_exprs[colnames(Dorsolateral_Prefrontal_Cortex_control_exprs)%in%rownames(female_samples)]
Superior_Parietal_Lobule_control_exprs_F<-Superior_Parietal_Lobule_control_exprs[colnames(Superior_Parietal_Lobule_control_exprs)%in%rownames(female_samples)]
Prefrontal_Cortex_control_exprs_F<-Prefrontal_Cortex_control_exprs[colnames(Prefrontal_Cortex_control_exprs)%in%rownames(female_samples)]
Parahippocampal_Gyrus_control_exprs_F<-Parahippocampal_Gyrus_control_exprs[colnames(Parahippocampal_Gyrus_control_exprs)%in%rownames(female_samples)]
Hippocampus_control_exprs_F<-Hippocampus_control_exprs[colnames(Hippocampus_control_exprs)%in%rownames(female_samples)]
Inferior_Temporal_Gyrus_control_exprs_F<-Inferior_Temporal_Gyrus_control_exprs[colnames(Inferior_Temporal_Gyrus_control_exprs)%in%rownames(female_samples)]
Middle_Temporal_Gyrus_control_exprs_F<-Middle_Temporal_Gyrus_control_exprs[colnames(Middle_Temporal_Gyrus_control_exprs)%in%rownames(female_samples)]
Superior_Temporal_Gyrus_control_exprs_F<-Superior_Temporal_Gyrus_control_exprs[colnames(Superior_Temporal_Gyrus_control_exprs)%in%rownames(female_samples)]
Temporal_Pole_control_exprs_F<-Temporal_Pole_control_exprs[colnames(Temporal_Pole_control_exprs)%in%rownames(female_samples)]

Frontal_Pole_case_exprs_M<-Frontal_Pole_case_exprs[colnames(Frontal_Pole_case_exprs)%in%rownames(male_samples)]
Precentral_Gyrus_case_exprs_M<-Precentral_Gyrus_case_exprs[colnames(Precentral_Gyrus_case_exprs)%in%rownames(male_samples)]
Inferior_Frontal_Gyrus_case_exprs_M<-Inferior_Frontal_Gyrus_case_exprs[colnames(Inferior_Frontal_Gyrus_case_exprs)%in%rownames(male_samples)]
Dorsolateral_Prefrontal_Cortex_case_exprs_M<-Dorsolateral_Prefrontal_Cortex_case_exprs[colnames(Dorsolateral_Prefrontal_Cortex_case_exprs)%in%rownames(male_samples)]
Superior_Parietal_Lobule_case_exprs_M<-Superior_Parietal_Lobule_case_exprs[colnames(Superior_Parietal_Lobule_case_exprs)%in%rownames(male_samples)]
Prefrontal_Cortex_case_exprs_M<-Prefrontal_Cortex_case_exprs[colnames(Prefrontal_Cortex_case_exprs)%in%rownames(male_samples)]
Parahippocampal_Gyrus_case_exprs_M<-Parahippocampal_Gyrus_case_exprs[colnames(Parahippocampal_Gyrus_case_exprs)%in%rownames(male_samples)]
Hippocampus_case_exprs_M<-Hippocampus_case_exprs[colnames(Hippocampus_case_exprs)%in%rownames(male_samples)]
Inferior_Temporal_Gyrus_case_exprs_M<-Inferior_Temporal_Gyrus_case_exprs[colnames(Inferior_Temporal_Gyrus_case_exprs)%in%rownames(male_samples)]
Middle_Temporal_Gyrus_case_exprs_M<-Middle_Temporal_Gyrus_case_exprs[colnames(Middle_Temporal_Gyrus_case_exprs)%in%rownames(male_samples)]
Superior_Temporal_Gyrus_case_exprs_M<-Superior_Temporal_Gyrus_case_exprs[colnames(Superior_Temporal_Gyrus_case_exprs)%in%rownames(male_samples)]
Temporal_Pole_case_exprs_M<-Temporal_Pole_case_exprs[colnames(Temporal_Pole_case_exprs)%in%rownames(male_samples)]
Frontal_Pole_control_exprs_M<-Frontal_Pole_control_exprs[colnames(Frontal_Pole_control_exprs)%in%rownames(male_samples)]
Precentral_Gyrus_control_exprs_M<-Precentral_Gyrus_control_exprs[colnames(Precentral_Gyrus_control_exprs)%in%rownames(male_samples)]
Inferior_Frontal_Gyrus_control_exprs_M<-Inferior_Frontal_Gyrus_control_exprs[colnames(Inferior_Frontal_Gyrus_control_exprs)%in%rownames(male_samples)]
Dorsolateral_Prefrontal_Cortex_control_exprs_M<-Dorsolateral_Prefrontal_Cortex_control_exprs[colnames(Dorsolateral_Prefrontal_Cortex_control_exprs)%in%rownames(male_samples)]
Superior_Parietal_Lobule_control_exprs_M<-Superior_Parietal_Lobule_control_exprs[colnames(Superior_Parietal_Lobule_control_exprs)%in%rownames(male_samples)]
Prefrontal_Cortex_control_exprs_M<-Prefrontal_Cortex_control_exprs[colnames(Prefrontal_Cortex_control_exprs)%in%rownames(male_samples)]
Parahippocampal_Gyrus_control_exprs_M<-Parahippocampal_Gyrus_control_exprs[colnames(Parahippocampal_Gyrus_control_exprs)%in%rownames(male_samples)]
Hippocampus_control_exprs_M<-Hippocampus_control_exprs[colnames(Hippocampus_control_exprs)%in%rownames(male_samples)]
Inferior_Temporal_Gyrus_control_exprs_M<-Inferior_Temporal_Gyrus_control_exprs[colnames(Inferior_Temporal_Gyrus_control_exprs)%in%rownames(male_samples)]
Middle_Temporal_Gyrus_control_exprs_M<-Middle_Temporal_Gyrus_control_exprs[colnames(Middle_Temporal_Gyrus_control_exprs)%in%rownames(male_samples)]
Superior_Temporal_Gyrus_control_exprs_M<-Superior_Temporal_Gyrus_control_exprs[colnames(Superior_Temporal_Gyrus_control_exprs)%in%rownames(male_samples)]
Temporal_Pole_control_exprs_M<-Temporal_Pole_control_exprs[colnames(Temporal_Pole_control_exprs)%in%rownames(male_samples)]

# check all data frames
head(Frontal_Pole_case_exprs_F)[1:5]
head(Precentral_Gyrus_case_exprs_F)[1:5]
head(Inferior_Frontal_Gyrus_case_exprs_F)[1:5]
head(Dorsolateral_Prefrontal_Cortex_case_exprs_F)[1:5]
head(Superior_Parietal_Lobule_case_exprs_F)[1:5]
head(Prefrontal_Cortex_case_exprs_F)[1:5]
head(Parahippocampal_Gyrus_case_exprs_F)[1:5]
head(Hippocampus_case_exprs_F)[1:5]
head(Inferior_Temporal_Gyrus_case_exprs_F)[1:5]
head(Middle_Temporal_Gyrus_case_exprs_F)[1:5]
head(Superior_Temporal_Gyrus_case_exprs_F)[1:5]
head(Temporal_Pole_case_exprs_F)[1:5]
head(Frontal_Pole_control_exprs_F)[1:5]
head(Precentral_Gyrus_control_exprs_F)[1:5]
head(Inferior_Frontal_Gyrus_control_exprs_F)[1:5]
head(Dorsolateral_Prefrontal_Cortex_control_exprs_F)[1:5]
head(Superior_Parietal_Lobule_control_exprs_F)[1:5]
head(Prefrontal_Cortex_control_exprs_F)[1:5]
head(Parahippocampal_Gyrus_control_exprs_F)[1:5]
head(Hippocampus_control_exprs_F)[1:5]
head(Inferior_Temporal_Gyrus_control_exprs_F)[1:5]
head(Middle_Temporal_Gyrus_control_exprs_F)[1:5]
head(Superior_Temporal_Gyrus_control_exprs_F)[1:5]
head(Temporal_Pole_control_exprs_F)[1:5]


head(Frontal_Pole_case_exprs_M)[1:5]
head(Precentral_Gyrus_case_exprs_M)[1:5]
head(Inferior_Frontal_Gyrus_case_exprs_M)[1:5]
head(Dorsolateral_Prefrontal_Cortex_case_exprs_M)[1:5]
head(Superior_Parietal_Lobule_case_exprs_M)[1:5]
head(Prefrontal_Cortex_case_exprs_M)[1:5]
head(Parahippocampal_Gyrus_case_exprs_M)[1:5]
head(Hippocampus_case_exprs_M)[1:5]
head(Inferior_Temporal_Gyrus_case_exprs_M)[1:5]
head(Middle_Temporal_Gyrus_case_exprs_M)[1:5]
head(Superior_Temporal_Gyrus_case_exprs_M)[1:5]
head(Temporal_Pole_case_exprs_M)[1:5]
head(Frontal_Pole_control_exprs_M)[1:5]
head(Precentral_Gyrus_control_exprs_M)
head(Inferior_Frontal_Gyrus_control_exprs_M)
head(Dorsolateral_Prefrontal_Cortex_control_exprs_M)[1:5]
head(Superior_Parietal_Lobule_control_exprs_M)
head(Prefrontal_Cortex_control_exprs_M)
head(Parahippocampal_Gyrus_control_exprs_M)
head(Hippocampus_control_exprs_M)
head(Inferior_Temporal_Gyrus_control_exprs_M)
head(Middle_Temporal_Gyrus_control_exprs_M)
head(Superior_Temporal_Gyrus_control_exprs_M)
head(Temporal_Pole_control_exprs_M)

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

# apply function 

Frontal_Pole_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Frontal_Pole_case_exprs_F, 0.9, 0.8)
length(Frontal_Pole_case_exprs_F_expressed_probes_list)

Precentral_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Precentral_Gyrus_case_exprs_F, 0.9, 0.8)
length(Precentral_Gyrus_case_exprs_F_expressed_probes_list)

Inferior_Frontal_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Inferior_Frontal_Gyrus_case_exprs_F, 0.9, 0.8)
length(Inferior_Frontal_Gyrus_case_exprs_F_expressed_probes_list)

Dorsolateral_Prefrontal_Cortex_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Dorsolateral_Prefrontal_Cortex_case_exprs_F, 0.9, 0.8)
length(Dorsolateral_Prefrontal_Cortex_case_exprs_F_expressed_probes_list)

Superior_Parietal_Lobule_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Superior_Parietal_Lobule_case_exprs_F, 0.9, 0.8)
length(Superior_Parietal_Lobule_case_exprs_F_expressed_probes_list)

Prefrontal_Cortex_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Prefrontal_Cortex_case_exprs_F, 0.9, 0.8)
length(Prefrontal_Cortex_case_exprs_F_expressed_probes_list)

Parahippocampal_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Parahippocampal_Gyrus_case_exprs_F, 0.9, 0.8)
length(Parahippocampal_Gyrus_case_exprs_F_expressed_probes_list)

Hippocampus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Hippocampus_case_exprs_F, 0.9, 0.8)
length(Hippocampus_case_exprs_F_expressed_probes_list)

Inferior_Temporal_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Inferior_Temporal_Gyrus_case_exprs_F, 0.9, 0.8)
length(Inferior_Temporal_Gyrus_case_exprs_F_expressed_probes_list)

Middle_Temporal_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Middle_Temporal_Gyrus_case_exprs_F, 0.9, 0.8)
length(Middle_Temporal_Gyrus_case_exprs_F_expressed_probes_list)

Superior_Temporal_Gyrus_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Superior_Temporal_Gyrus_case_exprs_F, 0.9, 0.8)
length(Superior_Temporal_Gyrus_case_exprs_F_expressed_probes_list)

Temporal_Pole_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Temporal_Pole_case_exprs_F, 0.9, 0.8)
length(Temporal_Pole_case_exprs_F_expressed_probes_list)

Frontal_Pole_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Frontal_Pole_control_exprs_F, 0.9, 0.8)
length(Frontal_Pole_control_exprs_F_expressed_probes_list)

Precentral_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Precentral_Gyrus_control_exprs_F, 0.9, 0.8)
length(Precentral_Gyrus_control_exprs_F_expressed_probes_list)

Inferior_Frontal_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Inferior_Frontal_Gyrus_control_exprs_F, 0.9, 0.8)
length(Inferior_Frontal_Gyrus_control_exprs_F_expressed_probes_list)

Dorsolateral_Prefrontal_Cortex_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Dorsolateral_Prefrontal_Cortex_control_exprs_F, 0.9, 0.8)
length(Dorsolateral_Prefrontal_Cortex_control_exprs_F_expressed_probes_list)

Superior_Parietal_Lobule_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Superior_Parietal_Lobule_control_exprs_F, 0.9, 0.8)
length(Superior_Parietal_Lobule_control_exprs_F_expressed_probes_list)

Prefrontal_Cortex_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Prefrontal_Cortex_control_exprs_F, 0.9, 0.8)
length(Prefrontal_Cortex_control_exprs_F_expressed_probes_list)

Parahippocampal_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Parahippocampal_Gyrus_control_exprs_F, 0.9, 0.8)
length(Parahippocampal_Gyrus_control_exprs_F_expressed_probes_list)

Hippocampus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Hippocampus_control_exprs_F, 0.9, 0.8)
length(Hippocampus_control_exprs_F_expressed_probes_list)

Inferior_Temporal_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Inferior_Temporal_Gyrus_control_exprs_F, 0.9, 0.8)
length(Inferior_Temporal_Gyrus_control_exprs_F_expressed_probes_list)

Middle_Temporal_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Middle_Temporal_Gyrus_control_exprs_F, 0.9, 0.8)
length(Middle_Temporal_Gyrus_control_exprs_F_expressed_probes_list)

Superior_Temporal_Gyrus_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Superior_Temporal_Gyrus_control_exprs_F, 0.9, 0.8)
length(Superior_Temporal_Gyrus_control_exprs_F_expressed_probes_list)

Temporal_Pole_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Temporal_Pole_control_exprs_F, 0.9, 0.8)
length(Temporal_Pole_control_exprs_F_expressed_probes_list)

Frontal_Pole_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Frontal_Pole_case_exprs_M, 0.9, 0.8)
length(Frontal_Pole_case_exprs_M_expressed_probes_list)

Precentral_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Precentral_Gyrus_case_exprs_M, 0.9, 0.8)
length(Precentral_Gyrus_case_exprs_M_expressed_probes_list)

Inferior_Frontal_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Inferior_Frontal_Gyrus_case_exprs_M, 0.9, 0.8)
length(Inferior_Frontal_Gyrus_case_exprs_M_expressed_probes_list)

Dorsolateral_Prefrontal_Cortex_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Dorsolateral_Prefrontal_Cortex_case_exprs_M, 0.9, 0.8)
length(Dorsolateral_Prefrontal_Cortex_case_exprs_M_expressed_probes_list)

Superior_Parietal_Lobule_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Superior_Parietal_Lobule_case_exprs_M, 0.9, 0.8)
length(Superior_Parietal_Lobule_case_exprs_M_expressed_probes_list)

Prefrontal_Cortex_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Prefrontal_Cortex_case_exprs_M, 0.9, 0.8)
length(Prefrontal_Cortex_case_exprs_M_expressed_probes_list)

Parahippocampal_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Parahippocampal_Gyrus_case_exprs_M, 0.9, 0.8)
length(Parahippocampal_Gyrus_case_exprs_M_expressed_probes_list)

Hippocampus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Hippocampus_case_exprs_M, 0.9, 0.8)
length(Hippocampus_case_exprs_M_expressed_probes_list)

Inferior_Temporal_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Inferior_Temporal_Gyrus_case_exprs_M, 0.9, 0.8)
length(Inferior_Temporal_Gyrus_case_exprs_M_expressed_probes_list)

Middle_Temporal_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Middle_Temporal_Gyrus_case_exprs_M, 0.9, 0.8)
length(Middle_Temporal_Gyrus_case_exprs_M_expressed_probes_list)

Superior_Temporal_Gyrus_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Superior_Temporal_Gyrus_case_exprs_M, 0.9, 0.8)
length(Superior_Temporal_Gyrus_case_exprs_M_expressed_probes_list)

Temporal_Pole_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Temporal_Pole_case_exprs_M, 0.9, 0.8)
length(Temporal_Pole_case_exprs_M_expressed_probes_list)

Frontal_Pole_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Frontal_Pole_control_exprs_M, 0.9, 0.8)
length(Frontal_Pole_control_exprs_M_expressed_probes_list)

Precentral_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Precentral_Gyrus_control_exprs_M, 0.9, 0.8)
length(Precentral_Gyrus_control_exprs_M_expressed_probes_list)

Inferior_Frontal_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Inferior_Frontal_Gyrus_control_exprs_M, 0.9, 0.8)
length(Inferior_Frontal_Gyrus_control_exprs_M_expressed_probes_list)

Dorsolateral_Prefrontal_Cortex_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Dorsolateral_Prefrontal_Cortex_control_exprs_M, 0.9, 0.8)
length(Dorsolateral_Prefrontal_Cortex_control_exprs_M_expressed_probes_list)

Superior_Parietal_Lobule_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Superior_Parietal_Lobule_control_exprs_M, 0.9, 0.8)
length(Superior_Parietal_Lobule_control_exprs_M_expressed_probes_list)

Prefrontal_Cortex_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Prefrontal_Cortex_control_exprs_M, 0.9, 0.8)
length(Prefrontal_Cortex_control_exprs_M_expressed_probes_list)

Parahippocampal_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Parahippocampal_Gyrus_control_exprs_M, 0.9, 0.8)
length(Parahippocampal_Gyrus_control_exprs_M_expressed_probes_list)

Hippocampus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Hippocampus_control_exprs_M, 0.9, 0.8)
length(Hippocampus_control_exprs_M_expressed_probes_list)

Inferior_Temporal_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Inferior_Temporal_Gyrus_control_exprs_M, 0.9, 0.8)
length(Inferior_Temporal_Gyrus_control_exprs_M_expressed_probes_list)

Middle_Temporal_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Middle_Temporal_Gyrus_control_exprs_M, 0.9, 0.8)
length(Middle_Temporal_Gyrus_control_exprs_M_expressed_probes_list)

Superior_Temporal_Gyrus_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Superior_Temporal_Gyrus_control_exprs_M, 0.9, 0.8)
length(Superior_Temporal_Gyrus_control_exprs_M_expressed_probes_list)

Temporal_Pole_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Temporal_Pole_control_exprs_M, 0.9, 0.8)
length(Temporal_Pole_control_exprs_M_expressed_probes_list)

# merge list of good probes from both case + control + and all brain regions, sort and keep unique values

good_probe_list<-unique(sort(c(Frontal_Pole_case_exprs_F_expressed_probes_list,
                               Precentral_Gyrus_case_exprs_F_expressed_probes_list,
                               Inferior_Frontal_Gyrus_case_exprs_F_expressed_probes_list,
                               Dorsolateral_Prefrontal_Cortex_case_exprs_F_expressed_probes_list,
                               Superior_Parietal_Lobule_case_exprs_F_expressed_probes_list,
                               Prefrontal_Cortex_case_exprs_F_expressed_probes_list,
                               Parahippocampal_Gyrus_case_exprs_F_expressed_probes_list,
                               Hippocampus_case_exprs_F_expressed_probes_list,
                               Inferior_Temporal_Gyrus_case_exprs_F_expressed_probes_list,
                               Middle_Temporal_Gyrus_case_exprs_F_expressed_probes_list,
                               Superior_Temporal_Gyrus_case_exprs_F_expressed_probes_list,
                               Temporal_Pole_case_exprs_F_expressed_probes_list,
                               Frontal_Pole_control_exprs_F_expressed_probes_list,
                               Precentral_Gyrus_control_exprs_F_expressed_probes_list,
                               Inferior_Frontal_Gyrus_control_exprs_F_expressed_probes_list,
                               Dorsolateral_Prefrontal_Cortex_control_exprs_F_expressed_probes_list,
                               Superior_Parietal_Lobule_control_exprs_F_expressed_probes_list,
                               Prefrontal_Cortex_control_exprs_F_expressed_probes_list,
                               Parahippocampal_Gyrus_control_exprs_F_expressed_probes_list,
                               Hippocampus_control_exprs_F_expressed_probes_list,
                               Inferior_Temporal_Gyrus_control_exprs_F_expressed_probes_list,
                               Middle_Temporal_Gyrus_control_exprs_F_expressed_probes_list,
                               Superior_Temporal_Gyrus_control_exprs_F_expressed_probes_list,
                               Temporal_Pole_control_exprs_F_expressed_probes_list,
                               Frontal_Pole_case_exprs_M_expressed_probes_list,
                               Precentral_Gyrus_case_exprs_M_expressed_probes_list,
                               Inferior_Frontal_Gyrus_case_exprs_M_expressed_probes_list,
                               Dorsolateral_Prefrontal_Cortex_case_exprs_M_expressed_probes_list,
                               Superior_Parietal_Lobule_case_exprs_M_expressed_probes_list,
                               Prefrontal_Cortex_case_exprs_M_expressed_probes_list,
                               Parahippocampal_Gyrus_case_exprs_M_expressed_probes_list,
                               Hippocampus_case_exprs_M_expressed_probes_list,
                               Inferior_Temporal_Gyrus_case_exprs_M_expressed_probes_list,
                               Middle_Temporal_Gyrus_case_exprs_M_expressed_probes_list,
                               Superior_Temporal_Gyrus_case_exprs_M_expressed_probes_list,
                               Temporal_Pole_case_exprs_M_expressed_probes_list,
                               Frontal_Pole_control_exprs_M_expressed_probes_list,
                               Precentral_Gyrus_control_exprs_M_expressed_probes_list,
                               Inferior_Frontal_Gyrus_control_exprs_M_expressed_probes_list,
                               Dorsolateral_Prefrontal_Cortex_control_exprs_M_expressed_probes_list,
                               Superior_Parietal_Lobule_control_exprs_M_expressed_probes_list,
                               Prefrontal_Cortex_control_exprs_M_expressed_probes_list,
                               Parahippocampal_Gyrus_control_exprs_M_expressed_probes_list,
                               Hippocampus_control_exprs_M_expressed_probes_list,
                               Inferior_Temporal_Gyrus_control_exprs_M_expressed_probes_list,
                               Middle_Temporal_Gyrus_control_exprs_M_expressed_probes_list,
                               Superior_Temporal_Gyrus_control_exprs_M_expressed_probes_list,
                               Temporal_Pole_control_exprs_M_expressed_probes_list)))

length(good_probe_list)
dim(brain_data_normalised_as_data_frame)

# extract good probes from dataset

brain_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]
dim(brain_exprs_good_probes)

Frontal_Pole_case_exprs_good_probes<-Frontal_Pole_case_exprs[rownames(Frontal_Pole_case_exprs)%in%good_probe_list,]
Precentral_Gyrus_case_exprs_good_probes<-Precentral_Gyrus_case_exprs[rownames(Precentral_Gyrus_case_exprs)%in%good_probe_list,]
Inferior_Frontal_Gyrus_case_exprs_good_probes<-Inferior_Frontal_Gyrus_case_exprs[rownames(Inferior_Frontal_Gyrus_case_exprs)%in%good_probe_list,]
Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes<-Dorsolateral_Prefrontal_Cortex_case_exprs[rownames(Dorsolateral_Prefrontal_Cortex_case_exprs)%in%good_probe_list,]
Superior_Parietal_Lobule_case_exprs_good_probes<-Superior_Parietal_Lobule_case_exprs[rownames(Superior_Parietal_Lobule_case_exprs)%in%good_probe_list,]
Prefrontal_Cortex_case_exprs_good_probes<-Prefrontal_Cortex_case_exprs[rownames(Prefrontal_Cortex_case_exprs)%in%good_probe_list,]
Parahippocampal_Gyrus_case_exprs_good_probes<-Parahippocampal_Gyrus_case_exprs[rownames(Parahippocampal_Gyrus_case_exprs)%in%good_probe_list,]
Hippocampus_case_exprs_good_probes<-Hippocampus_case_exprs[rownames(Hippocampus_case_exprs)%in%good_probe_list,]
Inferior_Temporal_Gyrus_case_exprs_good_probes<-Inferior_Temporal_Gyrus_case_exprs[rownames(Inferior_Temporal_Gyrus_case_exprs)%in%good_probe_list,]
Middle_Temporal_Gyrus_case_exprs_good_probes<-Middle_Temporal_Gyrus_case_exprs[rownames(Middle_Temporal_Gyrus_case_exprs)%in%good_probe_list,]
Superior_Temporal_Gyrus_case_exprs_good_probes<-Superior_Temporal_Gyrus_case_exprs[rownames(Superior_Temporal_Gyrus_case_exprs)%in%good_probe_list,]
Temporal_Pole_case_exprs_good_probes<-Temporal_Pole_case_exprs[rownames(Temporal_Pole_case_exprs)%in%good_probe_list,]

Frontal_Pole_control_exprs_good_probes<-Frontal_Pole_control_exprs[rownames(Frontal_Pole_control_exprs)%in%good_probe_list,]
Precentral_Gyrus_control_exprs_good_probes<-Precentral_Gyrus_control_exprs[rownames(Precentral_Gyrus_control_exprs)%in%good_probe_list,]
Inferior_Frontal_Gyrus_control_exprs_good_probes<-Inferior_Frontal_Gyrus_control_exprs[rownames(Inferior_Frontal_Gyrus_control_exprs)%in%good_probe_list,]
Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes<-Dorsolateral_Prefrontal_Cortex_control_exprs[rownames(Dorsolateral_Prefrontal_Cortex_control_exprs)%in%good_probe_list,]
Superior_Parietal_Lobule_control_exprs_good_probes<-Superior_Parietal_Lobule_control_exprs[rownames(Superior_Parietal_Lobule_control_exprs)%in%good_probe_list,]
Prefrontal_Cortex_control_exprs_good_probes<-Prefrontal_Cortex_control_exprs[rownames(Prefrontal_Cortex_control_exprs)%in%good_probe_list,]
Parahippocampal_Gyrus_control_exprs_good_probes<-Parahippocampal_Gyrus_control_exprs[rownames(Parahippocampal_Gyrus_control_exprs)%in%good_probe_list,]
Hippocampus_control_exprs_good_probes<-Hippocampus_control_exprs[rownames(Hippocampus_control_exprs)%in%good_probe_list,]
Inferior_Temporal_Gyrus_control_exprs_good_probes<-Inferior_Temporal_Gyrus_control_exprs[rownames(Inferior_Temporal_Gyrus_control_exprs)%in%good_probe_list,]
Middle_Temporal_Gyrus_control_exprs_good_probes<-Middle_Temporal_Gyrus_control_exprs[rownames(Middle_Temporal_Gyrus_control_exprs)%in%good_probe_list,]
Superior_Temporal_Gyrus_control_exprs_good_probes<-Superior_Temporal_Gyrus_control_exprs[rownames(Superior_Temporal_Gyrus_control_exprs)%in%good_probe_list,]
Temporal_Pole_control_exprs_good_probes<-Temporal_Pole_control_exprs[rownames(Temporal_Pole_control_exprs)%in%good_probe_list,]


# check dataframe - all should be 2613

head(Frontal_Pole_case_exprs_good_probes)[1:5]
dim(Frontal_Pole_case_exprs_good_probes)

head(Precentral_Gyrus_case_exprs_good_probes)[1:5]
dim(Precentral_Gyrus_case_exprs_good_probes)

head(Inferior_Frontal_Gyrus_case_exprs_good_probes)[1:5]
dim(Inferior_Frontal_Gyrus_case_exprs_good_probes)

head(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes)[1:5]
dim(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes)

head(Superior_Parietal_Lobule_case_exprs_good_probes)[1:5]
dim(Superior_Parietal_Lobule_case_exprs_good_probes)

head(Prefrontal_Cortex_case_exprs_good_probes)[1:5]
dim(Prefrontal_Cortex_case_exprs_good_probes)

head(Parahippocampal_Gyrus_case_exprs_good_probes)[1:5]
dim(Parahippocampal_Gyrus_case_exprs_good_probes)

head(Hippocampus_case_exprs_good_probes)[1:5]
dim(Hippocampus_case_exprs_good_probes)

head(Inferior_Temporal_Gyrus_case_exprs_good_probes)[1:5]
dim(Inferior_Temporal_Gyrus_case_exprs_good_probes)

head(Middle_Temporal_Gyrus_case_exprs_good_probes)[1:5]
dim(Middle_Temporal_Gyrus_case_exprs_good_probes)

head(Superior_Temporal_Gyrus_case_exprs_good_probes)[1:5]
dim(Superior_Temporal_Gyrus_case_exprs_good_probes)

head(Temporal_Pole_case_exprs_good_probes)[1:5]
dim(Temporal_Pole_case_exprs_good_probes)

head(Frontal_Pole_control_exprs_good_probes)[1:5]
dim(Frontal_Pole_control_exprs_good_probes)

head(Precentral_Gyrus_control_exprs_good_probes)[1:5]
dim(Precentral_Gyrus_control_exprs_good_probes)

head(Inferior_Frontal_Gyrus_control_exprs_good_probes)[1:5]
dim(Inferior_Frontal_Gyrus_control_exprs_good_probes)

head(Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes)[1:5]
dim(Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes)

head(Superior_Parietal_Lobule_control_exprs_good_probes)[1:5]
dim(Superior_Parietal_Lobule_control_exprs_good_probes)

head(Prefrontal_Cortex_control_exprs_good_probes)[1:5]
dim(Prefrontal_Cortex_control_exprs_good_probes)

head(Parahippocampal_Gyrus_control_exprs_good_probes)[1:5]
dim(Parahippocampal_Gyrus_control_exprs_good_probes)

head(Hippocampus_control_exprs_good_probes)[1:5]
dim(Hippocampus_control_exprs_good_probes)

head(Inferior_Temporal_Gyrus_control_exprs_good_probes)[1:5]
dim(Inferior_Temporal_Gyrus_control_exprs_good_probes)

head(Middle_Temporal_Gyrus_control_exprs_good_probes)[1:5]
dim(Middle_Temporal_Gyrus_control_exprs_good_probes)

head(Superior_Temporal_Gyrus_control_exprs_good_probes)[1:5]
dim(Superior_Temporal_Gyrus_control_exprs_good_probes)

head(Temporal_Pole_control_exprs_good_probes)[1:5]
dim(Temporal_Pole_control_exprs_good_probes)

dim(brain_data_normalised)

# create dataframe with good probes - contain all brain regions, case and control

brain_data_normalised_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]
dim(brain_data_normalised_exprs_good_probes)
head(brain_data_normalised_exprs_good_probes)[1:5]

##### GENDER SPECIFIC PROBE PLOTS #####

# get gene symbol list for chip
Gene_symbols_probes <- mappedkeys(hgu133aSYMBOL)

# Convert to a list
Gene_symbols <- as.data.frame(hgu133aSYMBOL[Gene_symbols_probes])

head(Gene_symbols)
dim(Gene_symbols)

#xist gene - 
XIST_probe_ID<-subset(Gene_symbols, symbol=="XIST")
XIST_probe_ID

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
                 NLGN4Y_probe_ID,
                 TMSB4Y_probe_ID,
                 USP9Y_probe_ID,
                 UTY_probe_ID)

#create function to plot
plot_gender_specific_genes<-function(Expression_table, gender_info, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[genes_to_extract[,1],]))
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
  Expression_table_t<-as.data.frame(t(Expression_table))
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

pdf("Gender_specific_gene_plot_and_detectable_cut_off_threshold_used.pdf")
plot_gender_specific_genes(Frontal_Pole_case_exprs, gender_comparison[1], gene_list, 0.9, "Frontal_Pole_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Precentral_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Precentral_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Inferior_Frontal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Inferior_Frontal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Dorsolateral_Prefrontal_Cortex_case_exprs, gender_comparison[1], gene_list, 0.9, "Dorsolateral_Prefrontal_Cortex_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Superior_Parietal_Lobule_case_exprs, gender_comparison[1], gene_list, 0.9, "Superior_Parietal_Lobule_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Prefrontal_Cortex_case_exprs, gender_comparison[1], gene_list, 0.9, "Prefrontal_Cortex_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Parahippocampal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Parahippocampal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Hippocampus_case_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Inferior_Temporal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Inferior_Temporal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Middle_Temporal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Middle_Temporal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Superior_Temporal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Superior_Temporal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Temporal_Pole_case_exprs, gender_comparison[1], gene_list, 0.9, "Temporal_Pole_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Frontal_Pole_control_exprs, gender_comparison[1], gene_list, 0.9, "Frontal_Pole_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Precentral_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Precentral_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Inferior_Frontal_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Inferior_Frontal_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Dorsolateral_Prefrontal_Cortex_control_exprs, gender_comparison[1], gene_list, 0.9, "Dorsolateral_Prefrontal_Cortex_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Superior_Parietal_Lobule_control_exprs, gender_comparison[1], gene_list, 0.9, "Superior_Parietal_Lobule_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Prefrontal_Cortex_control_exprs, gender_comparison[1], gene_list, 0.9, "Prefrontal_Cortex_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Parahippocampal_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Parahippocampal_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Hippocampus_control_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Inferior_Temporal_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Inferior_Temporal_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Middle_Temporal_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Middle_Temporal_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Superior_Temporal_Gyrus_control_exprs, gender_comparison[1], gene_list, 0.9, "Superior_Temporal_Gyrus_control_exprs_gender_specific_genes")
plot_gender_specific_genes(Temporal_Pole_control_exprs, gender_comparison[1], gene_list, 0.9, "Temporal_Pole_control_exprs_gender_specific_genes")
dev.off()

# plot gender missmatches

pdf("Gender_missmatch.pdf")
plot_gender_specific_genes(Parahippocampal_Gyrus_case_exprs, gender_comparison[1], gene_list, 0.9, "Parahippocampal_Gyrus_case_exprs_gender_specific_genes")
plot_gender_specific_genes(Parahippocampal_Gyrus_case_exprs, gender_comparison[2], gene_list, 0.9, "Parahippocampal_Gyrus_case_exprs_gender_specific_genes")
dev.off()

setwd(work_dir)

##### PCA ####

# calculate pca

pca<-prcomp(brain_data_normalised_exprs_good_probes)

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

head(phenotype_data_brain_region)

# phenotype and sample iD
Diagnosis_lookup<-phenotype_data_brain_region[c(5,17,18)]
dim(Diagnosis_lookup)
head(Diagnosis_lookup)
Diagnosis_lookup<-unique(Diagnosis_lookup)
dim(Diagnosis_lookup)
head(Diagnosis_lookup)

#subset to those in expression dataframe
Diagnosis_lookup_temp<-subset(Diagnosis_lookup, individual_ID %in% colnames(brain_data_normalised_exprs_good_probes))
table(Diagnosis_lookup_temp$Diagnosis)
anyDuplicated(Diagnosis_lookup_temp$individual_ID)

rownames(Diagnosis_lookup_temp)<-Diagnosis_lookup_temp$individual_ID

Diagnosis_lookup_temp$individual_ID<-NULL
head(Diagnosis_lookup_temp)

# order of samples in expression data
Diagnosis_temp<-colnames(brain_data_normalised_exprs_good_probes)
head(Diagnosis_temp)

# match order
Diagnosis_pca<-Diagnosis_lookup[match(Diagnosis_temp, rownames(Diagnosis_lookup_temp)),]

# assign color to group
Diagnosis_pca_color<-labels2colors(as.character(Diagnosis_pca))

# pca plot - color by disease - case/control
plot(pca$rotation[,1:2], main=" PCA plot coloured by Disease Status",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomleft', unique(Diagnosis_pca$Diagnosis), fill=unique(Diagnosis_pca_color))

#pca plot - color by brain region

# phenotype and sample iD
brain_region_lookup<-Diagnosis_lookup_temp
brain_region_lookup$Diagnosis<-NULL
colnames(brain_region_lookup)<-"brain_region"
brain_region_lookup

# order of samples in expression data
brain_region_temp<-colnames(brain_data_normalised_exprs_good_probes)
brain_region_temp

# match order
brain_region_pca<-brain_region_lookup[match(brain_region_temp, rownames(brain_region_lookup)),]

# assign color to group
brain_region_pca_color<-labels2colors(as.character(brain_region_pca))

# pca plot - color by disease - case/control
plot(pca$rotation[,1:2], main=" PCA plot coloured by Brain Region",col="black", pch=21,bg=brain_region_pca_color)
legend('bottomright', as.character(unique(brain_region_pca)), fill=unique(brain_region_pca_color), cex = 0.5)

#plot to pdf

##### SVA ####

# add diagnosis in

Frontal_Pole_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Frontal_Pole_case_exprs_good_probes)))
Precentral_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Precentral_Gyrus_case_exprs_good_probes)))
Inferior_Frontal_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Inferior_Frontal_Gyrus_case_exprs_good_probes)))
Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes)))
Superior_Parietal_Lobule_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Superior_Parietal_Lobule_case_exprs_good_probes)))
Prefrontal_Cortex_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Prefrontal_Cortex_case_exprs_good_probes)))
Parahippocampal_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Parahippocampal_Gyrus_case_exprs_good_probes)))
Hippocampus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Hippocampus_case_exprs_good_probes)))
Inferior_Temporal_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Inferior_Temporal_Gyrus_case_exprs_good_probes)))
Middle_Temporal_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Middle_Temporal_Gyrus_case_exprs_good_probes)))
Superior_Temporal_Gyrus_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Superior_Temporal_Gyrus_case_exprs_good_probes)))
Temporal_Pole_case_exprs_good_probes<-cbind(Diagnosis="AD", as.data.frame(t(Temporal_Pole_case_exprs_good_probes)))

Frontal_Pole_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Frontal_Pole_control_exprs_good_probes)))
Precentral_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Precentral_Gyrus_control_exprs_good_probes)))
Inferior_Frontal_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Inferior_Frontal_Gyrus_control_exprs_good_probes)))
Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes)))
Superior_Parietal_Lobule_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Superior_Parietal_Lobule_control_exprs_good_probes)))
Prefrontal_Cortex_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Prefrontal_Cortex_control_exprs_good_probes)))
Parahippocampal_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Parahippocampal_Gyrus_control_exprs_good_probes)))
Hippocampus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Hippocampus_control_exprs_good_probes)))
Inferior_Temporal_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Inferior_Temporal_Gyrus_control_exprs_good_probes)))
Middle_Temporal_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Middle_Temporal_Gyrus_control_exprs_good_probes)))
Superior_Temporal_Gyrus_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Superior_Temporal_Gyrus_control_exprs_good_probes)))
Temporal_Pole_control_exprs_good_probes<-cbind(Diagnosis="CONTROL", as.data.frame(t(Temporal_Pole_control_exprs_good_probes)))

# check dataframe

head(Frontal_Pole_case_exprs_good_probes)[1:5]
head(Precentral_Gyrus_case_exprs_good_probes)[1:5]
head(Inferior_Frontal_Gyrus_case_exprs_good_probes)[1:5]
head(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes)[1:5]
head(Superior_Parietal_Lobule_case_exprs_good_probes)[1:5]
head(Prefrontal_Cortex_case_exprs_good_probes)[1:5]
head(Parahippocampal_Gyrus_case_exprs_good_probes)[1:5]
head(Hippocampus_case_exprs_good_probes)[1:5]
head(Inferior_Temporal_Gyrus_case_exprs_good_probes)[1:5]
head(Middle_Temporal_Gyrus_case_exprs_good_probes)[1:5]
head(Superior_Temporal_Gyrus_case_exprs_good_probes)[1:5]
head(Temporal_Pole_case_exprs_good_probes)[1:5]

head(Frontal_Pole_control_exprs_good_probes)[1:5]
head(Precentral_Gyrus_control_exprs_good_probes)[1:5]
head(Inferior_Frontal_Gyrus_control_exprs_good_probes)[1:5]
head(Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes)[1:5]
head(Superior_Parietal_Lobule_control_exprs_good_probes)[1:5]
head(Prefrontal_Cortex_control_exprs_good_probes)[1:5]
head(Parahippocampal_Gyrus_control_exprs_good_probes)[1:5]
head(Hippocampus_control_exprs_good_probes)[1:5]
head(Inferior_Temporal_Gyrus_control_exprs_good_probes)[1:5]
head(Middle_Temporal_Gyrus_control_exprs_good_probes)[1:5]
head(Superior_Temporal_Gyrus_control_exprs_good_probes)[1:5]
head(Temporal_Pole_control_exprs_good_probes)[1:5]

# check colnames same in case and control datasets
any(colnames(Frontal_Pole_case_exprs_good_probes)==colnames(Frontal_Pole_control_exprs_good_probes))==F
any(colnames(Precentral_Gyrus_case_exprs_good_probes)==colnames(Precentral_Gyrus_control_exprs_good_probes))==F
any(colnames(Inferior_Frontal_Gyrus_case_exprs_good_probes)==colnames(Inferior_Frontal_Gyrus_control_exprs_good_probes))==F
any(colnames(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes)==colnames(Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes))==F
any(colnames(Superior_Parietal_Lobule_case_exprs_good_probes)==colnames(Superior_Parietal_Lobule_control_exprs_good_probes))==F
any(colnames(Prefrontal_Cortex_case_exprs_good_probes)==colnames(Prefrontal_Cortex_control_exprs_good_probes))==F
any(colnames(Parahippocampal_Gyrus_case_exprs_good_probes)==colnames(Parahippocampal_Gyrus_control_exprs_good_probes))==F
any(colnames(Hippocampus_case_exprs_good_probes)==colnames(Hippocampus_control_exprs_good_probes))==F
any(colnames(Inferior_Temporal_Gyrus_case_exprs_good_probes)==colnames(Inferior_Temporal_Gyrus_control_exprs_good_probes))==F
any(colnames(Middle_Temporal_Gyrus_case_exprs_good_probes)==colnames(Middle_Temporal_Gyrus_control_exprs_good_probes))==F
any(colnames(Superior_Temporal_Gyrus_case_exprs_good_probes)==colnames(Superior_Temporal_Gyrus_control_exprs_good_probes))==F
any(colnames(Temporal_Pole_case_exprs_good_probes)==colnames(Temporal_Pole_control_exprs_good_probes))==F

# rbind datasets togther

Frontal_Pole_exprs_good_probes<-rbind(Frontal_Pole_case_exprs_good_probes, Frontal_Pole_control_exprs_good_probes)
Precentral_Gyrus_exprs_good_probes<-rbind(Precentral_Gyrus_case_exprs_good_probes, Precentral_Gyrus_control_exprs_good_probes)
Inferior_Frontal_Gyrus_exprs_good_probes<-rbind(Inferior_Frontal_Gyrus_case_exprs_good_probes, Inferior_Frontal_Gyrus_control_exprs_good_probes)
Dorsolateral_Prefrontal_Cortex_exprs_good_probes<-rbind(Dorsolateral_Prefrontal_Cortex_case_exprs_good_probes, Dorsolateral_Prefrontal_Cortex_control_exprs_good_probes)
Superior_Parietal_Lobule_exprs_good_probes<-rbind(Superior_Parietal_Lobule_case_exprs_good_probes, Superior_Parietal_Lobule_control_exprs_good_probes)
Prefrontal_Cortex_exprs_good_probes<-rbind(Prefrontal_Cortex_case_exprs_good_probes, Prefrontal_Cortex_control_exprs_good_probes)
Parahippocampal_Gyrus_exprs_good_probes<-rbind(Parahippocampal_Gyrus_case_exprs_good_probes, Parahippocampal_Gyrus_control_exprs_good_probes)
Hippocampus_exprs_good_probes<-rbind(Hippocampus_case_exprs_good_probes, Hippocampus_control_exprs_good_probes)
Inferior_Temporal_Gyrus_exprs_good_probes<-rbind(Inferior_Temporal_Gyrus_case_exprs_good_probes, Inferior_Temporal_Gyrus_control_exprs_good_probes)
Middle_Temporal_Gyrus_exprs_good_probes<-rbind(Middle_Temporal_Gyrus_case_exprs_good_probes, Middle_Temporal_Gyrus_control_exprs_good_probes)
Superior_Temporal_Gyrus_exprs_good_probes<-rbind(Superior_Temporal_Gyrus_case_exprs_good_probes, Superior_Temporal_Gyrus_control_exprs_good_probes)
Temporal_Pole_exprs_good_probes<-rbind(Temporal_Pole_case_exprs_good_probes, Temporal_Pole_control_exprs_good_probes)

# check dataframe

head(Frontal_Pole_exprs_good_probes)[1:5]
dim(Frontal_Pole_exprs_good_probes)

head(Precentral_Gyrus_exprs_good_probes)[1:5]
dim(Precentral_Gyrus_exprs_good_probes)

head(Inferior_Frontal_Gyrus_exprs_good_probes)[1:5]
dim(Inferior_Frontal_Gyrus_exprs_good_probes)

head(Dorsolateral_Prefrontal_Cortex_exprs_good_probes)[1:5]
dim(Dorsolateral_Prefrontal_Cortex_exprs_good_probes)

head(Superior_Parietal_Lobule_exprs_good_probes)[1:5]
dim(Superior_Parietal_Lobule_exprs_good_probes)

head(Prefrontal_Cortex_exprs_good_probes)[1:5]
dim(Prefrontal_Cortex_exprs_good_probes)

head(Parahippocampal_Gyrus_exprs_good_probes)[1:5]
dim(Parahippocampal_Gyrus_exprs_good_probes)

head(Hippocampus_exprs_good_probes)[1:5]
dim(Hippocampus_exprs_good_probes)

head(Inferior_Temporal_Gyrus_exprs_good_probes)[1:5]
dim(Inferior_Temporal_Gyrus_exprs_good_probes)

head(Middle_Temporal_Gyrus_exprs_good_probes)[1:5]
dim(Middle_Temporal_Gyrus_exprs_good_probes)

head(Superior_Temporal_Gyrus_exprs_good_probes)[1:5]
dim(Superior_Temporal_Gyrus_exprs_good_probes)

head(Temporal_Pole_exprs_good_probes)[1:5]
dim(Temporal_Pole_exprs_good_probes)

# add gender in 

Frontal_Pole_exprs_good_probes<-merge(gender_comparison[1], Frontal_Pole_exprs_good_probes, by="row.names")
rownames(Frontal_Pole_exprs_good_probes)<-Frontal_Pole_exprs_good_probes$Row.names
Frontal_Pole_exprs_good_probes$Row.names<-NULL

Precentral_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Precentral_Gyrus_exprs_good_probes, by="row.names")
rownames(Precentral_Gyrus_exprs_good_probes)<-Precentral_Gyrus_exprs_good_probes$Row.names
Precentral_Gyrus_exprs_good_probes$Row.names<-NULL

Inferior_Frontal_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Inferior_Frontal_Gyrus_exprs_good_probes, by="row.names")
rownames(Inferior_Frontal_Gyrus_exprs_good_probes)<-Inferior_Frontal_Gyrus_exprs_good_probes$Row.names
Inferior_Frontal_Gyrus_exprs_good_probes$Row.names<-NULL

Dorsolateral_Prefrontal_Cortex_exprs_good_probes<-merge(gender_comparison[1], Dorsolateral_Prefrontal_Cortex_exprs_good_probes, by="row.names")
rownames(Dorsolateral_Prefrontal_Cortex_exprs_good_probes)<-Dorsolateral_Prefrontal_Cortex_exprs_good_probes$Row.names
Dorsolateral_Prefrontal_Cortex_exprs_good_probes$Row.names<-NULL

Superior_Parietal_Lobule_exprs_good_probes<-merge(gender_comparison[1], Superior_Parietal_Lobule_exprs_good_probes, by="row.names")
rownames(Superior_Parietal_Lobule_exprs_good_probes)<-Superior_Parietal_Lobule_exprs_good_probes$Row.names
Superior_Parietal_Lobule_exprs_good_probes$Row.names<-NULL

Prefrontal_Cortex_exprs_good_probes<-merge(gender_comparison[1], Prefrontal_Cortex_exprs_good_probes, by="row.names")
rownames(Prefrontal_Cortex_exprs_good_probes)<-Prefrontal_Cortex_exprs_good_probes$Row.names
Prefrontal_Cortex_exprs_good_probes$Row.names<-NULL

Parahippocampal_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Parahippocampal_Gyrus_exprs_good_probes, by="row.names")
rownames(Parahippocampal_Gyrus_exprs_good_probes)<-Parahippocampal_Gyrus_exprs_good_probes$Row.names
Parahippocampal_Gyrus_exprs_good_probes$Row.names<-NULL

Hippocampus_exprs_good_probes<-merge(gender_comparison[1], Hippocampus_exprs_good_probes, by="row.names")
rownames(Hippocampus_exprs_good_probes)<-Hippocampus_exprs_good_probes$Row.names
Hippocampus_exprs_good_probes$Row.names<-NULL

Inferior_Temporal_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Inferior_Temporal_Gyrus_exprs_good_probes, by="row.names")
rownames(Inferior_Temporal_Gyrus_exprs_good_probes)<-Inferior_Temporal_Gyrus_exprs_good_probes$Row.names
Inferior_Temporal_Gyrus_exprs_good_probes$Row.names<-NULL

Middle_Temporal_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Middle_Temporal_Gyrus_exprs_good_probes, by="row.names")
rownames(Middle_Temporal_Gyrus_exprs_good_probes)<-Middle_Temporal_Gyrus_exprs_good_probes$Row.names
Middle_Temporal_Gyrus_exprs_good_probes$Row.names<-NULL

Superior_Temporal_Gyrus_exprs_good_probes<-merge(gender_comparison[1], Superior_Temporal_Gyrus_exprs_good_probes, by="row.names")
rownames(Superior_Temporal_Gyrus_exprs_good_probes)<-Superior_Temporal_Gyrus_exprs_good_probes$Row.names
Superior_Temporal_Gyrus_exprs_good_probes$Row.names<-NULL

Temporal_Pole_exprs_good_probes<-merge(gender_comparison[1], Temporal_Pole_exprs_good_probes, by="row.names")
rownames(Temporal_Pole_exprs_good_probes)<-Temporal_Pole_exprs_good_probes$Row.names
Temporal_Pole_exprs_good_probes$Row.names<-NULL

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

# check SV in data

check_SV_in_data(Frontal_Pole_exprs_good_probes)
check_SV_in_data(Precentral_Gyrus_exprs_good_probes)
check_SV_in_data(Inferior_Frontal_Gyrus_exprs_good_probes)
check_SV_in_data(Dorsolateral_Prefrontal_Cortex_exprs_good_probes)
check_SV_in_data(Superior_Parietal_Lobule_exprs_good_probes)
check_SV_in_data(Prefrontal_Cortex_exprs_good_probes)
check_SV_in_data(Parahippocampal_Gyrus_exprs_good_probes)
check_SV_in_data(Hippocampus_exprs_good_probes)
check_SV_in_data(Inferior_Temporal_Gyrus_exprs_good_probes)
check_SV_in_data(Middle_Temporal_Gyrus_exprs_good_probes)
check_SV_in_data(Superior_Temporal_Gyrus_exprs_good_probes)
check_SV_in_data(Temporal_Pole_exprs_good_probes)

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

# apply function

Frontal_Pole_exprs_good_probes_sva_adjusted<-adjust_for_sva(Frontal_Pole_exprs_good_probes)
Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted<-adjust_for_sva(Inferior_Frontal_Gyrus_exprs_good_probes)
Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted<-adjust_for_sva(Dorsolateral_Prefrontal_Cortex_exprs_good_probes)
Prefrontal_Cortex_exprs_good_probes_sva_adjusted<-adjust_for_sva(Prefrontal_Cortex_exprs_good_probes)
Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted<-adjust_for_sva(Inferior_Temporal_Gyrus_exprs_good_probes)


# changed name from non sva adjusted data

Precentral_Gyrus_exprs_good_probes_sva_adjusted<-Precentral_Gyrus_exprs_good_probes
Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted<-Superior_Parietal_Lobule_exprs_good_probes
Hippocampus_exprs_good_probes_sva_adjusted<-Hippocampus_exprs_good_probes
Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted<-Middle_Temporal_Gyrus_exprs_good_probes
Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted<-Superior_Temporal_Gyrus_exprs_good_probes
Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted<-Parahippocampal_Gyrus_exprs_good_probes
Temporal_Pole_exprs_good_probes_sva_adjusted<-Temporal_Pole_exprs_good_probes

# re check SV

check_SV_in_data(Frontal_Pole_exprs_good_probes_sva_adjusted)
check_SV_in_data(Precentral_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted)
check_SV_in_data(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted)
check_SV_in_data(Prefrontal_Cortex_exprs_good_probes_sva_adjusted)
check_SV_in_data(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Hippocampus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted)
check_SV_in_data(Temporal_Pole_exprs_good_probes_sva_adjusted)

##### SAMPLE NETWORK PLOT #####

# separate case and control for sample netwrok analysis

Frontal_Pole_exprs_good_probes_sva_adjusted_case<-Frontal_Pole_exprs_good_probes_sva_adjusted[Frontal_Pole_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Precentral_Gyrus_exprs_good_probes_sva_adjusted_case<-Precentral_Gyrus_exprs_good_probes_sva_adjusted[Precentral_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case<-Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted[Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case<-Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted[Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case<-Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted[Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case<-Prefrontal_Cortex_exprs_good_probes_sva_adjusted[Prefrontal_Cortex_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case<-Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted[Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Hippocampus_exprs_good_probes_sva_adjusted_case<-Hippocampus_exprs_good_probes_sva_adjusted[Hippocampus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case<-Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case<-Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case<-Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Temporal_Pole_exprs_good_probes_sva_adjusted_case<-Temporal_Pole_exprs_good_probes_sva_adjusted[Temporal_Pole_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]

Frontal_Pole_exprs_good_probes_sva_adjusted_control<-Frontal_Pole_exprs_good_probes_sva_adjusted[Frontal_Pole_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Precentral_Gyrus_exprs_good_probes_sva_adjusted_control<-Precentral_Gyrus_exprs_good_probes_sva_adjusted[Precentral_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control<-Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted[Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control<-Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted[Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control<-Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted[Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control<-Prefrontal_Cortex_exprs_good_probes_sva_adjusted[Prefrontal_Cortex_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control<-Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted[Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Hippocampus_exprs_good_probes_sva_adjusted_control<-Hippocampus_exprs_good_probes_sva_adjusted[Hippocampus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control<-Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control<-Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control<-Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted[Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Temporal_Pole_exprs_good_probes_sva_adjusted_control<-Temporal_Pole_exprs_good_probes_sva_adjusted[Temporal_Pole_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]

table(Frontal_Pole_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Precentral_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Hippocampus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case$Diagnosis)
table(Temporal_Pole_exprs_good_probes_sva_adjusted_case$Diagnosis)

table(Frontal_Pole_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Precentral_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Hippocampus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control$Diagnosis)
table(Temporal_Pole_exprs_good_probes_sva_adjusted_control$Diagnosis)

# sample plot function - taken from steve expression pipeline

sampleNetwork_plot <- function(dataset) {
  datExprs<-t(dataset[3:dim(dataset)[2]])
  diagnosis<-dataset[2]
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
  datExprs<-t(dataset[3:dim(dataset)[2]])
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

# apply function to cases

Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Frontal_Pole_exprs_good_probes_sva_adjusted_case, -3)
Precentral_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Precentral_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case, -3)
Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case, -3)
Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case, -3)
Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
Hippocampus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Hippocampus_exprs_good_probes_sva_adjusted_case, -2)
Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -1.96)
Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -2)
Temporal_Pole_exprs_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Temporal_Pole_exprs_good_probes_sva_adjusted_case, -3)

#apply function to controls

Frontal_Pole_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Frontal_Pole_exprs_good_probes_sva_adjusted_control, -3)
Precentral_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Precentral_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control, -3)
Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control, -3)
Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control, -3)
Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Hippocampus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Hippocampus_exprs_good_probes_sva_adjusted_control, -3)
Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
Temporal_Pole_exprs_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Temporal_Pole_exprs_good_probes_sva_adjusted_control, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("Frontal_Pole_sample_network_analysis_cases.pdf")
run_sample_network_plot(Frontal_Pole_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Precentral_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Precentral_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Inferior_Frontal_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Dorsolateral_Prefrontal_Cortex_sample_network_analysis_cases.pdf")
run_sample_network_plot(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Superior_Parietal_Lobule_sample_network_analysis_cases.pdf")
run_sample_network_plot(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Prefrontal_Cortex_sample_network_analysis_cases.pdf")
run_sample_network_plot(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Parahippocampal_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Hippocampus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Hippocampus_exprs_good_probes_sva_adjusted_case, -2)
dev.off()

pdf("Inferior_Temporal_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Middle_Temporal_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -1.96)
dev.off()

pdf("Superior_Temporal_Gyrus_sample_network_analysis_cases.pdf")
run_sample_network_plot(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case, -2)
dev.off()

pdf("Temporal_Pole_sample_network_analysis_cases.pdf")
run_sample_network_plot(Temporal_Pole_exprs_good_probes_sva_adjusted_case, -3)
dev.off()


pdf("Frontal_Pole_sample_network_analysis_controls.pdf")
run_sample_network_plot(Frontal_Pole_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Precentral_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Precentral_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Inferior_Frontal_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Dorsolateral_Prefrontal_Cortex_sample_network_analysis_controls.pdf")
run_sample_network_plot(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Superior_Parietal_Lobule_sample_network_analysis_controls.pdf")
run_sample_network_plot(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Prefrontal_Cortex_sample_network_analysis_controls.pdf")
run_sample_network_plot(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Parahippocampal_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Hippocampus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Hippocampus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Inferior_Temporal_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Middle_Temporal_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Superior_Temporal_Gyrus_sample_network_analysis_controls.pdf")
run_sample_network_plot(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Temporal_Pole_sample_network_analysis_controls.pdf")
run_sample_network_plot(Temporal_Pole_exprs_good_probes_sva_adjusted_control, -3)
dev.off()

setwd(work_dir)

##### CREATE QC'd DATASET #####

# check colnames in dataframes

any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Precentral_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Hippocampus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Temporal_Pole_exprs_good_probes_sva_adjusted_case_QC))==F

any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Precentral_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Hippocampus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC)==colnames(Temporal_Pole_exprs_good_probes_sva_adjusted_control_QC))==F

# rbind all QC dataframes together

brain_data_QCd<-rbind(Frontal_Pole_exprs_good_probes_sva_adjusted_case_QC,
                      Precentral_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC,
                      Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_case_QC,
                      Prefrontal_Cortex_exprs_good_probes_sva_adjusted_case_QC,
                      Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Hippocampus_exprs_good_probes_sva_adjusted_case_QC,
                      Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_case_QC,
                      Temporal_Pole_exprs_good_probes_sva_adjusted_case_QC,
                      Frontal_Pole_exprs_good_probes_sva_adjusted_control_QC,
                      Precentral_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Inferior_Frontal_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Dorsolateral_Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC,
                      Superior_Parietal_Lobule_exprs_good_probes_sva_adjusted_control_QC,
                      Prefrontal_Cortex_exprs_good_probes_sva_adjusted_control_QC,
                      Parahippocampal_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Hippocampus_exprs_good_probes_sva_adjusted_control_QC,
                      Inferior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Middle_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Superior_Temporal_Gyrus_exprs_good_probes_sva_adjusted_control_QC,
                      Temporal_Pole_exprs_good_probes_sva_adjusted_control_QC)

dim(brain_data_normalised_exprs_good_probes)
dim(brain_data_QCd)

head(brain_data_QCd)[1:5]

##### PCA ON CLEAN DATA ####

# calculate pca

pca_QCd<-prcomp(t(brain_data_QCd[3:dim(brain_data_QCd)[2]]))

# plot variance
plot(pca_QCd, type="l")

# sumary of pca
summary_pca_QCd<-summary(pca_QCd)

# pca importance
pca_importance_var_exp_QCd<-summary_pca_QCd$importance[2,]
pca_importance_var_exp_cum_QCd<-summary_pca_QCd$importance[3,]

par(mfrow=c(1,2))

#plot variance explained
plot(pca_importance_var_exp_QCd, ylab="PCA Proportion of Variance Explained after QC", type="b", col="blue")

#plot variance explained cumalative
plot(pca_importance_var_exp_cum_QCd, ylab="PCA Cumulative Proportion of Variance Explained after QC", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)

par(dev.off)
dev.off()

# pca matrix plot

#color

# order of samples in expression data
Diagnosis_temp_QCd<-rownames(brain_data_QCd)
Diagnosis_temp_QCd

# keep only samples in Diagnosis_lookup that are in QC'd data
Diagnosis_lookup_QCd<-subset(Diagnosis_lookup, individual_ID %in% Diagnosis_temp_QCd)

# match order
Diagnosis_pca_QCd<-Diagnosis_lookup_QCd[match(Diagnosis_temp_QCd, Diagnosis_lookup_QCd$individual_ID),]

head(Diagnosis_pca_QCd)

# assign color to group
Diagnosis_pca_color_QCd<-labels2colors(as.character(Diagnosis_pca_QCd$Diagnosis))

# pca plot - color by disease - case/control
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Disease Status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomleft', unique(Diagnosis_pca_QCd$Diagnosis), fill=unique(Diagnosis_pca_color_QCd))

#pca plot - color by brain region

# order of samples in expression data
brain_region_temp_QCd<-rownames(brain_data_QCd)
brain_region_temp_QCd

# match order
brain_region_pca_QCd<-brain_region_lookup[match(brain_region_temp_QCd, rownames(brain_region_lookup)),]

# assign color to group
brain_region_pca_color_QCd<-labels2colors(as.character(brain_region_pca_QCd))

# pca plot - color by disease - case/control
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Brain Region after QC",col="black", pch=21,bg=brain_region_pca_color_QCd)
legend('bottomleft', as.character(unique(brain_region_pca_QCd)), fill=unique(brain_region_pca_color_QCd), cex=0.7)

#plot to pdf - before/after QC

setwd(pca_dir)

pdf("PCA_plot_prior_and_after_QC.pdf")
# pca plot by diagnosis
plot(pca$rotation[,1:2], main=" PCA plot coloured by Disease Status",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomright', unique(Diagnosis_pca$Diagnosis), fill=unique(Diagnosis_pca_color))
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Disease Status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomleft', unique(Diagnosis_pca_QCd$Diagnosis), fill=unique(Diagnosis_pca_color_QCd))
#pca plot by brain region
plot(pca$rotation[,1:2], main=" PCA plot coloured by Brain Region",col="black", pch=21,bg=brain_region_pca_color)
legend('bottomright', as.character(unique(brain_region_pca)), fill=unique(brain_region_pca_color), cex = 0.5)
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Brain Region after QC",col="black", pch=21,bg=brain_region_pca_color_QCd)
legend('bottomleft', as.character(unique(brain_region_pca_QCd)), fill=unique(brain_region_pca_color_QCd), cex=0.7)
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

brain_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(brain_data_QCd[3:dim(brain_data_QCd)[2]], hgu133a.db_mapping)
dim(brain_data_QCd)
dim(brain_data_QCd_entrez_id)
length(which(duplicated(colnames(brain_data_QCd_entrez_id))))

head(brain_data_QCd_entrez_id[100:110])

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

## apply function to dataframes - check number of probes - check duplicates

brain_data_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(brain_data_QCd_entrez_id)
dim(brain_data_QCd_entrez_id_unique)
length(which(duplicated(colnames(brain_data_QCd_entrez_id_unique))))

head(brain_data_QCd_entrez_id_unique[1:5])

##### ATTACH DIAGNOSIS AND BRAIN REGION #####

# "Tissue", "MMSE", "Clinical_Gender", "Predicted_Gender", Diagnosis_sub_cat", "Diagnosis", "BRAAK", "APOE", "Age"

head(phenotype_data_brain_region)

# read in APOE info

APOE_info<-read.table("/media/hamel/1TB/Projects/Brain_expression/old/1.old/AD/AMP_MSBB/Phenotype_info/AMP-AD_MSBB_MSSM_array_ApoE_genotype.txt", header = T)
APOE_info$APOE<-paste(APOE_info$ApoE1, APOE_info$ApoE2, sep=",")
APOE_info$ApoE1<-NULL
APOE_info$ApoE2<-NULL
colnames(APOE_info)[1]<-"BB_ID"

# merge APOE infor with pheno

phenotype<-merge(phenotype_data_brain_region, APOE_info, by="BB_ID")

head(phenotype)
dim(phenotype)

# merge clinial/predicted gender

head(gender_comparison)
gender_import<-gender_comparison
gender_import$individual_ID<-rownames(gender_import)
head(gender_import)

phenotype<-merge(phenotype, gender_import, by="individual_ID")

head(phenotype)

#subset phenotype to samples in brain data

dim(brain_data_QCd_entrez_id_unique)
dim(phenotype)

phenotype<-subset(phenotype, phenotype$individual_ID %in% rownames(brain_data_QCd_entrez_id_unique))
phenotype<-subset(phenotype, phenotype$ARRAY=="AffymetrixHG-U133B")

dim(brain_data_QCd_entrez_id_unique)
dim(phenotype)

# keep only "Tissue", "MMSE", "Gender", "Diagnosis_sub_cat", "Diagnosis", "BRAAK", "APOE", "Age" columns

rownames(phenotype)<-phenotype$individual_ID

head(phenotype)
colnames(phenotype)

phenotype_to_attach<-phenotype[c(6, 20, 21, 18, 13, 19, 11)]
head(phenotype_to_attach)

phenotype_to_attach$MMSE<-"Unknown"
phenotype_to_attach$Diagnosis_sub_cat<-"Unknown"

head(phenotype_to_attach)

colnames(phenotype_to_attach)<-c("Tissue", "Clinical_Gender", "Predicted_Gender", "Diagnosis", "BRAAK", "APOE", "Age", "MMSE", "Diagnosis_sub_cat")
head(phenotype_to_attach)

# create sub cats based on BRAAK/MMSE scores - use BRAAK- Control = 0-1, Incipient= 2-3, Moderate=4, Severe= 5-6
table(phenotype_to_attach$BRAAK, exclude=NULL)

phenotype_to_attach[phenotype_to_attach$BRAAK=="0" | phenotype_to_attach$BRAAK=="1", 9]<-"Control"
phenotype_to_attach[phenotype_to_attach$BRAAK=="2" | phenotype_to_attach$BRAAK=="3", 9]<-"Incipient_AD"
phenotype_to_attach[phenotype_to_attach$BRAAK=="4", 9]<-"Moderate_AD"
phenotype_to_attach[phenotype_to_attach$BRAAK=="5" | phenotype_to_attach$BRAAK=="6", 9]<-"Severe_AD"

table(phenotype_to_attach$Diagnosis_sub_cat, exclude=NULL)

# order pheno information

phenotype_to_attach<-phenotype_to_attach[,order(names(phenotype_to_attach), decreasing = T)]

# merge pheno info with expr

brain_data_QCd_entrez_id_unique_pheno<-merge(phenotype_to_attach, brain_data_QCd_entrez_id_unique, by="row.names", all.x=T)

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

rownames(brain_data_QCd_entrez_id_unique_pheno)<-brain_data_QCd_entrez_id_unique_pheno$Row.names

brain_data_QCd_entrez_id_unique_pheno$Row.names<-NULL

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

dim(brain_data_QCd_entrez_id_unique)

dim(phenotype_to_attach)

dim(brain_data_QCd_entrez_id_unique_pheno)

# lost 7 samples as no pheno available

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(brain_data_QCd_entrez_id_unique_pheno, file="AMP_MSBB_U133A_pre-processed_data2.txt", sep="\t")

##### SAVE PHENOTYPE DATAFRAME #####

setwd(clean_data_dir)

write.table(phenotype, file="AMP_MSBB_U133A_phenotype_data2.txt", sep="\t")

##### WRITE LIST OF CONTROL GENES EXPRESSED #####

setwd("/media/hamel/1TB/Projects/Brain_expression/2.Expression_across_brain_regions_in_control_datasets/1.Data")

write(Frontal_Pole_control_expressed_probes_list, file="AMP_MSBB_U133A_Frontal_Pole.txt") 
write(Precentral_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Precentral_Gyrus.txt") 
write(Inferior_Frontal_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Inferior_Frontal_Gyrus.txt") 
write(Dorsolateral_Prefrontal_Cortex_control_expressed_probes_list, file="AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex.txt") 
write(Superior_Parietal_Lobule_control_expressed_probes_list, file="AMP_MSBB_U133A_Superior_Parietal_Lobule.txt") 
write(Prefrontal_Cortex_control_expressed_probes_list, file="AMP_MSBB_U133A_Prefrontal_Cortex.txt") 
write(Parahippocampal_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Parahippocampal_Gyrus.txt") 
write(Hippocampus_control_expressed_probes_list, file="AMP_MSBB_U133A_Hippocampus.txt") 
write(Inferior_Temporal_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Inferior_Temporal_Gyrus.txt") 
write(Middle_Temporal_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Middle_Temporal_Gyrus.txt") 
write(Superior_Temporal_Gyrus_control_expressed_probes_list, file="AMP_MSBB_U133A_Superior_Temporal_Gyrus.txt") 
write(Temporal_Pole_control_expressed_probes_list, file="AMP_MSBB_U133A_Temporal_Pole.txt") 

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="AMP_MSBB_U133A_data_processing2.Rdata")
