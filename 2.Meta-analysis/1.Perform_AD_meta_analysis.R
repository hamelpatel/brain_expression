########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                                       META-ANALYSIS                                                  #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 11/08/2016

##### DESCRIPTION OF ANALYSIS ####
## brain expression datasets pre-processed - using clean data produced in DE analysis -i.e AD with <=3 BRAAK removed and control with 3>=BRAAK removed
## running MetaOmics - for DE and pathway analysis
##
##
## NOTE: latest gene set files files can be downloaded from http://software.broadinstitute.org/gsea/downloads.jsp
## Temporal lobe had to be partially run on cluster due to lack of ram. - only run with ind.tail="high" as "low" give 0.99 cor
## NO GENE ENTREZ ID CONVERSION TO GENE SYMBOL AS PACKAGE HAS ISSUES
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/Clean_Data"

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/4.AD_DE_Meta_Analysis/"

setwd(work_dir)

# create dir DE and patway analysis

# DE directory
dir.create(paste(work_dir,"Meta_DE", sep="/"))
Meta_DE_dir=paste(work_dir,"Meta_DE", sep="/")

# Pathway analysis directory
dir.create(paste(work_dir,"Meta_Pathway", sep="/"))
Meta_Pathway_dir=paste(work_dir,"Meta_Pathway", sep="/")

# QC plots
dir.create(paste(work_dir,"MetaQC_plots", sep="/"))
MetaQC_plots_dir=paste(work_dir,"MetaQC_plots", sep="/")

# data in metomics format
dir.create(paste(work_dir, "Data_in_metaomics_format", sep="/"))
Data_in_metaomics_format=paste(work_dir, "Data_in_metaomics_format", sep="/")

# DE heatmaps
dir.create(paste(Meta_DE_dir,"DE_Heatmaps", sep="/"))
DE_Heamaps_dir=paste(Meta_DE_dir,"DE_Heatmaps", sep="/")

# AW results
dir.create(paste(work_dir,"AW_results", sep="/"))
AW_results_dir=paste(work_dir,"AW_results", sep="/")

# AW results - with logfc
dir.create(paste(AW_results_dir,"AW_results_with_logFC", sep="/"))
AW_results_logfc_dir=paste(AW_results_dir,"AW_results_with_logFC", sep="/")

# directories for all brain regions
dir.create(paste(work_dir,"Temporal_Lobe", sep="/"))
Temporal_Lobe_dir=paste(work_dir,"Temporal_Lobe", sep="/")

dir.create(paste(work_dir,"Frontal_Lobe", sep="/"))
Frontal_Lobe_dir=paste(work_dir,"Frontal_Lobe", sep="/")

dir.create(paste(work_dir,"Parietal_Lobe", sep="/"))
Parietal_Lobe_dir=paste(work_dir,"Parietal_Lobe", sep="/")

dir.create(paste(work_dir,"Cerebellum", sep="/"))
Cerebellum_dir=paste(work_dir,"Cerebellum", sep="/")

dir.create(paste(work_dir,"All_brain_regions_combined", sep="/"))
All_brain_regions_combined_dir=paste(work_dir,"All_brain_regions_combined", sep="/")

##### LOAD LIBRARIES ####

# install.packages("MetaQC", dependencies=T)
# install.packages("MetaDE", dependencies=T)
# install.packages("MetaPath", dependencies=T)
# install.packages("MetaPCA", dependencies=T)

library(MetaQC)
library(MetaDE)
library(MetaPath)
library(MetaPCA)
library(org.Hs.eg.db)
library(dplyr)
library(knitr)
library(GSEABase)
library(rmeta)

##### READ DATA - 11 datasets #####

setwd(data_dir)

E_GEOD_28146_clean <- read.table(file="E_GEOD_28146_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_29378_clean <- read.table(file="E_GEOD_29378_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_36980_clean <- read.table(file="E_GEOD_36980_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_48350_clean <- read.table(file="E_GEOD_48350_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_1297_clean <- read.table(file="E_GEOD_1297_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_5281_clean <- read.table(file="E_GEOD_5281_clean.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MAYOeGWAS_Temporal_Cortex_clean <- read.table(file="AMP_MAYOeGWAS_Temporal_Cortex_clean", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MAYOeGWAS_Cerebellum_clean <- read.table(file="AMP_MAYOeGWAS_Cerebellum_clean", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MSBB_U133A_clean <- read.table(file="AMP_MSBB_U133A_clean", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MSBB_U133B_clean <- read.table(file="AMP_MSBB_U133B_clean", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
inhouse_data_clean <- read.table(file="inhouse_data_clean", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)

setwd(work_dir)

head(E_GEOD_28146_clean)[1:10]
head(E_GEOD_29378_clean)[1:10]
head(E_GEOD_36980_clean)[1:10]
head(E_GEOD_48350_clean)[1:10]
head(E_GEOD_1297_clean)[1:10]
head(E_GEOD_5281_clean)[1:10]
head(inhouse_data_clean)[1:10]
head(AMP_MSBB_U133A_clean)[1:10]
head(AMP_MSBB_U133B_clean)[1:10]
head(AMP_MAYOeGWAS_Temporal_Cortex_clean)[1:10]
head(AMP_MAYOeGWAS_Cerebellum_clean)[1:10]

dim(E_GEOD_28146_clean)
dim(E_GEOD_29378_clean)
dim(E_GEOD_36980_clean)
dim(E_GEOD_48350_clean)
dim(E_GEOD_1297_clean)
dim(E_GEOD_5281_clean)
dim(inhouse_data_clean)
dim(AMP_MSBB_U133A_clean)
dim(AMP_MSBB_U133B_clean)
dim(AMP_MAYOeGWAS_Temporal_Cortex_clean)
dim(AMP_MAYOeGWAS_Cerebellum_clean)

##### SEPARATE DATAFRAMES BY BRAIN REGION - ADD BRAIN REGION TITLE TO DATAFRAME #####

# E_GEOD_28146

E_GEOD_28146_Hippocampus_CA1<-E_GEOD_28146_clean

# E_GEOD_29378

E_GEOD_29378_Hippocampus_CA1<-E_GEOD_29378_clean[E_GEOD_29378_clean$Tissue=="Hippocampus_CA1",]
E_GEOD_29378_Hippocampus_CA3<-E_GEOD_29378_clean[E_GEOD_29378_clean$Tissue=="Hippocampus_CA3",]

# E_GEOD_36980

E_GEOD_36980_Frontal_Cortex<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Frontal_Cortex",]
E_GEOD_36980_Hippocampus<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Hippocampus",]
E_GEOD_36980_Temporal_Cortex<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Temporal_Cortex",]

# E_GEOD_48350

E_GEOD_48350_Entorhinal_Cortex<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Entorhinal_Cortex",]
E_GEOD_48350_Hippocampus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Hippocampus",]
E_GEOD_48350_Postcentral_Gyrus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Postcentral_Gyrus",]
E_GEOD_48350_Superior_Frontal_Gyrus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Superior_Frontal_Gyrus",]

# E_GEOD_1297

E_GEOD_1297_Hippocampus<-E_GEOD_1297_clean

# E_GEOD_5281

E_GEOD_5281_Entorhinal_Cortex<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Entorhinal_Cortex",]
E_GEOD_5281_Hippocampus_CA1<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Hippocampus_CA1",]
E_GEOD_5281_Medial_Temporal_Gyrus<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Medial_Temporal_Gyrus",]
E_GEOD_5281_Superior_Frontal_Gyrus<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Superior_Frontal_Gyrus",]
E_GEOD_5281_Posterior_Singulate<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Posterior_Singulate",]
E_GEOD_5281_Primary_Visual_Cortex<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Primary_Visual_Cortex",]

# inhouse data

inhouse_data_Entorhinal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Entorhinal_Cortex",]
inhouse_data_Cerebellum<-inhouse_data_clean[inhouse_data_clean$Tissue=="Cerebellum",]
inhouse_data_Frontal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Frontal_Cortex",]
inhouse_data_Temporal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Temporal_Cortex",]

# AMP MAYOeGWAS

AMP_MAYOeGWAS_Temporal_Cortex<-AMP_MAYOeGWAS_Temporal_Cortex_clean
AMP_MAYOeGWAS_Cerebellum<-AMP_MAYOeGWAS_Cerebellum_clean

# AMP_MSBB_U133A

AMP_MSBB_U133A_Frontal_Pole<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Frontal_Pole",]
AMP_MSBB_U133A_Precentral_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Precentral_Gyrus",]
AMP_MSBB_U133A_Inferior_Frontal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Inferior_Frontal_Gyrus",]
AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Dorsolateral_Prefrontal_Cortex",]
AMP_MSBB_U133A_Superior_Parietal_Lobule<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Superior_Parietal_Lobule",]
AMP_MSBB_U133A_Prefrontal_Cortex<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Prefrontal_Cortex",]
AMP_MSBB_U133A_Parahippocampal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Parahippocampal_Gyrus",]
AMP_MSBB_U133A_Hippocampus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Hippocampus",]
AMP_MSBB_U133A_Inferior_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Inferior_Temporal_Gyrus",]
AMP_MSBB_U133A_Middle_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Middle_Temporal_Gyrus",]
AMP_MSBB_U133A_Superior_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Superior_Temporal_Gyrus",]
AMP_MSBB_U133A_Temporal_Pole<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Temporal_Pole",]

# AMP_MSBB_U133B

AMP_MSBB_U133B_Frontal_Pole<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Frontal_Pole",]
AMP_MSBB_U133B_Precentral_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Precentral_Gyrus",]
AMP_MSBB_U133B_Inferior_Frontal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Inferior_Frontal_Gyrus",]
AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Dorsolateral_Prefrontal_Cortex",]
AMP_MSBB_U133B_Superior_Parietal_Lobule<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Superior_Parietal_Lobule",]
AMP_MSBB_U133B_Prefrontal_Cortex<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Prefrontal_Cortex",]
AMP_MSBB_U133B_Parahippocampal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Parahippocampal_Gyrus",]
AMP_MSBB_U133B_Hippocampus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Hippocampus",]
AMP_MSBB_U133B_Inferior_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Inferior_Temporal_Gyrus",]
AMP_MSBB_U133B_Middle_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Middle_Temporal_Gyrus",]
AMP_MSBB_U133B_Superior_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Superior_Temporal_Gyrus",]
AMP_MSBB_U133B_Temporal_Pole<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Temporal_Pole",]

##### GROUP BRAIN REGIONS - total 47 datasets #####

#Frontal_Lobe
head(E_GEOD_36980_Frontal_Cortex)[1:10]
head(E_GEOD_48350_Superior_Frontal_Gyrus)[1:10]
head(E_GEOD_5281_Superior_Frontal_Gyrus)[1:10]
head(inhouse_data_Frontal_Cortex)[1:10]
head(AMP_MSBB_U133A_Inferior_Frontal_Gyrus)[1:10]
head(AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex)[1:10]
head(AMP_MSBB_U133B_Inferior_Frontal_Gyrus)[1:10]
head(AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex)[1:10]
head(AMP_MSBB_U133A_Frontal_Pole)[1:10]
head(AMP_MSBB_U133B_Frontal_Pole)[1:10]
head(AMP_MSBB_U133A_Precentral_Gyrus)[1:10]
head(AMP_MSBB_U133B_Precentral_Gyrus)[1:10]
head(AMP_MSBB_U133A_Prefrontal_Cortex)[1:10]
head(AMP_MSBB_U133B_Prefrontal_Cortex)[1:10]

#Parietal_Lobe
head(AMP_MSBB_U133A_Superior_Parietal_Lobule)[1:10]
head(AMP_MSBB_U133B_Superior_Parietal_Lobule)[1:10]
head(E_GEOD_48350_Postcentral_Gyrus)[1:10]
head(E_GEOD_5281_Posterior_Singulate)[1:10]

#Temporal_Lobe
head(E_GEOD_28146_Hippocampus_CA1)[1:10]
head(E_GEOD_29378_Hippocampus_CA1)[1:10]
head(E_GEOD_29378_Hippocampus_CA3)[1:10]
head(E_GEOD_36980_Hippocampus)[1:10]
head(E_GEOD_36980_Temporal_Cortex)[1:10]
head(E_GEOD_48350_Entorhinal_Cortex)[1:10]
head(E_GEOD_48350_Hippocampus)[1:10]
head(E_GEOD_1297_Hippocampus)[1:10]
head(E_GEOD_5281_Entorhinal_Cortex)[1:10]
head(E_GEOD_5281_Hippocampus_CA1)[1:10]
head(E_GEOD_5281_Medial_Temporal_Gyrus)[1:10]
head(inhouse_data_Entorhinal_Cortex)[1:10]
head(inhouse_data_Temporal_Cortex)[1:10]
head(AMP_MAYOeGWAS_Temporal_Cortex)[1:10]
head(AMP_MSBB_U133A_Temporal_Pole)[1:10]
head(AMP_MSBB_U133B_Temporal_Pole)[1:10]
head(AMP_MSBB_U133A_Hippocampus)[1:10]
head(AMP_MSBB_U133A_Inferior_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133A_Middle_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133A_Superior_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133B_Parahippocampal_Gyrus)[1:10]
head(AMP_MSBB_U133B_Hippocampus)[1:10]
head(AMP_MSBB_U133B_Inferior_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133B_Middle_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133B_Superior_Temporal_Gyrus)[1:10]
head(AMP_MSBB_U133A_Parahippocampal_Gyrus)[1:10]

#Occiptal_Lobe
head(E_GEOD_5281_Primary_Visual_Cortex)[1:10]

#Cerebellum
head(inhouse_data_Cerebellum)[1:10]
head(AMP_MAYOeGWAS_Cerebellum)[1:10]

##### METAQC - CREATE FORMAT AND WRITE #####

# convert data to "MetaOmics" format - keep only diagnosis- create function for conversion of data

convert_data_to_MetaOmics_format<-function(dataset){
  # remove unwanted pheno - keep only diagnosis + expression (merge by rownames)
  dataset<-merge(dataset[5], dataset[10:dim(dataset)[2]], by="row.names")
  rownames(dataset)<-dataset$Row.names
  dataset$Row.names<-NULL
  # recode diagnosis column - case=1, control=0
  dataset[dataset$Diagnosis=="Control",1]<-"0"
  dataset[dataset$Diagnosis=="AD",1]<-"1"
  # re-label Diagnosis to Label
  colnames(dataset)[1]<-"label"
  # chane dataframe to numeric and transpose
  dataset<-as.data.frame(t(data.matrix(dataset)))
  print(head(dataset)[1:5])
  return(dataset)
}

# apply function to all datasets - leave out occipital lobe

#Frontal Lobe
E_GEOD_36980_Frontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_36980_Frontal_Cortex)
E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_48350_Superior_Frontal_Gyrus)
E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_5281_Superior_Frontal_Gyrus)
inhouse_data_Frontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(inhouse_data_Frontal_Cortex)
AMP_MSBB_U133A_Inferior_Frontal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Inferior_Frontal_Gyrus)
AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex)
AMP_MSBB_U133B_Inferior_Frontal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Inferior_Frontal_Gyrus)
AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex)
AMP_MSBB_U133A_Frontal_Pole_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Frontal_Pole)
AMP_MSBB_U133B_Frontal_Pole_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Frontal_Pole)
AMP_MSBB_U133A_Precentral_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Precentral_Gyrus)
AMP_MSBB_U133B_Precentral_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Precentral_Gyrus)
AMP_MSBB_U133A_Prefrontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Prefrontal_Cortex)
AMP_MSBB_U133B_Prefrontal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Prefrontal_Cortex)

#Parietal_Lobe_MetaOmics_format<-convert_data_to_MetaOmics_format(#Parietal_Lobe)
AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Superior_Parietal_Lobule)
AMP_MSBB_U133B_Superior_Parietal_Lobule_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Superior_Parietal_Lobule)
E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_48350_Postcentral_Gyrus)
E_GEOD_5281_Posterior_Singulate_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_5281_Posterior_Singulate)

#Temporal_Lobe_MetaOmics_format<-convert_data_to_MetaOmics_format(#Temporal_Lobe)
E_GEOD_28146_Hippocampus_CA1_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_28146_Hippocampus_CA1)
E_GEOD_29378_Hippocampus_CA1_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_29378_Hippocampus_CA1)
E_GEOD_29378_Hippocampus_CA3_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_29378_Hippocampus_CA3)
E_GEOD_36980_Hippocampus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_36980_Hippocampus)
E_GEOD_36980_Temporal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_36980_Temporal_Cortex)
E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_48350_Entorhinal_Cortex)
E_GEOD_48350_Hippocampus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_48350_Hippocampus)
E_GEOD_1297_Hippocampus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_1297_Hippocampus)
E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_5281_Entorhinal_Cortex)
E_GEOD_5281_Hippocampus_CA1_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_5281_Hippocampus_CA1)
E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(E_GEOD_5281_Medial_Temporal_Gyrus)
inhouse_data_Entorhinal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(inhouse_data_Entorhinal_Cortex)
inhouse_data_Temporal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(inhouse_data_Temporal_Cortex)
AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MAYOeGWAS_Temporal_Cortex)
AMP_MSBB_U133A_Temporal_Pole_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Temporal_Pole)
AMP_MSBB_U133B_Temporal_Pole_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Temporal_Pole)
AMP_MSBB_U133A_Hippocampus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Hippocampus)
AMP_MSBB_U133A_Inferior_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Inferior_Temporal_Gyrus)
AMP_MSBB_U133A_Middle_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Middle_Temporal_Gyrus)
AMP_MSBB_U133A_Superior_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Superior_Temporal_Gyrus)
AMP_MSBB_U133B_Parahippocampal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Parahippocampal_Gyrus)
AMP_MSBB_U133B_Hippocampus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Hippocampus)
AMP_MSBB_U133B_Inferior_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Inferior_Temporal_Gyrus)
AMP_MSBB_U133B_Middle_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Middle_Temporal_Gyrus)
AMP_MSBB_U133B_Superior_Temporal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133B_Superior_Temporal_Gyrus)
AMP_MSBB_U133A_Parahippocampal_Gyrus_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MSBB_U133A_Parahippocampal_Gyrus)

#Cerebellum_MetaOmics_format<-convert_data_to_MetaOmics_format(#Cerebellum)
inhouse_data_Cerebellum_MetaOmics_format<-convert_data_to_MetaOmics_format(inhouse_data_Cerebellum)
AMP_MAYOeGWAS_Cerebellum_MetaOmics_format<-convert_data_to_MetaOmics_format(AMP_MAYOeGWAS_Cerebellum)

# write out dataframe - ***NEED To CONVERT TO THIS FORMAT IN R

setwd(Data_in_metaomics_format)

# write out data

#Frontal lobe
write.table(E_GEOD_36980_Frontal_Cortex_MetaOmics_format, "E_GEOD_36980_Frontal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format, "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format, "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(inhouse_data_Frontal_Cortex_MetaOmics_format, "inhouse_data_Frontal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Inferior_Frontal_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_MetaOmics_format, "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Inferior_Frontal_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_MetaOmics_format, "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Frontal_Pole_MetaOmics_format, "AMP_MSBB_U133A_Frontal_Pole_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Frontal_Pole_MetaOmics_format, "AMP_MSBB_U133B_Frontal_Pole_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Precentral_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Precentral_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Precentral_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Precentral_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Prefrontal_Cortex_MetaOmics_format, "AMP_MSBB_U133A_Prefrontal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Prefrontal_Cortex_MetaOmics_format, "AMP_MSBB_U133B_Prefrontal_Cortex_MetaOmics_format.txt", sep="\t")

#Parietal_Lobe
write.table(AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format, "AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Superior_Parietal_Lobule_MetaOmics_format, "AMP_MSBB_U133B_Superior_Parietal_Lobule_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format, "E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5281_Posterior_Singulate_MetaOmics_format, "E_GEOD_5281_Posterior_Singulate_MetaOmics_format.txt", sep="\t")

#Temporal_Lobe
write.table(E_GEOD_28146_Hippocampus_CA1_MetaOmics_format, "E_GEOD_28146_Hippocampus_CA1_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_29378_Hippocampus_CA1_MetaOmics_format, "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_29378_Hippocampus_CA3_MetaOmics_format, "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_36980_Hippocampus_MetaOmics_format, "E_GEOD_36980_Hippocampus_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_36980_Temporal_Cortex_MetaOmics_format, "E_GEOD_36980_Temporal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format, "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_48350_Hippocampus_MetaOmics_format, "E_GEOD_48350_Hippocampus_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_1297_Hippocampus_MetaOmics_format, "E_GEOD_1297_Hippocampus_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format, "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5281_Hippocampus_CA1_MetaOmics_format, "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format, "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(inhouse_data_Entorhinal_Cortex_MetaOmics_format, "inhouse_data_Entorhinal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(inhouse_data_Temporal_Cortex_MetaOmics_format, "inhouse_data_Temporal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format, "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Temporal_Pole_MetaOmics_format, "AMP_MSBB_U133A_Temporal_Pole_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Temporal_Pole_MetaOmics_format, "AMP_MSBB_U133B_Temporal_Pole_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Hippocampus_MetaOmics_format, "AMP_MSBB_U133A_Hippocampus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Middle_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Middle_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Superior_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Superior_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Parahippocampal_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Parahippocampal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Hippocampus_MetaOmics_format, "AMP_MSBB_U133B_Hippocampus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Inferior_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Middle_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Middle_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133B_Superior_Temporal_Gyrus_MetaOmics_format, "AMP_MSBB_U133B_Superior_Temporal_Gyrus_MetaOmics_format.txt", sep="\t")
write.table(AMP_MSBB_U133A_Parahippocampal_Gyrus_MetaOmics_format, "AMP_MSBB_U133A_Parahippocampal_Gyrus_MetaOmics_format.txt", sep="\t")

#Cerebellum
write.table(inhouse_data_Cerebellum_MetaOmics_format, "inhouse_data_Cerebellum_MetaOmics_format.txt", sep="\t")
write.table(AMP_MAYOeGWAS_Cerebellum_MetaOmics_format, "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format.txt", sep="\t")

##### METAQC - READ IN CREATED FORMAT #####

setwd(Data_in_metaomics_format)

# read data by brain region

# Frontal Lobe

Frontal_Lobe<-MetaDE.Read(c("E_GEOD_36980_Frontal_Cortex_MetaOmics_format", 
                         "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format", 
                         "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format", 
                         "inhouse_data_Frontal_Cortex_MetaOmics_format", 
                         "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_MetaOmics_format", 
                         "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_MetaOmics_format", 
                         "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_MetaOmics_format", 
                         "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_MetaOmics_format", 
                         "AMP_MSBB_U133A_Frontal_Pole_MetaOmics_format", 
                         "AMP_MSBB_U133B_Frontal_Pole_MetaOmics_format", 
                         "AMP_MSBB_U133A_Precentral_Gyrus_MetaOmics_format", 
                         "AMP_MSBB_U133B_Precentral_Gyrus_MetaOmics_format", 
                         "AMP_MSBB_U133A_Prefrontal_Cortex_MetaOmics_format", 
                         "AMP_MSBB_U133B_Prefrontal_Cortex_MetaOmics_format"), # 14 datasets
                       via="txt", # txt format
                       skip=rep(1,14), # skip 1 row - 14 datasets
                       log=T, # data is log transformed
                       matched=T) #using entrez gene ID instea of gene symbols

Parietal_Lobe<-MetaDE.Read(c("AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format", 
                            "AMP_MSBB_U133B_Superior_Parietal_Lobule_MetaOmics_format", 
                            "E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
                            "E_GEOD_5281_Posterior_Singulate_MetaOmics_format"), # 4 datasets
                          via="txt", # txt format
                          skip=rep(1,4), # skip 1 row - 4 datasets
                          log=T, # data is log transformed
                          matched=T) #using entrez gene ID instea of gene symbols

Temporal_Lobe<-MetaDE.Read(c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                             "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                             "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                             "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                             "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                             "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                             "inhouse_data_Temporal_Cortex_MetaOmics_format", 
                             "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format", 
                             "AMP_MSBB_U133A_Temporal_Pole_MetaOmics_format", 
                             "AMP_MSBB_U133B_Temporal_Pole_MetaOmics_format", 
                             "AMP_MSBB_U133A_Hippocampus_MetaOmics_format", 
                             "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133A_Middle_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133A_Superior_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133B_Parahippocampal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133B_Hippocampus_MetaOmics_format", 
                             "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133B_Middle_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133B_Superior_Temporal_Gyrus_MetaOmics_format", 
                             "AMP_MSBB_U133A_Parahippocampal_Gyrus_MetaOmics_format"), # 26 datasets
                           via="txt", # txt format
                           skip=rep(1,26), # skip 1 row - 26 datasets
                           log=T, # data is log transformed
                           matched=T) #using entrez gene ID instea of gene symbols

Cerebellum<-MetaDE.Read(c("inhouse_data_Cerebellum_MetaOmics_format", 
                          "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"), # 2 datasets
                           via="txt", # txt format
                           skip=rep(1,2), # skip 1 row - 2 datasets
                           log=T, # data is log transformed
                        matched=T) #using entrez gene ID instea of gene symbols

setwd(work_dir)

##### EXTRACT COMMON PROBES - 20% #####

# required step for QC plot- will be removed during DE analysis

# no common probes - full data merge - minimum 20% common required due to knn.impute function

Frontal_Lobe_merged<-MetaDE.merge(Frontal_Lobe, MVperc=0.8)
Parietal_Lobe_merged<-MetaDE.merge(Parietal_Lobe, MVperc=0.8)
Temporal_Lobe_merged<-MetaDE.merge(Temporal_Lobe, MVperc=0.8)
Cerebellum_merged<-MetaDE.merge(Cerebellum, MVperc=0.8)

dim(Frontal_Lobe_merged[[1]][[1]])
dim(Parietal_Lobe_merged[[1]][[1]])
dim(Temporal_Lobe_merged[[1]][[1]])
dim(Cerebellum_merged[[1]][[1]])

##### FILTER GENES #####

#filter out 30% un-expressed genes and then 30% non-informative genes across datasets

Frontal_Lobe_merged_Filtered<-MetaDE.filter(Frontal_Lobe_merged,c(0.3,0.3)) 
Parietal_Lobe_merged_Filtered<-MetaDE.filter(Parietal_Lobe_merged,c(0.3,0.3)) 
Temporal_Lobe_merged_Filtered<-MetaDE.filter(Temporal_Lobe_merged,c(0.3,0.3)) 
Cerebellum_merged_Filtered<-MetaDE.filter(Cerebellum_merged,c(0.3,0.3)) 

dim(Frontal_Lobe_merged_Filtered[[1]][[1]])
dim(Parietal_Lobe_merged_Filtered[[1]][[1]])
dim(Temporal_Lobe_merged_Filtered[[1]][[1]])
dim(Cerebellum_merged_Filtered[[1]][[1]])

##### META QC - CREATE QC OBJECT #####

# meta QC

Frontal_Lobe_QC<-list()
Parietal_Lobe_QC<-list()
Temporal_Lobe_QC<-list()
Cerebellum_QC<-list()

# moves diagnosis label to column name - Frontal_Lobe_merged_Filtered - impute missing - 
for(i in 1:14){ 
  colnames(Frontal_Lobe_merged_Filtered[[i]][[1]])<-Frontal_Lobe_merged_Filtered[[i]][[2]]
  Frontal_Lobe_QC[[i]]<-impute.knn(Frontal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Frontal_Lobe_QC)<-names(Frontal_Lobe_merged_Filtered)

head(Frontal_Lobe_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)[,1:5]
dim(Frontal_Lobe_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)

# moves diagnosis label to column name - Parietal_Lobe_merged_Filtered - impute missing - 
for(i in 1:4){ 
  colnames(Parietal_Lobe_merged_Filtered[[i]][[1]])<-Parietal_Lobe_merged_Filtered[[i]][[2]]
  Parietal_Lobe_QC[[i]]<-impute.knn(Parietal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Parietal_Lobe_QC)<-names(Parietal_Lobe_merged_Filtered)

head(Parietal_Lobe_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)[,1:5]
dim(Parietal_Lobe_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)

# moves diagnosis label to column name - Temporal_Lobe_merged_Filtered - impute missing - 
for(i in 1:26){ 
  colnames(Temporal_Lobe_merged_Filtered[[i]][[1]])<-Temporal_Lobe_merged_Filtered[[i]][[2]]
  Temporal_Lobe_QC[[i]]<-impute.knn(Temporal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Temporal_Lobe_QC)<-names(Temporal_Lobe_merged_Filtered)

head(Temporal_Lobe_QC$E_GEOD_28146_Hippocampus_CA1_MetaOmics_format)[,1:5]
dim(Temporal_Lobe_QC$E_GEOD_28146_Hippocampus_CA1_MetaOmics_format)

# moves diagnosis label to column name - Cerebellum_merged_Filtered - impute missing - 
for(i in 1:2){ 
  colnames(Cerebellum_merged_Filtered[[i]][[1]])<-Cerebellum_merged_Filtered[[i]][[2]]
  Cerebellum_QC[[i]]<-impute.knn(Cerebellum_merged_Filtered[[i]][[1]])$data
}

names(Cerebellum_QC)<-names(Cerebellum_merged_Filtered)

head(Cerebellum_QC$inhouse_data_Cerebellum_MetaOmics_format)[,1:5]
dim(Cerebellum_QC$inhouse_data_Cerebellum_MetaOmics_format)

##### RUN QC CHECK #####

# create function to run quantitative quality control measures

run_meta_qc<-function(dataset) {
  MetaQC(dataset, # data
  "c2.cp.biocarta.v5.1.entrez.gmt", # latest version downloaded 02/09/2016 
  isParallel = T, # parallel processing
  nCores = 8, # n.o cores
  useCache = TRUE, #imported file to be used again?
  filterGenes = F, #gene filtered by missing - already done earlier
  #maxNApctAllowed=.3, # ratio of missing genes default .3
  #cutRatioByMean=.4, # filtering of genes if unexpressed
  #cutRatioByVar=.4, # filtering genes least sample wise expression variance
  minNumGenes=5, # min genes in pathway
  verbose = FALSE, 
  resp.type = "Twoclass") # c("Twoclass", "Multiclass", "Survival"))
}

Frontal_Lobe_QC_check<-run_meta_qc(Frontal_Lobe_QC)
Parietal_Lobe_QC_check<-run_meta_qc(Parietal_Lobe_QC)
Temporal_Lobe_QC_check<-run_meta_qc(Temporal_Lobe_QC)
Cerebellum_QC_check<-run_meta_qc(Cerebellum_QC)

Frontal_Lobe_QC_check
Parietal_Lobe_QC_check
Temporal_Lobe_QC_check
Cerebellum_QC_check

# create function to run 2nd part of QC

# runQC- **script freezes when 1st downloading fileForCQCp. downloaded manually and inserted into working folder http://software.broadinstitute.org/gsea/downloads.jsp

run_meta_qc2<-function(dataset){
  runQC(dataset, 
        nPath=NULL, # 
        B=1e4, # The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
        pvalCut=.05, #P-value threshold used for AQC calculation.
        #pvalAdjust=T, #whether to apply p-value adjustment due to multiple testing (B-H procedure is used).
        fileForCQCp="c2.all.v5.1.entrez.gmt") # latest version of c2.all.v5.1.symbols.gmt downloaded 02/09/2016
}
 
run_meta_qc2(Frontal_Lobe_QC_check)
run_meta_qc2(Parietal_Lobe_QC_check)
run_meta_qc2(Temporal_Lobe_QC_check)
#run_meta_qc2(Cerebellum_QC_check) # - "Error in { : task 1 failed - "not enough 'y' observations"

Frontal_Lobe_QC_check
Parietal_Lobe_QC_check
Temporal_Lobe_QC_check
#Cerebellum_QC_check

#plot PCA

plot(Frontal_Lobe_QC_check)
plot(Parietal_Lobe_QC_check)
plot(Temporal_Lobe_QC_check)
# plot(Cerebellum_QC_check)

##### PLOT QC PLOTS TO PDF #####

setwd(MetaQC_plots_dir)

pdf("Frontal_lobe_MetaQC.pdf")
plot(Frontal_Lobe_QC_check)
mtext("Frontal Lobe Datasets", cex=2, line = 1.5)
dev.off()

pdf("Parietal_lobe_MetaQC.pdf")
plot(Parietal_Lobe_QC_check)
mtext("Parietal Lobe Datasets", cex=2, line = 1.5)
dev.off()

pdf("Temporal_lobe_MetaQC.pdf")
plot(Temporal_Lobe_QC_check)
mtext("Temporal Lobe Datasets", cex=2, line = 1.5)
dev.off()

##### CREATE FUNCTION TO READ DATA IN METAOMICS FORMAT AND RUN QC #####

#some datasets may affect result. created function to read data in MetaOmics format and run qc

run_MetaQC<-function(datasets, MVperc, gene_filter){ 
  #datasets = list of dataset names, 
  #MVperv = merge by common probes, 0=100% common, 0.8=20% common
  #gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  #number of datasets
  number_of_datasets<-length(datasets) 
  #setwd to data
  setwd(Data_in_metaomics_format)
  #read in raw data
  raw_data<-MetaDE.Read(datasets,
                        via="txt", # txt format
                        skip=rep(1,number_of_datasets), # skip 1 row for all datasets
                        log=T, # data is log transformed
                        matched=T) # using entrez gene ID
  #setwd to work dir
  setwd(work_dir)
  # merge by common probes, 0=100% common, 0.8=20% common
  raw_data_merged<-MetaDE.merge(raw_data, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_merged[[1]][[1]])[1]), sep=" ")
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_merged_Filtered<-MetaDE.filter(raw_data_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_merged_Filtered[[1]][[1]])[1]), sep=" ")
  # create QC object 
  Data_QC<-list()
  # moves diagnosis label to column name - Frontal_Lobe_merged_Filtered 
  for(i in 1:number_of_datasets){ 
    colnames(raw_data_merged_Filtered[[i]][[1]])<-raw_data_merged_Filtered[[i]][[2]]
    Data_QC[[i]]<-impute.knn(raw_data_merged_Filtered[[i]][[1]])$data
  }
  names(Data_QC)<-names(raw_data_merged_Filtered)
  # run quantitative quality control measures
  Data_QC_measures<-MetaQC(Data_QC, # data
                           "c2.cp.biocarta.v5.1.entrez.gmt", # latest version downloaded 02/09/2016 
                           isParallel = T, # parallel processing
                           nCores = 8, # n.o cores
                           useCache = TRUE, #imported file to be used again?
                           filterGenes = F, #gene filtered by missing - already done earlier
                           minNumGenes=5, # min genes in pathway
                           verbose = FALSE, 
                           resp.type = "Twoclass") # c("Twoclass", "Multiclass", "Survival"))
  #run 2nd part of QC
  runQC(Data_QC_measures, 
        nPath=NULL, # 
        B=1e4, # The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
        pvalCut=.05, #P-value threshold used for AQC calculation.
        fileForCQCp="c2.all.v5.1.entrez.gmt") # latest version of c2.all.v5.1.symbols.gmt downloaded 02/09/2016
  plot(Data_QC_measures)
  return(Data_QC_measures)
}

# apply function to lobe regions without AMP MSBB datasets

Frontal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_36980_Frontal_Cortex_MetaOmics_format", 
                                              "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format", 
                                              "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format", 
                                              "inhouse_data_Frontal_Cortex_MetaOmics_format"),
                                   MVperc=0.1,
                                   gene_filter=c(0.3,0.3))

# Parietal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
#                                                "E_GEOD_5281_Posterior_Singulate_MetaOmics_format"),
#                                     MVperc=0.2,
#                                     gene_filter=c(0.3,0.3)) # error not enought 'y' observations

                                    
# Temporal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
#                                                "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
#                                                "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
#                                                "E_GEOD_36980_Hippocampus_MetaOmics_format", 
#                                                "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
#                                                "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
#                                                "E_GEOD_48350_Hippocampus_MetaOmics_format", 
#                                                "E_GEOD_1297_Hippocampus_MetaOmics_format", 
#                                                "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
#                                                "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
#                                                "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
#                                                "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
#                                                "inhouse_data_Temporal_Cortex_MetaOmics_format",
#                                                "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
#                                     MVperc=0.5,
#                                     gene_filter=c(0.3,0.3))

# Cerebellum_no_AMP_QC<-run_MetaQC(datasets=c("inhouse_data_Cerebellum_MetaOmics_format", 
#                                   "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"),
#                        MVperc=0.2,
#                        gene_filter=c(0.3,0.3)) # error not enought 'y' observations

# plot to pdf

setwd(MetaQC_plots_dir)

pdf("Frontal_lobe_MetaQC_NO_AMP.pdf")
plot(Frontal_lobe_no_AMP_QC)
mtext("Frontal Lobe Datasets (NO AMP)", cex=2, line = 1.5)
dev.off()

pdf("Temporal_lobe_MetaQC_NO_AMP.pdf")
plot(Temporal_lobe_no_AMP_QC)
mtext("Temporal Lobe Datasets (NO AMP)", cex=2, line = 1.5)
dev.off()

setwd(work_dir)

##### TEST - TEMPORAL LOBE MERGE USING DIFFERENT MISSING RATE PCA PLOT #####

test2<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                             "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                             "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                             "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                             "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                             "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                             "inhouse_data_Temporal_Cortex_MetaOmics_format",
                             "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                  MVperc=0.7,
                  gene_filter=c(0.3,0.3))

test3<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                             "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                             "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                             "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                             "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                             "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                             "inhouse_data_Temporal_Cortex_MetaOmics_format",
                             "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                  MVperc=0.6,
                  gene_filter=c(0.3,0.3))

test3.5<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                               "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                               "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                               "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                               "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                               "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                               "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                               "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                               "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                               "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                               "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                               "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                               "inhouse_data_Temporal_Cortex_MetaOmics_format",
                               "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                    MVperc=0.5,
                    gene_filter=c(0.3,0.3))


test4<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                             "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                             "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                             "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                             "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                             "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                             "inhouse_data_Temporal_Cortex_MetaOmics_format",
                             "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                  MVperc=0.4,
                  gene_filter=c(0.3,0.3))


test5<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                             "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                             "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                             "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                             "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                             "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                             "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                             "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                             "inhouse_data_Temporal_Cortex_MetaOmics_format",
                             "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                  MVperc=0.2,
                  gene_filter=c(0.3,0.3))


setwd(MetaQC_plots_dir)

pdf("Temporal_lobe_no_AMP_PCA_plot_test.pdf")
plot(test2)
mtext("Temporal Lobe Datasets (NO AMP) 30% missing", cex=1.5, line = 1)

plot(test3)
mtext("Temporal Lobe Datasets (NO AMP) 40% missing", cex=1.5, line = 1)

plot(test3.5)
mtext("Temporal Lobe Datasets (NO AMP) 50% missing", cex=1.5, line = 1)

plot(test4)
mtext("Temporal Lobe Datasets (NO AMP) 60% missing", cex=1.5, line = 1)

plot(test5)
mtext("Temporal Lobe Datasets (NO AMP) 80% missing", cex=1.5, line = 1)
dev.off()

##### RE-CREATE FILTERED OBJECT AFTER STUDY REMOVAL #####

#create function to create QC'd object - removing unwanted datasets
create_QC_object<-function(datasets, MVperc, gene_filter){ 
  #datasets = list of dataset names, 
  #MVperv = merge by common probes, 0=100% common, 0.8=20% common
  #gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  #number of datasets
  number_of_datasets<-length(datasets) 
  
  #setwd to data
  setwd(Data_in_metaomics_format)
  
  #read in raw data
  raw_data<-MetaDE.Read(datasets,
                        via="txt", # txt format
                        skip=rep(1,number_of_datasets), # skip 1 row for all datasets
                        log=T, # data is log transformed
                        matched=F) # matched format false - i.e need to match entrez id to gene symbol
  
  #setwd to work dir
  setwd(work_dir)
  
  # merge by common probes, 0=100% common, 0.8=20% common
  raw_data_merged<-MetaDE.merge(raw_data, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_merged[[1]][[1]])[1]), sep=" ")
  
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_merged_Filtered<-MetaDE.filter(raw_data_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_merged_Filtered[[1]][[1]])[1]), sep=" ")
  return(raw_data_merged_Filtered)
}

#apply function to each brain region - no merge by common probe and no filtering of probes

Frontal_lobe_no_AMP_QC_filtered<-create_QC_object(datasets=c("E_GEOD_36980_Frontal_Cortex_MetaOmics_format", 
                                                             "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format", 
                                                             "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format", 
                                                             "inhouse_data_Frontal_Cortex_MetaOmics_format"),
                                                  MVperc=0.8,
                                                  gene_filter=c(0,0))

Parietal_lobe_no_AMP_QC_filtered<-create_QC_object(datasets=c("E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
                                                              "E_GEOD_5281_Posterior_Singulate_MetaOmics_format"),
                                                   MVperc=0.8,
                                                   gene_filter=c(0,0))


Temporal_lobe_no_AMP_QC_filtered<-create_QC_object(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                                                              "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                                                              "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                                                              "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                                                              "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                                                              "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                                                              "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                                                              "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                                                              "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                                                              "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                                                              "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                                                              "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                                                              "inhouse_data_Temporal_Cortex_MetaOmics_format",
                                                              "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format"),
                                                   MVperc=0.8,
                                                   gene_filter=c(0,0))

Cerebellum_no_AMP_QC_filtered<-create_QC_object(datasets=c("inhouse_data_Cerebellum_MetaOmics_format", 
                                                           "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"),
                                                MVperc=0.8,
                                                gene_filter=c(0,0))

##### META DE ANALYSIS "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC" #####

# Frontal Lobe
Frontal_lobe_DE<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
                                ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                nperm=100,
                                miss.tol=0,
                                ind.tail="high")

tail(Frontal_lobe_DE$meta.analysis$stat)
tail(Frontal_lobe_DE$meta.analysis$FDR)
tail(Frontal_lobe_DE$meta.analysis$AW.weight)
tail(Frontal_lobe_DE$meta.analysis$pval)

#Parietal Lobe
Parietal_lobe_DE<-MetaDE.rawdata(Parietal_lobe_no_AMP_QC_filtered,
                                 ind.method = rep("modt", length(names(Parietal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                 meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                 nperm=100,
                                 miss.tol=0,
                                 ind.tail="high")

tail(Parietal_lobe_DE$meta.analysis$stat)
tail(Parietal_lobe_DE$meta.analysis$FDR)
tail(Parietal_lobe_DE$meta.analysis$AW.weight)
tail(Parietal_lobe_DE$meta.analysis$pval)

#Cerebellum
Cerebellum_DE<-MetaDE.rawdata(Cerebellum_no_AMP_QC_filtered,
                              ind.method = rep("modt", length(names(Cerebellum_no_AMP_QC_filtered))),# moderate t-test - same as limma
                              meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                              nperm=100,
                              miss.tol=0,
                              ind.tail="high")

tail(Cerebellum_DE$meta.analysis$stat)
tail(Cerebellum_DE$meta.analysis$FDR)
tail(Cerebellum_DE$meta.analysis$AW.weight)
tail(Cerebellum_DE$meta.analysis$pval)

#Temporal - requires 251.1G
# Temporal_lobe_DE<-MetaDE.rawdata(Temporal_lobe_no_AMP_QC_filtered,
#                                     ind.method = rep("modt", length(names(Temporal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                     meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
#                                     nperm=20,
#                                     miss.tol=0,
#                                     ind.tail="high")
# 


# temporal lobe fails to aquire 251 ram reqiured - save and run on cluster
setwd(work_dir)
save(Temporal_lobe_no_AMP_QC_filtered, file="Temporal_lobe_no_AMP_QC_filtered.Rdata")

#read in  temporal DE run on cluster

load("Meta_DE/Temporal_lobe_run_on_cluster/Temporal_lobe_DE.Rdata")
ls()
Temporal_lobe_DE<-Temporal_lobe_DE_up
rm(Temporal_lobe_DE_up)

tail(Temporal_lobe_DE$meta.analysis$stat)
tail(Temporal_lobe_DE$meta.analysis$FDR)
tail(Temporal_lobe_DE$meta.analysis$AW.weight)
tail(Temporal_lobe_DE$meta.analysis$pval)

##### META DE ANALYSIS ON COMBINED BRAIN REGIONS #####

# added 5/12/2016
#combine all brain regions into one meta-analysis

combined_brain_datasets_no_AMP_QC_filtered<-create_QC_object(datasets=c("E_GEOD_36980_Frontal_Cortex_MetaOmics_format", 
                                                                        "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format", 
                                                                        "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format", 
                                                                        "inhouse_data_Frontal_Cortex_MetaOmics_format",
                                                                        "E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
                                                                        "E_GEOD_5281_Posterior_Singulate_MetaOmics_format",
                                                                        "E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
                                                                        "E_GEOD_29378_Hippocampus_CA1_MetaOmics_format", 
                                                                        "E_GEOD_29378_Hippocampus_CA3_MetaOmics_format", 
                                                                        "E_GEOD_36980_Hippocampus_MetaOmics_format", 
                                                                        "E_GEOD_36980_Temporal_Cortex_MetaOmics_format", 
                                                                        "E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format", 
                                                                        "E_GEOD_48350_Hippocampus_MetaOmics_format", 
                                                                        "E_GEOD_1297_Hippocampus_MetaOmics_format", 
                                                                        "E_GEOD_5281_Entorhinal_Cortex_MetaOmics_format", 
                                                                        "E_GEOD_5281_Hippocampus_CA1_MetaOmics_format", 
                                                                        "E_GEOD_5281_Medial_Temporal_Gyrus_MetaOmics_format", 
                                                                        "inhouse_data_Entorhinal_Cortex_MetaOmics_format", 
                                                                        "inhouse_data_Temporal_Cortex_MetaOmics_format",
                                                                        "AMP_MAYOeGWAS_Temporal_Cortex_MetaOmics_format",
                                                                        "inhouse_data_Cerebellum_MetaOmics_format", 
                                                                        "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"),
                                                             MVperc=0.8,
                                                             gene_filter=c(0,0))
setwd(work_dir)
save(combined_brain_datasets_no_AMP_QC_filtered, file="combined_brain_datasets_no_AMP_QC_filtered.Rdata")

#run on rosalind - need ~200gb

# combined_brain_region_DE<-MetaDE.rawdata(combined_brain_datasets_no_AMP_QC_filtered,
#                               ind.method = rep("modt", length(names(combined_brain_datasets_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                               meta.method="AW.OC",
#                               nperm=100,
#                               miss.tol=0,
#                               ind.tail="high")

tail(combined_brain_region_DE$meta.analysis$stat)
tail(combined_brain_region_DE$meta.analysis$FDR)
tail(combined_brain_region_DE$meta.analysis$AW.weight)
tail(combined_brain_region_DE$meta.analysis$pval)

# ##### AW.OC DE METHOD TEST #####
# 
# # added 6/9/16
# # AW weights and limma logFC result discrepency -. testing different AW methods with different "ind.tail" argument
# 
# Frontal_lobe_DE_AW.OC_up<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
#                                 ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                 meta.method=c("AW.OC"),
#                                 nperm=100,
#                                 miss.tol=0,
#                                 ind.tail="high")
# 
# Frontal_lobe_DE_AW.OC_down<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
#                                          ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                          meta.method=c("AW.OC"),
#                                          nperm=100,
#                                          miss.tol=0,
#                                          ind.tail="low")
# 
# Frontal_lobe_DE_AW<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
#                                            ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                            meta.method=c("AW"),
#                                            nperm=100,
#                                            miss.tol=0,
#                                    ind.tail="abs")
# 
# 
# 
# tail(Frontal_lobe_DE_AW.OC_up$meta.analysis$stat)
# tail(Frontal_lobe_DE_AW.OC_up$meta.analysis$FDR)
# tail(Frontal_lobe_DE_AW.OC_up$meta.analysis$AW.weight)
# tail(Frontal_lobe_DE_AW.OC_up$meta.analysis$pval)
# count.DEnumber(Frontal_lobe_DE_AW.OC_up, p.cut=0.05, q.cut=0.05)
# 
# tail(Frontal_lobe_DE_AW.OC_down$meta.analysis$stat)
# tail(Frontal_lobe_DE_AW.OC_down$meta.analysis$FDR)
# tail(Frontal_lobe_DE_AW.OC_down$meta.analysis$AW.weight)
# tail(Frontal_lobe_DE_AW.OC_down$meta.analysis$pval)
# count.DEnumber(Frontal_lobe_DE_AW.OC_down, p.cut=0.05, q.cut=0.05)
# 
# tail(Frontal_lobe_DE_AW.OC$meta.analysis$stat)
# tail(Frontal_lobe_DE_AW.OC$meta.analysis$FDR)
# tail(Frontal_lobe_DE_AW.OC$meta.analysis$AW.weight)
# tail(Frontal_lobe_DE_AW.OC$meta.analysis$pval)
# count.DEnumber(Frontal_lobe_DE_AW.OC, p.cut=0.05, q.cut=0.05)

##### HEATMAP "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC" #####

# meta DE analysis done using "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC". plot heat map for each method
setwd(DE_Heamaps_dir)

pdf("Frontal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Frontal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe AW.OC DEG ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe maxP.OC DEG ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe minP.OC DEG ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe Fisher.OC DEG ", line = 1.5)
dev.off()

#parietal lobe
pdf("Parietal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Parietal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe AW.OC DEG ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe maxP.OC DEG ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe minP.OC DEG ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe Fisher.OC DEG ", line = 1.5)
dev.off()

# cerebellum
pdf("Cerebellum_DE_heatmaps.pdf")
heatmap.sig.genes(Cerebellum_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum AW.OC DEG ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum maxP.OC DEG ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum minP.OC DEG ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum Fisher.OC DEG ", line = 1.5)
dev.off()

#Temporal Lobe
pdf("Temporal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Temporal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe AW.OC DEG ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe maxP.OC DEG ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe minP.OC DEG ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe Fisher.OC DEG ", line = 1.5)
dev.off()

##### META DE ANALYSIS FEM REM #####

# Frontal Lobe - FEM
Frontal_lobe_DE_FEM<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
                                    ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                    meta.method=c("FEM"),
                                    nperm=1000,
                                    miss.tol=0,
                                    paired=rep(F,length(names(Frontal_lobe_no_AMP_QC_filtered))),
                                    ind.tail="abs")

#Parietal Lobe - FEM
Parietal_lobe_DE_FEM<-MetaDE.rawdata(Parietal_lobe_no_AMP_QC_filtered,
                                     ind.method = rep("modt", length(names(Parietal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                     meta.method=c("FEM"),
                                     nperm=1000,
                                     miss.tol=0,
                                     paired=rep(F,length(names(Parietal_lobe_no_AMP_QC_filtered))),
                                     ind.tail="abs")
# Temporal Lobe - FEM
Temporal_lobe_DE_FEM<-MetaDE.rawdata(Temporal_lobe_no_AMP_QC_filtered,
                                     ind.method = rep("modt", length(names(Temporal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                     meta.method=c("FEM"),
                                     nperm=1000,
                                     miss.tol=0,
                                     paired=rep(F,length(names(Temporal_lobe_no_AMP_QC_filtered))),
                                     ind.tail="abs")
# Cerebellum -FEM
Cerebellum_DE_FEM<-MetaDE.rawdata(Cerebellum_no_AMP_QC_filtered,
                                  ind.method = rep("modt", length(names(Cerebellum_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                  meta.method=c("FEM"),
                                  nperm=1000,
                                  miss.tol=0,
                                  paired=rep(F,length(names(Cerebellum_no_AMP_QC_filtered))),
                                  ind.tail="abs")


# Frontal Lobe - REM
Frontal_lobe_DE_REM<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
                                    ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                    meta.method=c("REM"),
                                    nperm=1000,
                                    miss.tol=0,
                                    paired=rep(F,length(names(Frontal_lobe_no_AMP_QC_filtered))),
                                    ind.tail="abs")

#Parietal Lobe - REM
Parietal_lobe_DE_REM<-MetaDE.rawdata(Parietal_lobe_no_AMP_QC_filtered,
                                     ind.method = rep("modt", length(names(Parietal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                     meta.method=c("REM"),
                                     nperm=1000,
                                     miss.tol=0,
                                     paired=rep(F,length(names(Parietal_lobe_no_AMP_QC_filtered))),
                                     ind.tail="abs")
# Temporal Lobe - REM
Temporal_lobe_DE_REM<-MetaDE.rawdata(Temporal_lobe_no_AMP_QC_filtered,
                                     ind.method = rep("modt", length(names(Temporal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                     meta.method=c("REM"),
                                     nperm=1000,
                                     miss.tol=0,
                                     paired=rep(F,length(names(Temporal_lobe_no_AMP_QC_filtered))),
                                     ind.tail="abs")
# Cerebellum -REM
Cerebellum_DE_REM<-MetaDE.rawdata(Cerebellum_no_AMP_QC_filtered,
                                  ind.method = rep("modt", length(names(Cerebellum_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                  meta.method=c("REM"),
                                  nperm=1000,
                                  miss.tol=0,
                                  paired=rep(F,length(names(Cerebellum_no_AMP_QC_filtered))),
                                  ind.tail="abs")

##### COUNT DE NUMBERS #####

count.DEnumber(Frontal_lobe_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Parietal_lobe_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Temporal_lobe_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Cerebellum_DE,
               p.cut=0.05,
               q.cut=0.05)

##### DE NUMBERS PER META_METHOD - PLOT #####

setwd(Meta_DE_dir)

pdf("DE_numbers_for_each_method.pdf")

draw.DEnumber(Frontal_lobe_DE,
              0.05,
              FDR=T)
mtext("Frontal Lobe DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Parietal_lobe_DE,
              0.05,
              FDR=T)
mtext("Parietal Lobe DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Cerebellum_DE,
              0.05,
              FDR=T)
mtext("Cerebellum DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Temporal_lobe_DE,
              0.05,
              FDR=T)
mtext("Temporal Lobe DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

dev.off()

##### EXTRACT DEG PER BRAIN REGION ####

# using  ind.tail="high" and ind.tail="low" gives same number of DEG using all meta-analysis methods - using only up-genes object

Frontal_lobe_DE_AW.OC<-(as.data.frame(subset(Frontal_lobe_DE$meta.analysis$FDR, Frontal_lobe_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Parietal_lobe_DE_AW.OC<-(as.data.frame(subset(Parietal_lobe_DE$meta.analysis$FDR, Parietal_lobe_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Cerebellum_DE_AW.OC<-(as.data.frame(subset(Cerebellum_DE$meta.analysis$FDR, Cerebellum_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Temporal_lobe_DE_AW.OC<-(as.data.frame(subset(Temporal_lobe_DE$meta.analysis$FDR, Temporal_lobe_DE$meta.analysis$FDR[,1]<=0.05)))[1]

dim(Temporal_lobe_DE_AW.OC)
head(Temporal_lobe_DE_AW.OC)

dim(Frontal_lobe_DE_AW.OC)
head(Frontal_lobe_DE_AW.OC)

dim(Parietal_lobe_DE_AW.OC)
head(Parietal_lobe_DE_AW.OC)

dim(Cerebellum_DE_AW.OC)
head(Cerebellum_DE_AW.OC)

# read in logC of each dataset

limma_results<-read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/3.DE_analysis/Differential_Expression/Differential_expression_on_full_expression_data_AD_vs_CONTROL.txt",
                          header=T,
                          as.is=T,
                          sep="\t")

head(limma_results)[1:5]

# extract fold change only
fold_change_columns<-grep("logFC", colnames(limma_results))
log_FC_results<-limma_results[,fold_change_columns]

head(log_FC_results)[1:5]
dim(log_FC_results)

# remove logFC name from colnames
colnames(log_FC_results)<-sub("_logFC", "", colnames(log_FC_results))

# create function - 
# 1.extract weights
# 2.subset limma results to match studies in meta analysis
# 3.according to AW weights - make any study where weight=0 to 0 in individual study stats
# 4.keep only values that have contributed to AW.OC meta-analysis.
# 5.make new column to count number of negative logFC values across studies per entrez gene id
# 6.make new column to count number of positive logFC values across studies per entrez gene id 
# 7.decide if genes are up/down regulated


calculate_AW_gene_regulation<-function(dataset, logFC_results){
  #extract sig meta P values
  sig_p<-(as.data.frame(subset(dataset$meta.analysis$FDR, dataset$meta.analysis$FDR[,1]<=0.05)))[1]
  # extract meta-DEG list AW weights - add colnames to match individual studies
  dataset_AW.OC_weight<-as.data.frame(dataset$meta.analysis$AW.weight[rownames(sig_p),])
  colnames(dataset_AW.OC_weight)<-colnames(dataset$ind.stat)
  #remove "_MetaOmics_format" from colnames
  colnames(dataset_AW.OC_weight)<-sub("_MetaOmics_format", "", colnames(dataset_AW.OC_weight))
  # remove "_data"
  colnames(dataset_AW.OC_weight)<-sub("_data", "", colnames(dataset_AW.OC_weight))
  # subset limma results to match dataset - columns
  log_FC_results_subset<-log_FC_results[,colnames(log_FC_results) %in% colnames(dataset_AW.OC_weight)]
  # subset limma rownames results to match dataset
  log_FC_results_subset<-log_FC_results_subset[rownames(log_FC_results_subset) %in% rownames(dataset_AW.OC_weight),]
  # sort limma results to match dataset order
  log_FC_results_subset<-log_FC_results_subset[order(rownames(log_FC_results_subset)),]
  dataset_AW.OC_weight<-dataset_AW.OC_weight[order(rownames(dataset_AW.OC_weight)),]
  # check number of studies match
  print(paste("Number of studies in limma results match meta analysis results: ", length(colnames(log_FC_results_subset))==length(colnames(dataset_AW.OC_weight)), sep=""))
  # check studies are ordered correctly across both dataframes
  print(paste("Studies correctly ordered by colnames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
  # check rownames are same
  print(paste("Studies correctly ordered by rownames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
  # print number of genes and studies
  print(paste("Number of Studies: ", dim(dataset_AW.OC_weight)[2], sep=""))
  print(paste("Number of Genes: ", dim(dataset_AW.OC_weight)[1], sep=""))
  # convert logFC to 0 if AW.OC weight==0
  for (x in 1:dim(dataset_AW.OC_weight)[1]){
    for (y in 1:dim(dataset_AW.OC_weight)[2])
      if (dataset_AW.OC_weight[x,y]==0){
        log_FC_results_subset[x,y]<-0
      }
  }
  # make column down, up, problem according to logFC signs.
  # count number of postive + negative values per row
  # create temp dataframe with rownames
  positive_temp<-log_FC_results_subset[1]
  negative_temp<-log_FC_results_subset[1]
  # rename columns
  colnames(positive_temp)[1]<-"positive"
  colnames(negative_temp)[1]<-"negative"
  # count number of columns positive/negative per row
  for (x in 1:dim(log_FC_results_subset)[1]){
    positive_temp[x,1]<-sum(log_FC_results_subset[x,]>0)
    negative_temp[x,1]<-sum(log_FC_results_subset[x,]<0)
  }
  #merge
  results<-cbind(log_FC_results_subset, positive_temp, negative_temp)
  # decide if gene is up/down regulated
  # if all columns for a row is negative == down
  # if all columns for a row is positive == up
  # if negative and positive values exist for a row == problem
  for (x in 1:dim(results)[1]){
    if ( results$positive[x]>0 & results$negative[x]==0 ) {
      results$outcome[x]<-"up"
    }
    else if ( results$positive[x]==0 & results$negative[x]>0 ) {
      results$outcome[x]<-"down"
    }
    else
      results$outcome[x]<-"problem"
  }
  print(table(results$outcome))
  return(results)
}

# apply function
Frontal_lobe_DEG_regulation<-calculate_AW_gene_regulation(Frontal_lobe_DE,log_FC_results)
Frontal_lobe_DEG_regulation[rownames(subset(Frontal_lobe_DEG_regulation, outcome=="problem")),]

Parietal_lobe_DEG_regulation<-calculate_AW_gene_regulation(Parietal_lobe_DE,log_FC_results)
Parietal_lobe_DEG_regulation[rownames(subset(Parietal_lobe_DEG_regulation, outcome=="problem")),]

Temporal_lobe_DEG_regulation<-calculate_AW_gene_regulation(Temporal_lobe_DE,log_FC_results)
Temporal_lobe_DEG_regulation[rownames(subset(Temporal_lobe_DEG_regulation, outcome=="problem")),]

Cerebellum_DEG_regulation<-calculate_AW_gene_regulation(Cerebellum_DE,log_FC_results)
Cerebellum_DEG_regulation[rownames(subset(Cerebellum_DEG_regulation, outcome=="problem")),]

# keep genes which are not problematic

Frontal_lobe_DEG_regulation_good_genes<-subset(Frontal_lobe_DEG_regulation, outcome!="problem")
Parietal_lobe_DEG_regulation_good_genes<-subset(Parietal_lobe_DEG_regulation, outcome!="problem")
Temporal_lobe_DEG_regulation_good_genes<-subset(Temporal_lobe_DEG_regulation, outcome!="problem")
Cerebellum_DEG_regulation_good_genes<-subset(Cerebellum_DEG_regulation, outcome!="problem")

dim(Frontal_lobe_DEG_regulation)
dim(Frontal_lobe_DEG_regulation_good_genes)

dim(Parietal_lobe_DEG_regulation)
dim(Parietal_lobe_DEG_regulation_good_genes)

dim(Temporal_lobe_DEG_regulation)
dim(Temporal_lobe_DEG_regulation_good_genes)

dim(Cerebellum_DEG_regulation)
dim(Cerebellum_DEG_regulation_good_genes)

# merge meta AW p value - FDR adjusted

Frontal_lobe_DEG_regulation_good_genes<-merge(Frontal_lobe_DEG_regulation_good_genes, Frontal_lobe_DE_AW.OC, by="row.names")
rownames(Frontal_lobe_DEG_regulation_good_genes)<-Frontal_lobe_DEG_regulation_good_genes$Row.names
Frontal_lobe_DEG_regulation_good_genes$Row.names<-NULL
dim(Frontal_lobe_DEG_regulation_good_genes)

Parietal_lobe_DEG_regulation_good_genes<-merge(Parietal_lobe_DEG_regulation_good_genes, Parietal_lobe_DE_AW.OC, by="row.names")
rownames(Parietal_lobe_DEG_regulation_good_genes)<-Parietal_lobe_DEG_regulation_good_genes$Row.names
Parietal_lobe_DEG_regulation_good_genes$Row.names<-NULL
dim(Parietal_lobe_DEG_regulation_good_genes)

Temporal_lobe_DEG_regulation_good_genes<-merge(Temporal_lobe_DEG_regulation_good_genes, Temporal_lobe_DE_AW.OC, by="row.names")
rownames(Temporal_lobe_DEG_regulation_good_genes)<-Temporal_lobe_DEG_regulation_good_genes$Row.names
Temporal_lobe_DEG_regulation_good_genes$Row.names<-NULL
dim(Temporal_lobe_DEG_regulation_good_genes)

Cerebellum_DEG_regulation_good_genes<-merge(Cerebellum_DEG_regulation_good_genes, Cerebellum_DE_AW.OC, by="row.names")
rownames(Cerebellum_DEG_regulation_good_genes)<-Cerebellum_DEG_regulation_good_genes$Row.names
Cerebellum_DEG_regulation_good_genes$Row.names<-NULL
dim(Cerebellum_DEG_regulation_good_genes)

# separate list into up/down regulated per brain region and write
setwd(AW_results_dir)

## FRONTAL LOBE

# extract up
Frontal_lobe_DEG_regulation_good_genes_up<-subset(Frontal_lobe_DEG_regulation_good_genes, outcome=="up")
# keep only AW column
Frontal_lobe_DEG_regulation_good_genes_up<-Frontal_lobe_DEG_regulation_good_genes_up[grep("AW.OC", colnames(Frontal_lobe_DEG_regulation_good_genes_up))]
#move rownames to entrez column
Frontal_lobe_DEG_regulation_good_genes_up$Entrez_Gene_ID<-rownames(Frontal_lobe_DEG_regulation_good_genes_up)
#re-arrange column- entrez id 1st then p val
Frontal_lobe_DEG_regulation_good_genes_up<-Frontal_lobe_DEG_regulation_good_genes_up[c(2,1)]
# rename columns
colnames(Frontal_lobe_DEG_regulation_good_genes_up)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Frontal_lobe_DEG_regulation_good_genes_up)
#write
write.table(Frontal_lobe_DEG_regulation_good_genes_up, file="Frontal_lobe_AW.OC_DEG_up", row.names=F, sep="\t", quote=F)

# extract down
Frontal_lobe_DEG_regulation_good_genes_down<-subset(Frontal_lobe_DEG_regulation_good_genes, outcome=="down")
# keep only AW column
Frontal_lobe_DEG_regulation_good_genes_down<-Frontal_lobe_DEG_regulation_good_genes_down[grep("AW.OC", colnames(Frontal_lobe_DEG_regulation_good_genes_down))]
#move rownames to entrez column
Frontal_lobe_DEG_regulation_good_genes_down$Entrez_Gene_ID<-rownames(Frontal_lobe_DEG_regulation_good_genes_down)
#re-arrange column- entrez id 1st then p val
Frontal_lobe_DEG_regulation_good_genes_down<-Frontal_lobe_DEG_regulation_good_genes_down[c(2,1)]
# rename columns
colnames(Frontal_lobe_DEG_regulation_good_genes_down)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Frontal_lobe_DEG_regulation_good_genes_down)
#write
write.table(Frontal_lobe_DEG_regulation_good_genes_down, file="Frontal_lobe_AW.OC_DEG_down", row.names=F, sep="\t", quote=F)


# check up/down add up to correct total

dim(Frontal_lobe_DEG_regulation_good_genes_up)[1]+dim(Frontal_lobe_DEG_regulation_good_genes_down)[1]==dim(Frontal_lobe_DEG_regulation_good_genes)[1]


## PARIETAL LOBE

# extract up
Parietal_lobe_DEG_regulation_good_genes_up<-subset(Parietal_lobe_DEG_regulation_good_genes, outcome=="up")
# keep only AW column
Parietal_lobe_DEG_regulation_good_genes_up<-Parietal_lobe_DEG_regulation_good_genes_up[grep("AW.OC", colnames(Parietal_lobe_DEG_regulation_good_genes_up))]
#move rownames to entrez column
Parietal_lobe_DEG_regulation_good_genes_up$Entrez_Gene_ID<-rownames(Parietal_lobe_DEG_regulation_good_genes_up)
#re-arrange column- entrez id 1st then p val
Parietal_lobe_DEG_regulation_good_genes_up<-Parietal_lobe_DEG_regulation_good_genes_up[c(2,1)]
# rename columns
colnames(Parietal_lobe_DEG_regulation_good_genes_up)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Parietal_lobe_DEG_regulation_good_genes_up)
#write
write.table(Parietal_lobe_DEG_regulation_good_genes_up, file="Parietal_lobe_AW.OC_DEG_up", row.names=F, sep="\t", quote=F)

# extract down
Parietal_lobe_DEG_regulation_good_genes_down<-subset(Parietal_lobe_DEG_regulation_good_genes, outcome=="down")
# keep only AW column
Parietal_lobe_DEG_regulation_good_genes_down<-Parietal_lobe_DEG_regulation_good_genes_down[grep("AW.OC", colnames(Parietal_lobe_DEG_regulation_good_genes_down))]
#move rownames to entrez column
Parietal_lobe_DEG_regulation_good_genes_down$Entrez_Gene_ID<-rownames(Parietal_lobe_DEG_regulation_good_genes_down)
#re-arrange column- entrez id 1st then p val
Parietal_lobe_DEG_regulation_good_genes_down<-Parietal_lobe_DEG_regulation_good_genes_down[c(2,1)]
# rename columns
colnames(Parietal_lobe_DEG_regulation_good_genes_down)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Parietal_lobe_DEG_regulation_good_genes_down)
#write
write.table(Parietal_lobe_DEG_regulation_good_genes_down, file="Parietal_lobe_AW.OC_DEG_down", row.names=F, sep="\t", quote=F)


# check up/down add up to correct total

dim(Parietal_lobe_DEG_regulation_good_genes_up)[1]+dim(Parietal_lobe_DEG_regulation_good_genes_down)[1]==dim(Parietal_lobe_DEG_regulation_good_genes)[1]


## TEMPORAL LOBE

# extract up
Temporal_lobe_DEG_regulation_good_genes_up<-subset(Temporal_lobe_DEG_regulation_good_genes, outcome=="up")
# keep only AW column
Temporal_lobe_DEG_regulation_good_genes_up<-Temporal_lobe_DEG_regulation_good_genes_up[grep("AW.OC", colnames(Temporal_lobe_DEG_regulation_good_genes_up))]
#move rownames to entrez column
Temporal_lobe_DEG_regulation_good_genes_up$Entrez_Gene_ID<-rownames(Temporal_lobe_DEG_regulation_good_genes_up)
#re-arrange column- entrez id 1st then p val
Temporal_lobe_DEG_regulation_good_genes_up<-Temporal_lobe_DEG_regulation_good_genes_up[c(2,1)]
# rename columns
colnames(Temporal_lobe_DEG_regulation_good_genes_up)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Temporal_lobe_DEG_regulation_good_genes_up)
#write
write.table(Temporal_lobe_DEG_regulation_good_genes_up, file="Temporal_lobe_AW.OC_DEG_up", row.names=F, sep="\t", quote=F)

# extract down
Temporal_lobe_DEG_regulation_good_genes_down<-subset(Temporal_lobe_DEG_regulation_good_genes, outcome=="down")
# keep only AW column
Temporal_lobe_DEG_regulation_good_genes_down<-Temporal_lobe_DEG_regulation_good_genes_down[grep("AW.OC", colnames(Temporal_lobe_DEG_regulation_good_genes_down))]
#move rownames to entrez column
Temporal_lobe_DEG_regulation_good_genes_down$Entrez_Gene_ID<-rownames(Temporal_lobe_DEG_regulation_good_genes_down)
#re-arrange column- entrez id 1st then p val
Temporal_lobe_DEG_regulation_good_genes_down<-Temporal_lobe_DEG_regulation_good_genes_down[c(2,1)]
# rename columns
colnames(Temporal_lobe_DEG_regulation_good_genes_down)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Temporal_lobe_DEG_regulation_good_genes_down)
#write
write.table(Temporal_lobe_DEG_regulation_good_genes_down, file="Temporal_lobe_AW.OC_DEG_down", row.names=F, sep="\t", quote=F)


# check up/down add up to correct total

dim(Temporal_lobe_DEG_regulation_good_genes_up)[1]+dim(Temporal_lobe_DEG_regulation_good_genes_down)[1]==dim(Temporal_lobe_DEG_regulation_good_genes)[1]


## CEREBELLUM

# extract up
Cerebellum_DEG_regulation_good_genes_up<-subset(Cerebellum_DEG_regulation_good_genes, outcome=="up")
# keep only AW column
Cerebellum_DEG_regulation_good_genes_up<-Cerebellum_DEG_regulation_good_genes_up[grep("AW.OC", colnames(Cerebellum_DEG_regulation_good_genes_up))]
#move rownames to entrez column
Cerebellum_DEG_regulation_good_genes_up$Entrez_Gene_ID<-rownames(Cerebellum_DEG_regulation_good_genes_up)
#re-arrange column- entrez id 1st then p val
Cerebellum_DEG_regulation_good_genes_up<-Cerebellum_DEG_regulation_good_genes_up[c(2,1)]
# rename columns
colnames(Cerebellum_DEG_regulation_good_genes_up)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Cerebellum_DEG_regulation_good_genes_up)
#write
write.table(Cerebellum_DEG_regulation_good_genes_up, file="Cerebellum_AW.OC_DEG_up", row.names=F, sep="\t", quote=F)

# extract down
Cerebellum_DEG_regulation_good_genes_down<-subset(Cerebellum_DEG_regulation_good_genes, outcome=="down")
# keep only AW column
Cerebellum_DEG_regulation_good_genes_down<-Cerebellum_DEG_regulation_good_genes_down[grep("AW.OC", colnames(Cerebellum_DEG_regulation_good_genes_down))]
#move rownames to entrez column
Cerebellum_DEG_regulation_good_genes_down$Entrez_Gene_ID<-rownames(Cerebellum_DEG_regulation_good_genes_down)
#re-arrange column- entrez id 1st then p val
Cerebellum_DEG_regulation_good_genes_down<-Cerebellum_DEG_regulation_good_genes_down[c(2,1)]
# rename columns
colnames(Cerebellum_DEG_regulation_good_genes_down)[2]<-"AW_FDR_adjusted_p_val"
#check
head(Cerebellum_DEG_regulation_good_genes_down)
#write
write.table(Cerebellum_DEG_regulation_good_genes_down, file="Cerebellum_AW.OC_DEG_down", row.names=F, sep="\t", quote=F)


# check up/down add up to correct total

dim(Cerebellum_DEG_regulation_good_genes_up)[1]+dim(Cerebellum_DEG_regulation_good_genes_down)[1]==dim(Cerebellum_DEG_regulation_good_genes)[1]

##### EXRACT ALL AW DEG RESULTS FOR EACH BRAIN REGION #####

#added 6/12/16
# for GSEA in CMap all genes are required from AW results

# ##### EXTRACT DEG PER BRAIN REGION ####
# 
# 
# 
# # using  ind.tail="high" and ind.tail="low" gives same number of DEG using all meta-analysis methods - using only up-genes object
# 
# Frontal_lobe_all_DE_AW.OC<-na.omit((as.data.frame(Frontal_lobe_DE$meta.analysis$FDR))[1])
# 
# Parietal_lobe_all_DE_AW.OC<-na.omit((as.data.frame(Parietal_lobe_DE$meta.analysis$FDR))[1])
# 
# Cerebellum_all_DE_AW.OC<-na.omit((as.data.frame(Cerebellum_DE$meta.analysis$FDR))[1])
# 
# Temporal_lobe_all_DE_AW.OC<-na.omit((as.data.frame(Temporal_lobe_DE$meta.analysis$FDR))[1])
# 
# dim(Temporal_lobe_all_DE_AW.OC)
# head(Temporal_lobe_all_DE_AW.OC)
# 
# dim(Frontal_lobe_all_DE_AW.OC)
# head(Frontal_lobe_all_DE_AW.OC)
# 
# dim(Parietal_lobe_all_DE_AW.OC)
# head(Parietal_lobe_all_DE_AW.OC)
# 
# dim(Cerebellum_all_DE_AW.OC)
# head(Cerebellum_all_DE_AW.OC)
# 
# # create function - 
# # 1.extract weights
# # 2.subset limma results to match studies in meta analysis
# # 3.according to AW weights - make any study where weight=0 to 0 in individual study stats
# # 4.keep only values that have contributed to AW.OC meta-analysis.
# # 5.make new column to count number of negative logFC values across studies per entrez gene id
# # 6.make new column to count number of positive logFC values across studies per entrez gene id 
# # 7.decide if genes are up/down regulated
# 
# calculate_gene_regulation_for_all_AW_results<-function(dataset, logFC_results){
#   #extract sig meta P values
#   sig_p<-(as.data.frame(subset(dataset$meta.analysis$FDR, dataset$meta.analysis$FDR[,1]<=0.05)))[1]
#   # extract meta-DEG list AW weights - add colnames to match individual studies
#   dataset_AW.OC_weight<-as.data.frame(dataset$meta.analysis$AW.weight[rownames(sig_p),])
#   colnames(dataset_AW.OC_weight)<-colnames(dataset$ind.stat)
#   #remove "_MetaOmics_format" from colnames
#   colnames(dataset_AW.OC_weight)<-sub("_MetaOmics_format", "", colnames(dataset_AW.OC_weight))
#   # remove "_data"
#   colnames(dataset_AW.OC_weight)<-sub("_data", "", colnames(dataset_AW.OC_weight))
#   # subset limma results to match dataset - columns
#   log_FC_results_subset<-log_FC_results[,colnames(log_FC_results) %in% colnames(dataset_AW.OC_weight)]
#   # subset limma rownames results to match dataset
#   log_FC_results_subset<-log_FC_results_subset[rownames(log_FC_results_subset) %in% rownames(dataset_AW.OC_weight),]
#   # sort limma results to match dataset order
#   log_FC_results_subset<-log_FC_results_subset[order(rownames(log_FC_results_subset)),]
#   dataset_AW.OC_weight<-dataset_AW.OC_weight[order(rownames(dataset_AW.OC_weight)),]
#   # check number of studies match
#   print(paste("Number of studies in limma results match meta analysis results: ", length(colnames(log_FC_results_subset))==length(colnames(dataset_AW.OC_weight)), sep=""))
#   # check studies are ordered correctly across both dataframes
#   print(paste("Studies correctly ordered by colnames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
#   # check rownames are same
#   print(paste("Studies correctly ordered by rownames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
#   # print number of genes and studies
#   print(paste("Number of Studies: ", dim(dataset_AW.OC_weight)[2], sep=""))
#   print(paste("Number of Genes: ", dim(dataset_AW.OC_weight)[1], sep=""))
#   # convert logFC to 0 if AW.OC weight==0
#   for (x in 1:dim(dataset_AW.OC_weight)[1]){
#     for (y in 1:dim(dataset_AW.OC_weight)[2])
#       if (dataset_AW.OC_weight[x,y]==0){
#         log_FC_results_subset[x,y]<-0
#       }
#   }
#   # make column down, up, problem according to logFC signs.
#   # count number of postive + negative values per row
#   # create temp dataframe with rownames
#   positive_temp<-log_FC_results_subset[1]
#   negative_temp<-log_FC_results_subset[1]
#   # rename columns
#   colnames(positive_temp)[1]<-"positive"
#   colnames(negative_temp)[1]<-"negative"
#   # count number of columns positive/negative per row
#   for (x in 1:dim(log_FC_results_subset)[1]){
#     positive_temp[x,1]<-sum(log_FC_results_subset[x,]>0)
#     negative_temp[x,1]<-sum(log_FC_results_subset[x,]<0)
#   }
#   #merge
#   results<-cbind(log_FC_results_subset, positive_temp, negative_temp)
#   # decide if gene is up/down regulated
#   # if all columns for a row is negative == down
#   # if all columns for a row is positive == up
#   # if negative and positive values exist for a row == problem
#   for (x in 1:dim(results)[1]){
#     if ( results$positive[x]>0 & results$negative[x]==0 ) {
#       results$outcome[x]<-"up"
#     }
#     else if ( results$positive[x]==0 & results$negative[x]>0 ) {
#       results$outcome[x]<-"down"
#     }
#     else
#       results$outcome[x]<-"problem"
#   }
#   print(table(results$outcome))
#   return(results)
# }
# 
# # apply function
# Frontal_lobe_all_DEG_regulation<-calculate_gene_regulation_for_all_AW_results(Frontal_lobe_DE,log_FC_results)
# Frontal_lobe_all_DEG_regulation[rownames(subset(Frontal_lobe_all_DEG_regulation, outcome=="problem")),]
# 
# Parietal_lobe_all_DEG_regulation<-calculate_gene_regulation_for_all_AW_results(Parietal_lobe_DE,log_FC_results)
# Parietal_lobe_all_DEG_regulation[rownames(subset(Parietal_lobe_all_DEG_regulation, outcome=="problem")),]
# 
# Temporal_lobe_all_DEG_regulation<-calculate_gene_regulation_for_all_AW_results(Temporal_lobe_DE,log_FC_results)
# Temporal_lobe_all_DEG_regulation[rownames(subset(Temporal_lobe_all_DEG_regulation, outcome=="problem")),]
# 
# Cerebellum_all_DEG_regulation<-calculate_gene_regulation_for_all_AW_results(Cerebellum_DE,log_FC_results)
# Cerebellum_all_DEG_regulation[rownames(subset(Cerebellum_all_DEG_regulation, outcome=="problem")),]
# 
# # keep genes which are not problematic
# 
# Frontal_lobe_all_DEG_regulation_good_genes<-subset(Frontal_lobe_all_DEG_regulation, outcome!="problem")
# Parietal_lobe_all_DEG_regulation_good_genes<-subset(Parietal_lobe_all_DEG_regulation, outcome!="problem")
# Temporal_lobe_all_DEG_regulation_good_genes<-subset(Temporal_lobe_all_DEG_regulation, outcome!="problem")
# Cerebellum_all_DEG_regulation_good_genes<-subset(Cerebellum_all_DEG_regulation, outcome!="problem")
# 
# dim(Frontal_lobe_all_DEG_regulation)
# dim(Frontal_lobe_all_DEG_regulation_good_genes)
# 
# dim(Parietal_lobe_all_DEG_regulation)
# dim(Parietal_lobe_all_DEG_regulation_good_genes)
# 
# dim(Temporal_lobe_all_DEG_regulation)
# dim(Temporal_lobe_all_DEG_regulation_good_genes)
# 
# dim(Cerebellum_all_DEG_regulation)
# dim(Cerebellum_all_DEG_regulation_good_genes)
# 
# # merge meta AW p value - FDR adjusted
# 
# Frontal_lobe_all_DEG_regulation_good_genes<-merge(Frontal_lobe_all_DEG_regulation_good_genes, Frontal_lobe_all_DE_AW.OC, by="row.names")
# rownames(Frontal_lobe_all_DEG_regulation_good_genes)<-Frontal_lobe_all_DEG_regulation_good_genes$Row.names
# Frontal_lobe_all_DEG_regulation_good_genes$Row.names<-NULL
# dim(Frontal_lobe_all_DEG_regulation_good_genes)
# 
# Parietal_lobe_all_DEG_regulation_good_genes<-merge(Parietal_lobe_all_DEG_regulation_good_genes, Parietal_lobe_all_DE_AW.OC, by="row.names")
# rownames(Parietal_lobe_all_DEG_regulation_good_genes)<-Parietal_lobe_all_DEG_regulation_good_genes$Row.names
# Parietal_lobe_all_DEG_regulation_good_genes$Row.names<-NULL
# dim(Parietal_lobe_all_DEG_regulation_good_genes)
# 
# Temporal_lobe_all_DEG_regulation_good_genes<-merge(Temporal_lobe_all_DEG_regulation_good_genes, Temporal_lobe_all_DE_AW.OC, by="row.names")
# rownames(Temporal_lobe_all_DEG_regulation_good_genes)<-Temporal_lobe_all_DEG_regulation_good_genes$Row.names
# Temporal_lobe_all_DEG_regulation_good_genes$Row.names<-NULL
# dim(Temporal_lobe_all_DEG_regulation_good_genes)
# 
# Cerebellum_all_DEG_regulation_good_genes<-merge(Cerebellum_all_DEG_regulation_good_genes, Cerebellum_all_DE_AW.OC, by="row.names")
# rownames(Cerebellum_all_DEG_regulation_good_genes)<-Cerebellum_all_DEG_regulation_good_genes$Row.names
# Cerebellum_all_DEG_regulation_good_genes$Row.names<-NULL
# dim(Cerebellum_all_DEG_regulation_good_genes)
# 
# ## FRONTAL LOBE - write Entrez ID, LogFC, FDR p val
# 
# # move rownames to Entrez gene column
# Frontal_lobe_all_DEG_regulation_good_genes$Entrez_Gene_ID<-rownames(Frontal_lobe_all_DEG_regulation_good_genes)
# # keep only entrez, outcome and AW.OC column
# Frontal_lobe_all_DEG_regulation_good_genes<-Frontal_lobe_all_DEG_regulation_good_genes[grep("Entrez_Gene_ID|outcome|AW.OC", colnames(Frontal_lobe_all_DEG_regulation_good_genes))]
# #sort by colnames
# Frontal_lobe_all_DEG_regulation_good_genes<-Frontal_lobe_all_DEG_regulation_good_genes[,order(names(Frontal_lobe_all_DEG_regulation_good_genes))]
# Frontal_lobe_all_DEG_regulation_good_genes<-Frontal_lobe_all_DEG_regulation_good_genes[c(2,3,1)]
# # change colnames
# colnames(Frontal_lobe_all_DEG_regulation_good_genes)<-c("Entrez_Gene_ID", "Regulation", "AW_FDR_adjusted_p_val")
# head(Frontal_lobe_all_DEG_regulation_good_genes)
# dim(Frontal_lobe_all_DEG_regulation_good_genes)
# 
# setwd(Frontal_Lobe_dir)
# write.table(Frontal_lobe_all_DEG_regulation_good_genes, file="Frontal_lobe_AW.OC_DEG.txt", row.names=F, sep="\t", quote=F)
# 
# ## PARIETAL LOBE
# 
# # move rownames to Entrez gene column
# Parietal_lobe_all_DEG_regulation_good_genes$Entrez_Gene_ID<-rownames(Parietal_lobe_all_DEG_regulation_good_genes)
# # keep only entrez, outcome and AW.OC column
# Parietal_lobe_all_DEG_regulation_good_genes<-Parietal_lobe_all_DEG_regulation_good_genes[grep("Entrez_Gene_ID|outcome|AW.OC", colnames(Parietal_lobe_all_DEG_regulation_good_genes))]
# #sort by colnames
# Parietal_lobe_all_DEG_regulation_good_genes<-Parietal_lobe_all_DEG_regulation_good_genes[,order(names(Parietal_lobe_all_DEG_regulation_good_genes))]
# Parietal_lobe_all_DEG_regulation_good_genes<-Parietal_lobe_all_DEG_regulation_good_genes[c(2,3,1)]
# # change colnames
# colnames(Parietal_lobe_all_DEG_regulation_good_genes)<-c("Entrez_Gene_ID", "Regulation", "AW_FDR_adjusted_p_val")
# head(Parietal_lobe_all_DEG_regulation_good_genes)
# dim(Parietal_lobe_all_DEG_regulation_good_genes)
# 
# setwd(Parietal_Lobe_dir)
# write.table(Parietal_lobe_all_DEG_regulation_good_genes, file="Parietal_lobe_AW.OC_DEG.txt", row.names=F, sep="\t", quote=F)
# 
# 
# ## TEMPORAL LOBE
# 
# # move rownames to Entrez gene column
# Temporal_lobe_all_DEG_regulation_good_genes$Entrez_Gene_ID<-rownames(Temporal_lobe_all_DEG_regulation_good_genes)
# # keep only entrez, outcome and AW.OC column
# Temporal_lobe_all_DEG_regulation_good_genes<-Temporal_lobe_all_DEG_regulation_good_genes[grep("Entrez_Gene_ID|outcome|AW.OC", colnames(Temporal_lobe_all_DEG_regulation_good_genes))]
# #sort by colnames
# Temporal_lobe_all_DEG_regulation_good_genes<-Temporal_lobe_all_DEG_regulation_good_genes[,order(names(Temporal_lobe_all_DEG_regulation_good_genes))]
# Temporal_lobe_all_DEG_regulation_good_genes<-Temporal_lobe_all_DEG_regulation_good_genes[c(2,3,1)]
# # change colnames
# colnames(Temporal_lobe_all_DEG_regulation_good_genes)<-c("Entrez_Gene_ID", "Regulation", "AW_FDR_adjusted_p_val")
# head(Temporal_lobe_all_DEG_regulation_good_genes)
# dim(Temporal_lobe_all_DEG_regulation_good_genes)
# 
# setwd(Temporal_Lobe_dir)
# write.table(Temporal_lobe_all_DEG_regulation_good_genes, file="Temporal_lobe_AW.OC_DEG.txt", row.names=F, sep="\t", quote=F)
# 
# ## CEREBELLUM
# 
# # move rownames to Entrez gene column
# Cerebellum_all_DEG_regulation_good_genes$Entrez_Gene_ID<-rownames(Cerebellum_all_DEG_regulation_good_genes)
# # keep only entrez, outcome and AW.OC column
# Cerebellum_all_DEG_regulation_good_genes<-Cerebellum_all_DEG_regulation_good_genes[grep("Entrez_Gene_ID|outcome|AW.OC", colnames(Cerebellum_all_DEG_regulation_good_genes))]
# #sort by colnames
# Cerebellum_all_DEG_regulation_good_genes<-Cerebellum_all_DEG_regulation_good_genes[,order(names(Cerebellum_all_DEG_regulation_good_genes))]
# Cerebellum_all_DEG_regulation_good_genes<-Cerebellum_all_DEG_regulation_good_genes[c(2,3,1)]
# # change colnames
# colnames(Cerebellum_all_DEG_regulation_good_genes)<-c("Entrez_Gene_ID", "Regulation", "AW_FDR_adjusted_p_val")
# head(Cerebellum_all_DEG_regulation_good_genes)
# dim(Cerebellum_all_DEG_regulation_good_genes)
# 
# setwd(Cerebellum_dir)
# write.table(Cerebellum_all_DEG_regulation_good_genes, file="Cerebellum_AW.OC_DEG.txt", row.names=F, sep="\t", quote=F)

##### EXTRACT DEG PER BRAIN REGION - ADD META SUMMARIES (LOGFC) ####

# 6/4/17 modified to include logFC results for SPIA pathway analysis

# raw AW results
dim(Temporal_lobe_DE_AW.OC)
head(Temporal_lobe_DE_AW.OC)

dim(Frontal_lobe_DE_AW.OC)
head(Frontal_lobe_DE_AW.OC)

dim(Parietal_lobe_DE_AW.OC)
head(Parietal_lobe_DE_AW.OC)

dim(Cerebellum_DE_AW.OC)
head(Cerebellum_DE_AW.OC)

#read limma results with SE
limma_results_with_SE<-read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/3.AD_single_dataset_DE_analysis/Differential_Expression/Differential_expression_on_full_expression_data_AD_vs_CONTROL_with_SE.txt",
                                  header=T,
                                  as.is=T,
                                  sep="\t")

head(limma_results_with_SE)[1:10]

# extract fold change and SE
new_limma_logFC<-limma_results_with_SE[,grep("_logFC", colnames(limma_results_with_SE))]
new_limma_SE<-limma_results_with_SE[,grep("_SE", colnames(limma_results_with_SE))]

head(new_limma_logFC)[1:5]
head(new_limma_SE)[1.:5]

# remove logFC name from colnames
colnames(new_limma_logFC)<-sub("_logFC", "", colnames(new_limma_logFC))
colnames(new_limma_SE)<-sub("_SE", "", colnames(new_limma_SE))

all(colnames(new_limma_SE)==colnames(new_limma_logFC))==T


# create function - 
# 1.extract weights
# 2.subset limma results to match studies in meta analysis
# 3.according to AW weights - make any study where weight=0 to 0 in individual study stats
# 4.keep only values that have contributed to AW.OC meta-analysis.
# 5.make new column to count number of negative logFC values across studies per entrez gene id
# 6.make new column to count number of positive logFC values across studies per entrez gene id 
# 7.calculate meta summaries from logfc and SE wusing rmeta package


calculate_AW_gene_regulation_with_logFC<-function(dataset, AW_results){
  #extract sig meta P values
  sig_p<-(as.data.frame(subset(dataset$meta.analysis$FDR, dataset$meta.analysis$FDR[,1]<=0.05)))[1]
  # extract meta-DEG list AW weights - add colnames to match individual studies
  dataset_AW.OC_weight<-as.data.frame(dataset$meta.analysis$AW.weight[rownames(sig_p),])
  colnames(dataset_AW.OC_weight)<-colnames(dataset$ind.stat)
  #remove "_MetaOmics_format" from colnames
  colnames(dataset_AW.OC_weight)<-sub("_MetaOmics_format", "", colnames(dataset_AW.OC_weight))
  # remove "_data"
  colnames(dataset_AW.OC_weight)<-sub("_data", "", colnames(dataset_AW.OC_weight))
  # subset limma results to match dataset - columns
  log_FC_results_subset<-new_limma_logFC[,colnames(new_limma_logFC) %in% colnames(dataset_AW.OC_weight)]
  # subset limma logFC rownames results to match dataset
  log_FC_results_subset<-log_FC_results_subset[rownames(log_FC_results_subset) %in% rownames(dataset_AW.OC_weight),]
  #subset limma SE
  SE_results_subset<-new_limma_SE[,colnames(new_limma_SE) %in% colnames(dataset_AW.OC_weight)]
  # subset limma logFC rownames results to match dataset
  SE_results_subset<-SE_results_subset[rownames(SE_results_subset) %in% rownames(dataset_AW.OC_weight),]
  # sort limma results to match dataset order
  log_FC_results_subset<-log_FC_results_subset[order(rownames(log_FC_results_subset)),]
  SE_results_subset<-SE_results_subset[order(rownames(SE_results_subset)),]
  dataset_AW.OC_weight<-dataset_AW.OC_weight[order(rownames(dataset_AW.OC_weight)),]
  # check limma results for logFC + SE orered corectly
  print(paste("Limma colnames same: ", all(colnames(log_FC_results_subset)==colnames(SE_results_subset))==T, sep=""))
  print(paste("Limma rownames same: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
  # check number of studies match
  print(paste("Number of studies in limma results match meta analysis results: ", length(colnames(log_FC_results_subset))==length(colnames(dataset_AW.OC_weight)), sep=""))
  # check studies are ordered correctly across both dataframes
  print(paste("Studies correctly ordered by colnames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
  # check rownames are same
  print(paste("Studies correctly ordered by rownames: ", all(colnames(log_FC_results_subset)==colnames(dataset_AW.OC_weight))==T, sep=""))
  
  
  # print number of genes and studies
  print(paste("Number of Studies: ", dim(dataset_AW.OC_weight)[2], sep=""))
  print(paste("Number of Genes: ", dim(dataset_AW.OC_weight)[1], sep=""))
  # convert logFC + SE to 0 if AW.OC weight==0
  for (x in 1:dim(dataset_AW.OC_weight)[1]){
    for (y in 1:dim(dataset_AW.OC_weight)[2])
      if (dataset_AW.OC_weight[x,y]==0){
        log_FC_results_subset[x,y]<-0
        SE_results_subset[x,y]<-0
      }
  }
  # make column down, up, problem according to logFC signs.
  # count number of postive + negative values per row
  # create temp dataframe with rownames
  positive_temp<-log_FC_results_subset[1]
  negative_temp<-log_FC_results_subset[1]
  # rename columns
  colnames(positive_temp)[1]<-"positive"
  colnames(negative_temp)[1]<-"negative"
  # count number of columns positive/negative per row
  for (x in 1:dim(log_FC_results_subset)[1]){
    positive_temp[x,1]<-sum(log_FC_results_subset[x,]>0)
    negative_temp[x,1]<-sum(log_FC_results_subset[x,]<0)
  }
  # add dummy summaries column
  log_FC_results_subset$meta_summary<-0
  # loop through logFC results and add meta summaries
  for (x in 1:nrow(log_FC_results_subset)) {
    # vector with logFC + remove zero
    d<-as.matrix(log_FC_results_subset[x,])
    d<-as.vector(d[d!=0])
    # vector with SE + remove zero
    se<-as.matrix(SE_results_subset[x,])
    se<-as.vector(se[se!=0])
    # meta summaries
    log_FC_results_subset$meta_summary[x]<-(meta.summaries(d, se, method="fixed"))$summary
  }
  #merge
  results<-cbind(log_FC_results_subset, positive_temp, negative_temp)
  # decide if gene is up/down regulated
  # if all columns for a row is negative == down
  # if all columns for a row is positive == up
  # if negative and positive values exist for a row == problem
  for (x in 1:dim(results)[1]){
    if ( results$positive[x]>0 & results$negative[x]==0 ) {
      results$outcome[x]<-"up"
    }
    else if ( results$positive[x]==0 & results$negative[x]>0 ) {
      results$outcome[x]<-"down"
    }
    else
      results$outcome[x]<-"problem"
  }
  print(table(results$outcome))
  # keep non problematic genes only
  results<-subset(results, outcome!="problem")
  #merge with AW results
  results_with_AW<-merge(results, AW_results, by="row.names")
  rownames(results_with_AW)<-results_with_AW$Row.names
  results_with_AW$Row.names<-NULL
  # keep only p value and logFC
  results_with_AW<-results_with_AW[grep("AW.OC|meta_summary", colnames(results_with_AW))]
  #move rownames to entrez column
  results_with_AW$Entrez_Gene_ID<-rownames(results_with_AW)
  #re-arrange column- entrez id 1st then logFC then Pval
  results_with_AW<-results_with_AW[c(3,1,2)]
  # rename columns
  colnames(results_with_AW)[3]<-"AW_FDR_adjusted_p_val"
  #rename rownames
  rownames(results_with_AW)<-1:nrow(results_with_AW)
  return(results_with_AW)
}

# apply function
Frontal_lobe_sig_AW_results_with_logFC<-calculate_AW_gene_regulation_with_logFC(Frontal_lobe_DE, Frontal_lobe_DE_AW.OC)
head(Frontal_lobe_sig_AW_results_with_logFC)
tail(Frontal_lobe_sig_AW_results_with_logFC)
#compare to original results
dim(Frontal_lobe_sig_AW_results_with_logFC)
dim(Frontal_lobe_DEG_regulation_good_genes)


Parietal_lobe_sig_AW_results_with_logFC<-calculate_AW_gene_regulation_with_logFC(Parietal_lobe_DE, Parietal_lobe_DE_AW.OC)
head(Parietal_lobe_sig_AW_results_with_logFC)
tail(Parietal_lobe_sig_AW_results_with_logFC)
#compare to original results
dim(Parietal_lobe_sig_AW_results_with_logFC)
dim(Parietal_lobe_DEG_regulation_good_genes)


Temporal_lobe_sig_AW_results_with_logFC<-calculate_AW_gene_regulation_with_logFC(Temporal_lobe_DE, Temporal_lobe_DE_AW.OC)
head(Temporal_lobe_sig_AW_results_with_logFC)
tail(Temporal_lobe_sig_AW_results_with_logFC)
#compare to original results
dim(Temporal_lobe_sig_AW_results_with_logFC)
dim(Temporal_lobe_DEG_regulation_good_genes)


Cerebellum_DEG_sig_AW_results_with_logFC<-calculate_AW_gene_regulation_with_logFC(Cerebellum_DE, Cerebellum_DE_AW.OC)
head(Cerebellum_DEG_sig_AW_results_with_logFC)
tail(Cerebellum_DEG_sig_AW_results_with_logFC)
#compare to original results
dim(Cerebellum_DEG_sig_AW_results_with_logFC)
dim(Cerebellum_DEG_regulation_good_genes)

# write results

setwd(AW_results_logfc_dir)

write.table(Frontal_lobe_sig_AW_results_with_logFC, file="Frontal_lobe_sig_AW_results_with_logFC.txt", row.names=F, sep="\t", quote=F)
write.table(Parietal_lobe_sig_AW_results_with_logFC, file="Parietal_lobe_sig_AW_results_with_logFC.txt", row.names=F, sep="\t", quote=F)
write.table(Temporal_lobe_sig_AW_results_with_logFC, file="Temporal_lobe_sig_AW_results_with_logFC.txt", row.names=F, sep="\t", quote=F)
write.table(Cerebellum_DEG_sig_AW_results_with_logFC, file="Cerebellum_DEG_sig_AW_results_with_logFC.txt", row.names=F, sep="\t", quote=F)

##### META PATHWAY ANALYSIS #####

setwd(work_dir)

# pathway database - file downloaded manually from website ("http://software.broadinstitute.org/gsea/downloads.jsp")
pathway_database<-getGmt("c2.all.v5.1.entrez.gmt") 

setwd(Meta_Pathway_dir)

# FRONTAL LOBE
# create empty list 
Frontal_Lobe_pathwaydata<-vector(mode ="list", length = length(Frontal_lobe_no_AMP_QC_filtered))

#make data in format for pathway analysis
for(t1 in 1:length(Frontal_lobe_no_AMP_QC_filtered)) {
  temp<-impute.knn(Frontal_lobe_no_AMP_QC_filtered[[t1]][[1]])$data
  Frontal_Lobe_pathwaydata[[t1]]=list(x=temp, y=Frontal_lobe_no_AMP_QC_filtered[[t1]][[2]],
                                      geneid=rownames(Frontal_lobe_no_AMP_QC_filtered[[1]][[1]]),
                                      samplename=paste("sample", 1:ncol(Frontal_lobe_no_AMP_QC_filtered[[t1]][[1]]), sep=""))
}

# run pathway analysis
Frontal_Lobe_pathway_MAPE<-MAPE(arraydata=Frontal_Lobe_pathwaydata, # data
                                pathway.DB=pathway_database, # pathway data
                                resp.type="twoclass", # case vs control
                                stat="Fisher", # minP, maxP, rth, Fisher
                                nperm=300, # n.o permutations
                                permutation="sample", # "sample" or "gene" permutation
                                size.min=5, # min pathway size considered
                                size.max=500) # max pathway size considered

head(Frontal_Lobe_pathway_MAPE)

# check significant pathways
subset(Frontal_Lobe_pathway_MAPE$qvalue, MAPE_I<0.9)

#heatmap of pathways
plotMAPE(Frontal_Lobe_pathway_MAPE, 
         cutoff=0.865, 
         MAPE.method="MAPE_I")

# PARIETAL LOBE
# create empty list 
Parietal_Lobe_pathwaydata<-vector(mode ="list", length = length(Parietal_lobe_no_AMP_QC_filtered))

#make data in format for pathway analysis
for(t1 in 1:length(Parietal_lobe_no_AMP_QC_filtered)) {
  temp<-impute.knn(Parietal_lobe_no_AMP_QC_filtered[[t1]][[1]])$data
  Parietal_Lobe_pathwaydata[[t1]]=list(x=temp, y=Parietal_lobe_no_AMP_QC_filtered[[t1]][[2]],
                                       geneid=rownames(Parietal_lobe_no_AMP_QC_filtered[[1]][[1]]),
                                       samplename=paste("sample", 1:ncol(Parietal_lobe_no_AMP_QC_filtered[[t1]][[1]]), sep=""))
}

# run pathway analysis
Parietal_Lobe_pathway_MAPE<-MAPE(arraydata=Parietal_Lobe_pathwaydata, # data
                                 pathway.DB=pathway_database, # pathway data
                                 resp.type="twoclass", # case vs control
                                 stat="Fisher", # minP, maxP, rth, Fisher
                                 nperm=300, # n.o permutations
                                 permutation="sample", # "sample" or "gene" permutation
                                 size.min=5, # min pathway size considered
                                 size.max=500) # max pathway size considered

head(Parietal_Lobe_pathway_MAPE)

# check significant pathways
subset(Parietal_Lobe_pathway_MAPE$qvalue, MAPE_I<1)

#heatmap of pathways
plotMAPE(Parietal_Lobe_pathway_MAPE, 
         cutoff=1, 
         MAPE.method="MAPE_I")

# CEREBELLUM
# create empty list 
Cerebellum_pathwaydata<-vector(mode ="list", length = length(Cerebellum_no_AMP_QC_filtered))

#make data in format for pathway analysis
for(t1 in 1:length(Cerebellum_no_AMP_QC_filtered)) {
  temp<-impute.knn(Cerebellum_no_AMP_QC_filtered[[t1]][[1]])$data
  Cerebellum_pathwaydata[[t1]]=list(x=temp, y=Cerebellum_no_AMP_QC_filtered[[t1]][[2]],
                                    geneid=rownames(Cerebellum_no_AMP_QC_filtered[[1]][[1]]),
                                    samplename=paste("sample", 1:ncol(Cerebellum_no_AMP_QC_filtered[[t1]][[1]]), sep=""))
}

# run pathway analysis
Cerebellum_pathway_MAPE<-MAPE(arraydata=Cerebellum_pathwaydata, # data
                              pathway.DB=pathway_database, # pathway data
                              resp.type="twoclass", # case vs control
                              stat="Fisher", # minP, maxP, rth, Fisher
                              nperm=300, # n.o permutations
                              permutation="sample", # "sample" or "gene" permutation
                              size.min=5, # min pathway size considered
                              size.max=500) # max pathway size considered

head(Cerebellum_pathway_MAPE)

# check significant pathways
subset(Cerebellum_pathway_MAPE$qvalue, MAPE_I<1)
summary(Cerebellum_pathway_MAPE$qvalue[5])
summary(Cerebellum_pathway_MAPE$qvalue[4])


#heatmap of pathways
plotMAPE(Cerebellum_pathway_MAPE, 
         cutoff=1, 
         MAPE.method="MAPE_I")

# TEMPORAL LOBE
# create empty list 
Temporal_Lobe_pathwaydata<-vector(mode ="list", length = length(Temporal_lobe_no_AMP_QC_filtered))

#make data in format for pathway analysis
for(t1 in 1:length(Temporal_lobe_no_AMP_QC_filtered)) {
  temp<-impute.knn(Temporal_lobe_no_AMP_QC_filtered[[t1]][[1]])$data
  Temporal_Lobe_pathwaydata[[t1]]=list(x=temp, y=Temporal_lobe_no_AMP_QC_filtered[[t1]][[2]],
                                       geneid=rownames(Temporal_lobe_no_AMP_QC_filtered[[1]][[1]]),
                                       samplename=paste("sample", 1:ncol(Temporal_lobe_no_AMP_QC_filtered[[t1]][[1]]), sep=""))
}

# run pathway analysis
Temporal_Lobe_pathway_MAPE<-MAPE(arraydata=Temporal_Lobe_pathwaydata, # data
                                 pathway.DB=pathway_database, # pathway data
                                 resp.type="twoclass", # case vs control
                                 stat="Fisher", # minP, maxP, rth, Fisher
                                 nperm=300, # n.o permutations
                                 permutation="sample", # "sample" or "gene" permutation
                                 size.min=5, # min pathway size considered
                                 size.max=500) # max pathway size considered

head(Temporal_Lobe_pathway_MAPE)

# check significant pathways
dim(subset(Temporal_Lobe_pathway_MAPE$qvalue, MAPE_I<0.9))
dim(subset(Temporal_Lobe_pathway_MAPE$qvalue, MAPE_G<0.6))
dim(subset(Temporal_Lobe_pathway_MAPE$qvalue, MAPE_P<0.6))

#heatmap of pathways
plotMAPE(Temporal_Lobe_pathway_MAPE, 
         cutoff=0.9, 
         MAPE.method="MAPE_I")

subset(Temporal_Lobe_pathway_MAPE$qvalue, MAPE_I<0.9)

##### COMPILE REPORT #####

#rmarkdown::render("Meta_Analysis.R") # make sure this is hashed out and file is saved. run separetly on command line

##### SAVE IMAGE #####

setwd(work_dir)

#save.image("Meta_analysis.Rdata")
