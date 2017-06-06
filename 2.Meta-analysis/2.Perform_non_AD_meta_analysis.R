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
## DATE: 20/02/2017

##### DESCRIPTION OF ANALYSIS ####
## grouping non AD-neurodegenerative datasets into brain regions + diseases and running AW.OC meta analysis
## running MetaOmics - for DE 
##
##
## NOTE: latest gene set files files can be downloaded from http://software.broadinstitute.org/gsea/downloads.jsp
## 
## NO GENE ENTREZ ID CONVERSION TO GENE SYMBOL AS PACKAGE HAS ISSUES
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Clean_Data/"

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/5.non_AD_DE_Meta_Analysis"

setwd(work_dir)

# create dir DE and patway analysis

# DE directory
dir.create(paste(work_dir,"Meta_DE", sep="/"))
Meta_DE_dir=paste(work_dir,"Meta_DE", sep="/")

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

##### READ DATA - 21 datasets #####

setwd(data_dir)

E_GEOD_12649_BD_Frontal_Lobe<-read.table(file="E-GEOD-12649_Bipolar_Disorder_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_12649_Schizophrenia_Frontal_Lobe<-read.table(file="E-GEOD-12649_Schizophrenia_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_17612_Schizophrenia_Frontal_Lobe<-read.table(file="E-GEOD-17612_Schizophrenia_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_20168_PD_Frontal_Lobe<-read.table(file="E-GEOD-20168_PD_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_21138_Schizophrenia_Frontal_Lobe<-read.table(file="E-GEOD-21138_Schizophrenia_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_21935_Schizophrenia_Temporal_Lobe<-read.table(file="E-GEOD-21935_Schizophrenia_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_35978_BD_Cerebellum<-read.table(file="E-GEOD-35978_Bipolar_Disorder_Cerebellum_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_35978_BD_Parietal_Lobe<-read.table(file="E-GEOD-35978_Bipolar_Disorder_Parietal_Lobe_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_35978_Schizophrenia_Cerebellum<-read.table(file="E-GEOD-35978_Cerebellum_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_35978_Schizophrenia_Parietal_Lobe<-read.table(file="E-GEOD-35978_Parietal_lobe_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_3790_HD_Cerebellum_Affy_U133A<-read.table(file="E-GEOD-3790_HD_Affy_U133A_Cerebellum_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A<-read.table(file="E-GEOD-3790_HD_Affy_U133A_Frontal_Lobe_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_3790_HD_Cerebellum_Affy_U133B<-read.table(file="E-GEOD-3790_HD_Affy_U133B_Cerebellum_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B<-read.table(file="E-GEOD-3790_HD_Affy_U133B_Frontal_Lobe_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_5388_BD_Frontal_Lobe<-read.table(file="E-GEOD-5388_Bipolar_Disorder_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_BD_Temporal_Lobe<-read.table(file="E-GEOD-53987_Bipolar_Disorder_Hippocampus_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_BD_Frontal_Lobe<-read.table(file="E-GEOD-53987_Bipolar_Disorder_Prefrontal_cortex_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_MDD_Temporal_Lobe<-read.table(file="E-GEOD-53987_MDD_Hippocampus_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_MDD_Frontal_Lobe<-read.table(file="E-GEOD-53987_MDD_Prefrontal_cortex_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_Schizophrenia_Frontal_Lobe<-read.table(file="E-GEOD-53987_Schizophrenia_Prefrontal_cortex_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_53987_Schizophrenia_Temporal_Lobe<-read.table(file="E-GEOD-53987_Scizophrenia_Hippocampus_processed_data.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)

setwd(work_dir)

head(E_GEOD_12649_BD_Frontal_Lobe)[1:10]
head(E_GEOD_12649_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_17612_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_20168_PD_Frontal_Lobe)[1:10]
head(E_GEOD_21138_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_21935_Schizophrenia_Temporal_Lobe)[1:10]
head(E_GEOD_35978_Schizophrenia_Cerebellum)[1:10]
head(E_GEOD_35978_Schizophrenia_Parietal_Lobe)[1:10]
head(E_GEOD_35978_BD_Cerebellum)[1:10]
head(E_GEOD_35978_BD_Parietal_Lobe)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133B)[1:10]
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B)[1:10]
head(E_GEOD_5388_BD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_BD_Temporal_Lobe)[1:10]
head(E_GEOD_53987_BD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_MDD_Temporal_Lobe)[1:10]
head(E_GEOD_53987_MDD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Temporal_Lobe)[1:10]

##### GROUP BRAIN REGIONS - total 21 datasets #####

# FRONTAL LOBE
head(E_GEOD_12649_BD_Frontal_Lobe)[1:10]
head(E_GEOD_12649_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_17612_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_20168_PD_Frontal_Lobe)[1:10]
head(E_GEOD_21138_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B)[1:10]
head(E_GEOD_5388_BD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_BD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_MDD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Frontal_Lobe)[1:10]

# TEMPORAL LOBE
head(E_GEOD_21935_Schizophrenia_Temporal_Lobe)[1:10]
head(E_GEOD_53987_BD_Temporal_Lobe)[1:10]
head(E_GEOD_53987_MDD_Temporal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Temporal_Lobe)[1:10]

# CEREBELLUM
head(E_GEOD_35978_BD_Cerebellum)[1:10]
head(E_GEOD_35978_Schizophrenia_Cerebellum)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133B)[1:10]

# PARIETAL LOBE
head(E_GEOD_35978_BD_Parietal_Lobe)[1:10]
head(E_GEOD_35978_Schizophrenia_Parietal_lobe)[1:10]

##### GROUP BY DISEASE #####

# BIPOLAR DISORDER
head(E_GEOD_12649_BD_Frontal_Lobe)[1:10]
head(E_GEOD_5388_BD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_BD_Frontal_Lobe)[1:10]
head(E_GEOD_35978_BD_Cerebellum)[1:10]
head(E_GEOD_35978_BD_Parietal_Lobe)[1:10]
head(E_GEOD_53987_BD_Temporal_Lobe)[1:10]

# SCHIZOPHRENIA
head(E_GEOD_12649_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_17612_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_21138_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Frontal_Lobe)[1:10]
head(E_GEOD_21935_Schizophrenia_Temporal_Lobe)[1:10]
head(E_GEOD_53987_Schizophrenia_Temporal_Lobe)[1:10]
head(E_GEOD_35978_Schizophrenia_Cerebellum)[1:10]
head(E_GEOD_35978_Schizophrenia_Parietal_Lobe)[1:10]

# PARKINONS DISEASE
head(E_GEOD_20168_PD_Frontal_Lobe)[1:10]

# HUNTINDONS DISEASE
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133A)[1:10]
head(E_GEOD_3790_HD_Cerebellum_Affy_U133B)[1:10]

# MAJOR DEPRESSIVE DISORDER
head(E_GEOD_53987_MDD_Frontal_Lobe)[1:10]
head(E_GEOD_53987_MDD_Temporal_Lobe)[1:10]

##### METAQC - CREATE FORMAT AND WRITE #####

# convert data to "MetaOmics" format - keep only diagnosis- create function for conversion of data

convert_data_to_MetaOmics_format<-function(dataset){
  # remove unwanted pheno - keep only diagnosis + expression (merge by rownames)
  dataset<-merge(dataset[1], dataset[7:dim(dataset)[2]], by="row.names")
  rownames(dataset)<-dataset$Row.names
  dataset$Row.names<-NULL
  # recode diagnosis column - case=1, control=0
  dataset[dataset$Diagnosis=="control",1]<-"0"
  dataset[dataset$Diagnosis=="case",1]<-"1"
  # re-label Diagnosis to Label
  colnames(dataset)[1]<-"label"
  # change dataframe to numeric and transpose
  dataset<-as.data.frame(t(data.matrix(dataset)))
  print(head(dataset)[1:5])
  return(dataset)
}

# apply function to all datasets - leave out occipital lobe

# FRONTAL LOBE
E_GEOD_12649_BD_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_12649_BD_Frontal_Lobe)
E_GEOD_12649_Schizophrenia_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_12649_Schizophrenia_Frontal_Lobe)
E_GEOD_17612_Schizophrenia_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_17612_Schizophrenia_Frontal_Lobe)
E_GEOD_20168_PD_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_20168_PD_Frontal_Lobe)
E_GEOD_21138_Schizophrenia_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_21138_Schizophrenia_Frontal_Lobe)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B)
E_GEOD_5388_BD_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_5388_BD_Frontal_Lobe)
E_GEOD_53987_BD_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_BD_Frontal_Lobe)
E_GEOD_53987_MDD_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_MDD_Frontal_Lobe)
E_GEOD_53987_Schizophrenia_Frontal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_Schizophrenia_Frontal_Lobe)

# TEMPORAL LOBE
E_GEOD_21935_Schizophrenia_Temporal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_21935_Schizophrenia_Temporal_Lobe)
E_GEOD_53987_BD_Temporal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_BD_Temporal_Lobe)
E_GEOD_53987_MDD_Temporal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_MDD_Temporal_Lobe)
E_GEOD_53987_Schizophrenia_Temporal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_53987_Schizophrenia_Temporal_Lobe)

# CEREBELLUM
E_GEOD_35978_BD_Cerebellum_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_35978_BD_Cerebellum)
E_GEOD_35978_Schizophrenia_Cerebellum_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_35978_Schizophrenia_Cerebellum)
E_GEOD_3790_HD_Cerebellum_Affy_U133A_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_3790_HD_Cerebellum_Affy_U133A)
E_GEOD_3790_HD_Cerebellum_Affy_U133B_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_3790_HD_Cerebellum_Affy_U133B)

# PARIETAL LOBE
E_GEOD_35978_BD_Parietal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_35978_BD_Parietal_Lobe)
E_GEOD_35978_Schizophrenia_Parietal_Lobe_MetaOmics_format <- convert_data_to_MetaOmics_format(E_GEOD_35978_Schizophrenia_Parietal_Lobe)

# write out dataframe - ***NEED To CONVERT TO THIS FORMAT IN R

setwd(Data_in_metaomics_format)

# write out data

# FRONTAL LOBE
write.table(E_GEOD_12649_BD_Frontal_Lobe_MetaOmics_format, "E_GEOD_12649_BD_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_12649_Schizophrenia_Frontal_Lobe_MetaOmics_format, "E_GEOD_12649_Schizophrenia_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_17612_Schizophrenia_Frontal_Lobe_MetaOmics_format, "E_GEOD_17612_Schizophrenia_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_20168_PD_Frontal_Lobe_MetaOmics_format, "E_GEOD_20168_PD_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_21138_Schizophrenia_Frontal_Lobe_MetaOmics_format, "E_GEOD_21138_Schizophrenia_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_MetaOmics_format, "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_MetaOmics_format, "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_5388_BD_Frontal_Lobe_MetaOmics_format, "E_GEOD_5388_BD_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_BD_Frontal_Lobe_MetaOmics_format, "E_GEOD_53987_BD_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_MDD_Frontal_Lobe_MetaOmics_format, "E_GEOD_53987_MDD_Frontal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_Schizophrenia_Frontal_Lobe_MetaOmics_format, "E_GEOD_53987_Schizophrenia_Frontal_Lobe_MetaOmics_format.txt", sep="\t")

# TEMPORAL LOBE
write.table(E_GEOD_21935_Schizophrenia_Temporal_Lobe_MetaOmics_format, "E_GEOD_21935_Schizophrenia_Temporal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_BD_Temporal_Lobe_MetaOmics_format, "E_GEOD_53987_BD_Temporal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_MDD_Temporal_Lobe_MetaOmics_format, "E_GEOD_53987_MDD_Temporal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_53987_Schizophrenia_Temporal_Lobe_MetaOmics_format, "E_GEOD_53987_Schizophrenia_Temporal_Lobe_MetaOmics_format.txt", sep="\t")

# CEREBELLUM
write.table(E_GEOD_35978_BD_Cerebellum_MetaOmics_format, "E_GEOD_35978_BD_Cerebellum_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_35978_Schizophrenia_Cerebellum_MetaOmics_format, "E_GEOD_35978_Schizophrenia_Cerebellum_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_3790_HD_Cerebellum_Affy_U133A_MetaOmics_format, "E_GEOD_3790_HD_Cerebellum_Affy_U133A_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_3790_HD_Cerebellum_Affy_U133B_MetaOmics_format, "E_GEOD_3790_HD_Cerebellum_Affy_U133B_MetaOmics_format.txt", sep="\t")

# PARIETAL LOBE
write.table(E_GEOD_35978_BD_Parietal_Lobe_MetaOmics_format, "E_GEOD_35978_BD_Parietal_Lobe_MetaOmics_format.txt", sep="\t")
write.table(E_GEOD_35978_Schizophrenia_Parietal_Lobe_MetaOmics_format, "E_GEOD_35978_Schizophrenia_Parietal_Lobe_MetaOmics_format.txt", sep="\t")

##### METAQC - READ IN CREATED FORMAT #####

setwd(Data_in_metaomics_format)

# read data by brain region

# FRONTAL LOBE - 11 datasets

Frontal_Lobe<-MetaDE.Read(c("E_GEOD_12649_BD_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_12649_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_17612_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_20168_PD_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_21138_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_MetaOmics_format" ,
                            "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_MetaOmics_format" ,
                            "E_GEOD_5388_BD_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_53987_BD_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_53987_MDD_Frontal_Lobe_MetaOmics_format" ,
                            "E_GEOD_53987_Schizophrenia_Frontal_Lobe_MetaOmics_format"), 
                          via="txt", # txt format
                          skip=rep(1,11), # skip 1 row - 11 datasets
                          log=T, # data is log transformed
                          matched=T) #using entrez gene ID instead of gene symbols

# TEMPORAL LOBE - 4 datasets

Temporal_Lobe<-MetaDE.Read(c("E_GEOD_21935_Schizophrenia_Temporal_Lobe_MetaOmics_format" ,
                             "E_GEOD_53987_BD_Temporal_Lobe_MetaOmics_format" ,
                             "E_GEOD_53987_MDD_Temporal_Lobe_MetaOmics_format" ,
                             "E_GEOD_53987_Schizophrenia_Temporal_Lobe_MetaOmics_format"),
                           via="txt", # txt format
                           skip=rep(1,4), # skip 1 row - 4 datasets
                           log=T, # data is log transformed
                           matched=T) #using entrez gene ID instead of gene symbols


# CEREBELLUM - 4 datasets

Cerebellum<-MetaDE.Read(c("E_GEOD_35978_BD_Cerebellum_MetaOmics_format" ,
                          "E_GEOD_35978_Schizophrenia_Cerebellum_MetaOmics_format" ,
                          "E_GEOD_3790_HD_Cerebellum_Affy_U133A_MetaOmics_format" ,
                          "E_GEOD_3790_HD_Cerebellum_Affy_U133B_MetaOmics_format" ),
                        via="txt", # txt format
                        skip=rep(1,4), # skip 1 row - 4 datasets
                        log=T, # data is log transformed
                        matched=T) #using entrez gene ID instead of gene symbols

# PARIETAL LOBE - 2 datasets

Parietal_Lobe<-MetaDE.Read(c("E_GEOD_35978_BD_Parietal_Lobe_MetaOmics_format" ,
                             "E_GEOD_35978_Schizophrenia_Parietal_Lobe_MetaOmics_format"), 
                           via="txt", # txt format
                           skip=rep(1,2), # skip 1 row - 2 datasets
                           log=T, # data is log transformed
                           matched=T) #using entrez gene ID instead of gene symbols

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
for(i in 1:11){ 
  colnames(Frontal_Lobe_merged_Filtered[[i]][[1]])<-Frontal_Lobe_merged_Filtered[[i]][[2]]
  Frontal_Lobe_QC[[i]]<-impute.knn(Frontal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Frontal_Lobe_QC)<-names(Frontal_Lobe_merged_Filtered)

# moves diagnosis label to column name - Parietal_Lobe_merged_Filtered - impute missing - 
for(i in 1:2){ 
  colnames(Parietal_Lobe_merged_Filtered[[i]][[1]])<-Parietal_Lobe_merged_Filtered[[i]][[2]]
  Parietal_Lobe_QC[[i]]<-impute.knn(Parietal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Parietal_Lobe_QC)<-names(Parietal_Lobe_merged_Filtered)

# moves diagnosis label to column name - Temporal_Lobe_merged_Filtered - impute missing - 
for(i in 1:4){ 
  colnames(Temporal_Lobe_merged_Filtered[[i]][[1]])<-Temporal_Lobe_merged_Filtered[[i]][[2]]
  Temporal_Lobe_QC[[i]]<-impute.knn(Temporal_Lobe_merged_Filtered[[i]][[1]])$data
}

names(Temporal_Lobe_QC)<-names(Temporal_Lobe_merged_Filtered)

# moves diagnosis label to column name - Cerebellum_merged_Filtered - impute missing - 
for(i in 1:4){ 
  colnames(Cerebellum_merged_Filtered[[i]][[1]])<-Cerebellum_merged_Filtered[[i]][[2]]
  Cerebellum_QC[[i]]<-impute.knn(Cerebellum_merged_Filtered[[i]][[1]])$data
}

names(Cerebellum_QC)<-names(Cerebellum_merged_Filtered)

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

# runQC- **script freezes when 1st downloading fileForCQCp. downloaded manually and inserted into working folder. odwnloaded from: http://software.broadinstitute.org/gsea/downloads.jsp

run_meta_qc2<-function(dataset){
  runQC(dataset, 
        nPath=NULL, # 
        B=1e4, # The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
        pvalCut=.05, #P-value threshold used for AQC calculation.
        #pvalAdjust=T, #whether to apply p-value adjustment due to multiple testing (B-H procedure is used).
        fileForCQCp="c2.all.v5.1.entrez.gmt") # latest version of c2.all.v5.1.symbols.gmt downloaded 02/09/2016
}
 
run_meta_qc2(Frontal_Lobe_QC_check)
#run_meta_qc2(Parietal_Lobe_QC_check) # - "Error in { : task 1 failed - "not enough 'y' observations"
run_meta_qc2(Temporal_Lobe_QC_check)
run_meta_qc2(Cerebellum_QC_check) # - "Error in { : task 1 failed - "not enough 'y' observations"

Frontal_Lobe_QC_check
#Parietal_Lobe_QC_check
Temporal_Lobe_QC_check
Cerebellum_QC_check

#plot PCA

plot(Frontal_Lobe_QC_check)
#plot(Parietal_Lobe_QC_check)
plot(Temporal_Lobe_QC_check)
plot(Cerebellum_QC_check)

##### PLOT QC PLOTS TO PDF #####

setwd(MetaQC_plots_dir)

pdf("Frontal_lobe_MetaQC.pdf")
plot(Frontal_Lobe_QC_check)
mtext("Frontal Lobe Datasets", cex=2, line = 1.5)
dev.off()

pdf("Cerebellum_MetaQC.pdf")
plot(Cerebellum_QC_check)
mtext("Cerebellum Datasets", cex=2, line = 1.5)
dev.off()

pdf("Temporal_lobe_MetaQC.pdf")
plot(Temporal_Lobe_QC_check)
mtext("Temporal Lobe Datasets", cex=2, line = 1.5)
dev.off()

###############################################
#####                                     #####
##### ANALYSIS 1 - GROUP BY BRAIN REGION  #####
#####                                     #####
###############################################

##### RE-CREATE QC OBJECT - WITHOUT REMOVING GENES #####

#setwd to work dir
setwd(work_dir)

create_QC_object<-function(metaQC_object, MVperc, gene_filter){ 
  # merge by at least common probes, 0=100% common, 0.8=20% common
  raw_data_merged<-MetaDE.merge(metaQC_object, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_merged[[1]][[1]])[1]), sep=" ")
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_merged_Filtered<-MetaDE.filter(raw_data_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_merged_Filtered[[1]][[1]])[1]), sep=" ")
  return(raw_data_merged_Filtered)
}


Frontal_Lobe_QC<-create_QC_object(metaQC_object=Frontal_Lobe,
                                  MVperc=0.8,
                                  gene_filter=c(0,0))

Temporal_Lobe_QC<-create_QC_object(metaQC_object=Temporal_Lobe,
                                   MVperc=0.8,
                                   gene_filter=c(0,0))

Cerebellum_QC<-create_QC_object(metaQC_object=Cerebellum,
                                MVperc=0.8,
                                gene_filter=c(0,0))

Parietal_Lobe_QC<-create_QC_object(metaQC_object=Parietal_Lobe,
                                   MVperc=0.8,
                                   gene_filter=c(0,0))

##### META DE ANALYSIS "AW.OC"#####

# Frontal Lobe
Frontal_lobe_DE<-MetaDE.rawdata(Frontal_Lobe_QC,
                                ind.method = rep("modt", length(names(Frontal_Lobe_QC))),# moderate t-test - same as limma
                                meta.method="AW.OC",
                                nperm=100,
                                miss.tol=0,
                                ind.tail="high")

tail(Frontal_lobe_DE$meta.analysis$stat)
tail(Frontal_lobe_DE$meta.analysis$FDR)
tail(Frontal_lobe_DE$meta.analysis$AW.weight)
tail(Frontal_lobe_DE$meta.analysis$pval)

#Parietal Lobe
Parietal_lobe_DE<-MetaDE.rawdata(Parietal_Lobe_QC,
                                 ind.method = rep("modt", length(names(Parietal_Lobe_QC))),# moderate t-test - same as limma
                                 meta.method="AW.OC",
                                 nperm=100,
                                 miss.tol=0,
                                 ind.tail="high")

tail(Parietal_lobe_DE$meta.analysis$stat)
tail(Parietal_lobe_DE$meta.analysis$FDR)
tail(Parietal_lobe_DE$meta.analysis$AW.weight)
tail(Parietal_lobe_DE$meta.analysis$pval)

#Cerebellum
Cerebellum_DE<-MetaDE.rawdata(Cerebellum_QC,
                              ind.method = rep("modt", length(names(Cerebellum_QC))),# moderate t-test - same as limma
                              meta.method="AW.OC",
                              nperm=100,
                              miss.tol=0,
                              ind.tail="high")

tail(Cerebellum_DE$meta.analysis$stat)
tail(Cerebellum_DE$meta.analysis$FDR)
tail(Cerebellum_DE$meta.analysis$AW.weight)
tail(Cerebellum_DE$meta.analysis$pval)

#temporal Lobe
Temporal_lobe_DE<-MetaDE.rawdata(Temporal_Lobe_QC,
                                    ind.method = rep("modt", length(names(Temporal_Lobe_QC))),# moderate t-test - same as limma
                                    meta.method="AW.OC",
                                    nperm=100,
                                    miss.tol=0,
                                    ind.tail="high")


tail(Temporal_lobe_DE$meta.analysis$stat)
tail(Temporal_lobe_DE$meta.analysis$FDR)
tail(Temporal_lobe_DE$meta.analysis$AW.weight)
tail(Temporal_lobe_DE$meta.analysis$pval)

##### HEATMAP "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC" #####

setwd(DE_Heamaps_dir)

pdf("Frontal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Frontal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe AW.OC DEG ", line = 1.5)
dev.off()

#parietal lobe
pdf("Parietal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Parietal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe AW.OC DEG ", line = 1.5)
dev.off()

# cerebellum
pdf("Cerebellum_DE_heatmaps.pdf")
heatmap.sig.genes(Cerebellum_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum AW.OC DEG ", line = 1.5)
dev.off()

#Temporal Lobe
pdf("Temporal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Temporal_lobe_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe AW.OC DEG ", line = 1.5)
dev.off()

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

# read limma DE table for logFC values

E_GEOD_5388_BD_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-5388/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_12649_BD_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-12649/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_35978_BD_Cerebellum_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-35978/Cerebellum/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_35978_BD_Parietal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-35978/Parietal_Lobe/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_53987_BD_Temporal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-53987/Hippocampus/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_53987_BD_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Bipolar_Disorder/E-GEOD-53987/Prefrontal_cortex/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)

E_GEOD_3790_HD_Cerebellum_Affy_U133A_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Huntingdons_Disease/E-GEOD-3790/Affy_U133A/Cerebellum/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Huntingdons_Disease/E-GEOD-3790/Affy_U133A/Frontal_Lobe/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Huntingdons_Disease/E-GEOD-3790/Affy_U133B/Frontal_Lobe/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_3790_HD_Cerebellum_Affy_U133B_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Huntingdons_Disease/E-GEOD-3790/Affy_U133B/Cerebellum/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)

E_GEOD_53987_MDD_Temporal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Major_Depressive_Disorder/E-GEOD-53987/Hippocampus/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_53987_MDD_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Major_Depressive_Disorder/E-GEOD-53987/Prefrontal_Cortex/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)

E_GEOD_20168_PD_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Parkinsons/E-GEOD-20168/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)

E_GEOD_12649_Schizophrenia_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-12649/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_17612_Schizophrenia_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-17612/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_21138_Schizophrenia_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-21138/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_21935_Schizophrenia_Temporal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-21935/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_35978_Schizophrenia_Cerebellum_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-35978/Cerebellum/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_35978_Schizophrenia_Parietal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E-GEOD-35978/Parietal_Lobe/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_53987_Schizophrenia_Temporal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E_GEOD-53987/Hippocampus/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)
E_GEOD_53987_Schizophrenia_Frontal_Lobe_DE <- read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Schizophrenia/E_GEOD-53987/Prefrontal_Cortex/DE_results/Full_Differential_Expression_Results.txt", as.is=T, head=T)

# change colnames to include dataset name

colnames(E_GEOD_12649_BD_Frontal_Lobe_DE)<-c("E_GEOD_12649_BD_Frontal_Lobe_logFC", "E_GEOD_12649_BD_Frontal_Lobe_CI.L", "E_GEOD_12649_BD_Frontal_Lobe_CI.R", "E_GEOD_12649_BD_Frontal_Lobe_AveExpr", "E_GEOD_12649_BD_Frontal_Lobe_t", "E_GEOD_12649_BD_Frontal_Lobe_P.Value",  "E_GEOD_12649_BD_Frontal_Lobe_adj.P.Value", "E_GEOD_12649_BD_Frontal_Lobe_B")
colnames(E_GEOD_12649_Schizophrenia_Frontal_Lobe_DE)<-c("E_GEOD_12649_Schizophrenia_Frontal_Lobe_logFC", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_CI.L", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_CI.R", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_AveExpr", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_t", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_P.Value",  "E_GEOD_12649_Schizophrenia_Frontal_Lobe_adj.P.Value", "E_GEOD_12649_Schizophrenia_Frontal_Lobe_B")
colnames(E_GEOD_17612_Schizophrenia_Frontal_Lobe_DE)<-c("E_GEOD_17612_Schizophrenia_Frontal_Lobe_logFC", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_CI.L", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_CI.R", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_AveExpr", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_t", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_P.Value",  "E_GEOD_17612_Schizophrenia_Frontal_Lobe_adj.P.Value", "E_GEOD_17612_Schizophrenia_Frontal_Lobe_B")
colnames(E_GEOD_20168_PD_Frontal_Lobe_DE)<-c("E_GEOD_20168_PD_Frontal_Lobe_logFC", "E_GEOD_20168_PD_Frontal_Lobe_CI.L", "E_GEOD_20168_PD_Frontal_Lobe_CI.R", "E_GEOD_20168_PD_Frontal_Lobe_AveExpr", "E_GEOD_20168_PD_Frontal_Lobe_t", "E_GEOD_20168_PD_Frontal_Lobe_P.Value",  "E_GEOD_20168_PD_Frontal_Lobe_adj.P.Value", "E_GEOD_20168_PD_Frontal_Lobe_B")
colnames(E_GEOD_21138_Schizophrenia_Frontal_Lobe_DE)<-c("E_GEOD_21138_Schizophrenia_Frontal_Lobe_logFC", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_CI.L", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_CI.R", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_AveExpr", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_t", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_P.Value",  "E_GEOD_21138_Schizophrenia_Frontal_Lobe_adj.P.Value", "E_GEOD_21138_Schizophrenia_Frontal_Lobe_B")
colnames(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_DE)<-c("E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_logFC", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_CI.L", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_CI.R", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_AveExpr", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_t", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_P.Value",  "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_adj.P.Value", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_B")
colnames(E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_DE)<-c("E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_logFC", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_CI.L", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_CI.R", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_AveExpr", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_t", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_P.Value",  "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_adj.P.Value", "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_B")
colnames(E_GEOD_5388_BD_Frontal_Lobe_DE)<-c("E_GEOD_5388_BD_Frontal_Lobe_logFC", "E_GEOD_5388_BD_Frontal_Lobe_CI.L", "E_GEOD_5388_BD_Frontal_Lobe_CI.R", "E_GEOD_5388_BD_Frontal_Lobe_AveExpr", "E_GEOD_5388_BD_Frontal_Lobe_t", "E_GEOD_5388_BD_Frontal_Lobe_P.Value",  "E_GEOD_5388_BD_Frontal_Lobe_adj.P.Value", "E_GEOD_5388_BD_Frontal_Lobe_B")
colnames(E_GEOD_53987_BD_Frontal_Lobe_DE)<-c("E_GEOD_53987_BD_Frontal_Lobe_logFC", "E_GEOD_53987_BD_Frontal_Lobe_CI.L", "E_GEOD_53987_BD_Frontal_Lobe_CI.R", "E_GEOD_53987_BD_Frontal_Lobe_AveExpr", "E_GEOD_53987_BD_Frontal_Lobe_t", "E_GEOD_53987_BD_Frontal_Lobe_P.Value",  "E_GEOD_53987_BD_Frontal_Lobe_adj.P.Value", "E_GEOD_53987_BD_Frontal_Lobe_B")
colnames(E_GEOD_53987_MDD_Frontal_Lobe_DE)<-c("E_GEOD_53987_MDD_Frontal_Lobe_logFC", "E_GEOD_53987_MDD_Frontal_Lobe_CI.L", "E_GEOD_53987_MDD_Frontal_Lobe_CI.R", "E_GEOD_53987_MDD_Frontal_Lobe_AveExpr", "E_GEOD_53987_MDD_Frontal_Lobe_t", "E_GEOD_53987_MDD_Frontal_Lobe_P.Value",  "E_GEOD_53987_MDD_Frontal_Lobe_adj.P.Value", "E_GEOD_53987_MDD_Frontal_Lobe_B")
colnames(E_GEOD_53987_Schizophrenia_Frontal_Lobe_DE)<-c("E_GEOD_53987_Schizophrenia_Frontal_Lobe_logFC", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_CI.L", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_CI.R", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_AveExpr", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_t", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_P.Value",  "E_GEOD_53987_Schizophrenia_Frontal_Lobe_adj.P.Value", "E_GEOD_53987_Schizophrenia_Frontal_Lobe_B")

colnames(E_GEOD_21935_Schizophrenia_Temporal_Lobe_DE)<-c("E_GEOD_21935_Schizophrenia_Temporal_Lobe_logFC", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_CI.L", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_CI.R", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_AveExpr", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_t", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_P.Value",  "E_GEOD_21935_Schizophrenia_Temporal_Lobe_adj.P.Value", "E_GEOD_21935_Schizophrenia_Temporal_Lobe_B")
colnames(E_GEOD_53987_BD_Temporal_Lobe_DE)<-c("E_GEOD_53987_BD_Temporal_Lobe_logFC", "E_GEOD_53987_BD_Temporal_Lobe_CI.L", "E_GEOD_53987_BD_Temporal_Lobe_CI.R", "E_GEOD_53987_BD_Temporal_Lobe_AveExpr", "E_GEOD_53987_BD_Temporal_Lobe_t", "E_GEOD_53987_BD_Temporal_Lobe_P.Value",  "E_GEOD_53987_BD_Temporal_Lobe_adj.P.Value", "E_GEOD_53987_BD_Temporal_Lobe_B")
colnames(E_GEOD_53987_MDD_Temporal_Lobe_DE)<-c("E_GEOD_53987_MDD_Temporal_Lobe_logFC", "E_GEOD_53987_MDD_Temporal_Lobe_CI.L", "E_GEOD_53987_MDD_Temporal_Lobe_CI.R", "E_GEOD_53987_MDD_Temporal_Lobe_AveExpr", "E_GEOD_53987_MDD_Temporal_Lobe_t", "E_GEOD_53987_MDD_Temporal_Lobe_P.Value",  "E_GEOD_53987_MDD_Temporal_Lobe_adj.P.Value", "E_GEOD_53987_MDD_Temporal_Lobe_B")
colnames(E_GEOD_53987_Schizophrenia_Temporal_Lobe_DE)<-c("E_GEOD_53987_Schizophrenia_Temporal_Lobe_logFC", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_CI.L", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_CI.R", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_AveExpr", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_t", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_P.Value",  "E_GEOD_53987_Schizophrenia_Temporal_Lobe_adj.P.Value", "E_GEOD_53987_Schizophrenia_Temporal_Lobe_B")

colnames(E_GEOD_35978_BD_Cerebellum_DE)<-c("E_GEOD_35978_BD_Cerebellum_logFC", "E_GEOD_35978_BD_Cerebellum_CI.L", "E_GEOD_35978_BD_Cerebellum_CI.R", "E_GEOD_35978_BD_Cerebellum_AveExpr", "E_GEOD_35978_BD_Cerebellum_t", "E_GEOD_35978_BD_Cerebellum_P.Value",  "E_GEOD_35978_BD_Cerebellum_adj.P.Value", "E_GEOD_35978_BD_Cerebellum_B")
colnames(E_GEOD_35978_Schizophrenia_Cerebellum_DE)<-c("E_GEOD_35978_Schizophrenia_Cerebellum_logFC", "E_GEOD_35978_Schizophrenia_Cerebellum_CI.L", "E_GEOD_35978_Schizophrenia_Cerebellum_CI.R", "E_GEOD_35978_Schizophrenia_Cerebellum_AveExpr", "E_GEOD_35978_Schizophrenia_Cerebellum_t", "E_GEOD_35978_Schizophrenia_Cerebellum_P.Value",  "E_GEOD_35978_Schizophrenia_Cerebellum_adj.P.Value", "E_GEOD_35978_Schizophrenia_Cerebellum_B")
colnames(E_GEOD_3790_HD_Cerebellum_Affy_U133A_DE)<-c("E_GEOD_3790_HD_Cerebellum_Affy_U133A_logFC", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_CI.L", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_CI.R", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_AveExpr", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_t", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_P.Value",  "E_GEOD_3790_HD_Cerebellum_Affy_U133A_adj.P.Value", "E_GEOD_3790_HD_Cerebellum_Affy_U133A_B")
colnames(E_GEOD_3790_HD_Cerebellum_Affy_U133B_DE)<-c("E_GEOD_3790_HD_Cerebellum_Affy_U133B_logFC", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_CI.L", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_CI.R", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_AveExpr", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_t", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_P.Value",  "E_GEOD_3790_HD_Cerebellum_Affy_U133B_adj.P.Value", "E_GEOD_3790_HD_Cerebellum_Affy_U133B_B")

colnames(E_GEOD_35978_BD_Parietal_Lobe_DE)<-c("E_GEOD_35978_BD_Parietal_Lobe_logFC", "E_GEOD_35978_BD_Parietal_Lobe_CI.L", "E_GEOD_35978_BD_Parietal_Lobe_CI.R", "E_GEOD_35978_BD_Parietal_Lobe_AveExpr", "E_GEOD_35978_BD_Parietal_Lobe_t", "E_GEOD_35978_BD_Parietal_Lobe_P.Value",  "E_GEOD_35978_BD_Parietal_Lobe_adj.P.Value", "E_GEOD_35978_BD_Parietal_Lobe_B")
colnames(E_GEOD_35978_Schizophrenia_Parietal_Lobe_DE)<-c("E_GEOD_35978_Schizophrenia_Parietal_Lobe_logFC", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_CI.L", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_CI.R", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_AveExpr", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_t", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_P.Value",  "E_GEOD_35978_Schizophrenia_Parietal_Lobe_adj.P.Value", "E_GEOD_35978_Schizophrenia_Parietal_Lobe_B")

# check
head(E_GEOD_35978_Schizophrenia_Parietal_Lobe_DE)[1:8]

# merge function

MyMerge       <- function(x, y){
  df<- merge(x, y, by= "row.names", all=T)
  rownames(df)<- df$Row.names
  df$Row.names<- NULL
  return(df)
}

limma_results<-Reduce(MyMerge, list(E_GEOD_12649_BD_Frontal_Lobe_DE ,
                                        E_GEOD_12649_Schizophrenia_Frontal_Lobe_DE ,
                                        E_GEOD_17612_Schizophrenia_Frontal_Lobe_DE ,
                                        E_GEOD_20168_PD_Frontal_Lobe_DE ,
                                        E_GEOD_21138_Schizophrenia_Frontal_Lobe_DE ,
                                        E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_DE ,
                                        E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_DE ,
                                        E_GEOD_5388_BD_Frontal_Lobe_DE ,
                                        E_GEOD_53987_BD_Frontal_Lobe_DE ,
                                        E_GEOD_53987_MDD_Frontal_Lobe_DE ,
                                        E_GEOD_53987_Schizophrenia_Frontal_Lobe_DE ,
                                        E_GEOD_21935_Schizophrenia_Temporal_Lobe_DE ,
                                        E_GEOD_53987_BD_Temporal_Lobe_DE ,
                                        E_GEOD_53987_MDD_Temporal_Lobe_DE ,
                                        E_GEOD_53987_Schizophrenia_Temporal_Lobe_DE ,
                                        E_GEOD_35978_BD_Cerebellum_DE ,
                                        E_GEOD_35978_Schizophrenia_Cerebellum_DE ,
                                        E_GEOD_3790_HD_Cerebellum_Affy_U133A_DE ,
                                        E_GEOD_3790_HD_Cerebellum_Affy_U133B_DE ,
                                        E_GEOD_35978_BD_Parietal_Lobe_DE ,
                                        E_GEOD_35978_Schizophrenia_Parietal_Lobe_DE))
                                        
dim(limma_results)
head(limma_results)[100:110]

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


setwd(AW_results_dir)

## FRONTAL LOBE

#extract AW results and outcome
Frontal_lobe_AW <- Frontal_lobe_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Frontal_lobe_DEG_regulation_good_genes))]
#move rownames to entrez id column
Frontal_lobe_AW$Entrez_Gene_ID<-rownames(Frontal_lobe_AW)
#rearrange columns
Frontal_lobe_AW<-Frontal_lobe_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Frontal_lobe_AW)[3]<-"Gene_Regulation"
#order by p val
Frontal_lobe_AW<-Frontal_lobe_AW[order(Frontal_lobe_AW$AW.OC),]
#check
head(Frontal_lobe_AW)
nrow(Frontal_lobe_AW)

#write
write.table(Frontal_lobe_AW, file="Frontal_lobe_AW.OC", row.names=F, sep="\t", quote=F)

## PARIETAL LOBE

#extract AW results and outcome
Parietal_lobe_AW <- Parietal_lobe_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Parietal_lobe_DEG_regulation_good_genes))]
#move rownames to entrez id column
Parietal_lobe_AW$Entrez_Gene_ID<-rownames(Parietal_lobe_AW)
#rearrange columns
Parietal_lobe_AW<-Parietal_lobe_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Parietal_lobe_AW)[3]<-"Gene_Regulation"
#order by p val
Parietal_lobe_AW<-Parietal_lobe_AW[order(Parietal_lobe_AW$AW.OC),]
#check
head(Parietal_lobe_AW)
nrow(Parietal_lobe_AW)

#write
write.table(Parietal_lobe_AW, file="Parietal_lobe_AW.OC", row.names=F, sep="\t", quote=F)

## TEMPORAL LOBE

#extract AW results and outcome
Temporal_lobe_AW <- Temporal_lobe_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Temporal_lobe_DEG_regulation_good_genes))]
#move rownames to entrez id column
Temporal_lobe_AW$Entrez_Gene_ID<-rownames(Temporal_lobe_AW)
#rearrange columns
Temporal_lobe_AW<-Temporal_lobe_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Temporal_lobe_AW)[3]<-"Gene_Regulation"
#order by p val
Temporal_lobe_AW<-Temporal_lobe_AW[order(Temporal_lobe_AW$AW.OC),]
#check
head(Temporal_lobe_AW)
nrow(Temporal_lobe_AW)

#write
write.table(Temporal_lobe_AW, file="Temporal_lobe_AW.OC", row.names=F, sep="\t", quote=F)

## CEREBELLUM

#extract AW results and outcome
Cerebellum_AW <- Cerebellum_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Cerebellum_DEG_regulation_good_genes))]
#move rownames to entrez id column
Cerebellum_AW$Entrez_Gene_ID<-rownames(Cerebellum_AW)
#rearrange columns
Cerebellum_AW<-Cerebellum_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Cerebellum_AW)[3]<-"Gene_Regulation"
#order by p val
Cerebellum_AW<-Cerebellum_AW[order(Cerebellum_AW$AW.OC),]
#check
head(Cerebellum_AW)
nrow(Cerebellum_AW)

#write
write.table(Cerebellum_AW, file="Cerebellum_AW.OC", row.names=F, sep="\t", quote=F)

##### #########################################
#####                                     #####
#####    ANALYSIS 2 - GROUP BY DISEASE    #####
#####                                     #####
###############################################

##### METAQC - READ IN CREATED FORMAT #####

setwd(Data_in_metaomics_format)

# read data by brain region

# BIPOLAR DISORDER - 6 datasets

Bipolar_Disorder<-MetaDE.Read(c("E_GEOD_12649_BD_Frontal_Lobe_MetaOmics_format" ,
                                "E_GEOD_35978_BD_Parietal_Lobe_MetaOmics_format" ,
                                "E_GEOD_53987_BD_Temporal_Lobe_MetaOmics_format" ,
                                "E_GEOD_35978_BD_Cerebellum_MetaOmics_format" ,
                                "E_GEOD_5388_BD_Frontal_Lobe_MetaOmics_format" ,
                                "E_GEOD_53987_BD_Frontal_Lobe_MetaOmics_format"), 
                              via="txt", # txt format
                              skip=rep(1,6), # skip 1 row - 6 datasets
                              log=T, # data is log transformed
                              matched=T) #using entrez gene ID instead of gene symbols

# HUNTINGDON'S DISEASE - 4 datasets

Huntingdons_Disease<-MetaDE.Read(c("E_GEOD_3790_HD_Cerebellum_Affy_U133A_MetaOmics_format" ,
                                   "E_GEOD_3790_HD_Cerebellum_Affy_U133B_MetaOmics_format" ,
                                   "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133A_MetaOmics_format" ,
                                   "E_GEOD_3790_HD_Frontal_Lobe_Affy_U133B_MetaOmics_format") ,
                                 via="txt", # txt format
                                 skip=rep(1,4), # skip 1 row - 4 datasets
                                 log=T, # data is log transformed
                                 matched=T) #using entrez gene ID instead of gene symbols


# MAJOR DEPRESSIVE DISORDER - 2 datasets

Major_Depressive_Disorder<-MetaDE.Read(c("E_GEOD_53987_MDD_Temporal_Lobe_MetaOmics_format" ,
                                         "E_GEOD_53987_MDD_Frontal_Lobe_MetaOmics_format" ),
                                       via="txt", # txt format
                                       skip=rep(1,2), # skip 1 row - 2 datasets
                                       log=T, # data is log transformed
                                       matched=T) #using entrez gene ID instead of gene symbols

# SCHIZOPHRENIA - 8 datasets

Schizophrenia<-MetaDE.Read(c("E_GEOD_21935_Schizophrenia_Temporal_Lobe_MetaOmics_format" ,
                             "E_GEOD_12649_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                             "E_GEOD_17612_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                             "E_GEOD_21138_Schizophrenia_Frontal_Lobe_MetaOmics_format" ,
                             "E_GEOD_35978_Schizophrenia_Parietal_Lobe_MetaOmics_format",
                             "E_GEOD_53987_Schizophrenia_Temporal_Lobe_MetaOmics_format" ,
                             "E_GEOD_35978_Schizophrenia_Cerebellum_MetaOmics_format" ,
                             "E_GEOD_53987_Schizophrenia_Frontal_Lobe_MetaOmics_format"), 
                           via="txt", # txt format
                           skip=rep(1,8), # skip 1 row - 8 datasets
                           log=T, # data is log transformed
                           matched=T) #using entrez gene ID instead of gene symbols

setwd(work_dir)

##### RE-CREATE QC OBJECT - WITHOUT REMOVING GENES #####

#setwd to work dir
setwd(work_dir)

create_QC_object<-function(metaQC_object, MVperc, gene_filter){ 
  # merge by at least common probes, 0=100% common, 0.8=20% common
  raw_data_merged<-MetaDE.merge(metaQC_object, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_merged[[1]][[1]])[1]), sep=" ")
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_merged_Filtered<-MetaDE.filter(raw_data_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_merged_Filtered[[1]][[1]])[1]), sep=" ")
  return(raw_data_merged_Filtered)
}


Bipolar_Disorder_QC<-create_QC_object(metaQC_object=Bipolar_Disorder,
                                      MVperc=0.8,
                                      gene_filter=c(0,0))

Huntingdons_Disease_QC<-create_QC_object(metaQC_object=Huntingdons_Disease,
                                         MVperc=0.8,
                                         gene_filter=c(0,0))

Major_Depressive_Disorder_QC<-create_QC_object(metaQC_object=Major_Depressive_Disorder,
                                               MVperc=0.8,
                                               gene_filter=c(0,0))

Schizophrenia_QC<-create_QC_object(metaQC_object=Schizophrenia,
                                   MVperc=0.8,
                                   gene_filter=c(0,0))

##### META DE ANALYSIS "AW.OC"#####

# Bipolar_Disorder
Bipolar_Disorder_DE<-MetaDE.rawdata(Bipolar_Disorder_QC,
                                    ind.method = rep("modt", length(names(Bipolar_Disorder_QC))),# moderate t-test - same as limma
                                    meta.method="AW.OC",
                                    nperm=100,
                                    miss.tol=0,
                                    ind.tail="high")

tail(Bipolar_Disorder_DE$meta.analysis$stat)
tail(Bipolar_Disorder_DE$meta.analysis$FDR)
tail(Bipolar_Disorder_DE$meta.analysis$AW.weight)
tail(Bipolar_Disorder_DE$meta.analysis$pval)

#Huntingdons_Disease
Huntingdons_Disease_DE<-MetaDE.rawdata(Huntingdons_Disease_QC,
                                       ind.method = rep("modt", length(names(Huntingdons_Disease_QC))),# moderate t-test - same as limma
                                       meta.method="AW.OC",
                                       nperm=100,
                                       miss.tol=0,
                                       ind.tail="high")

tail(Huntingdons_Disease_DE$meta.analysis$stat)
tail(Huntingdons_Disease_DE$meta.analysis$FDR)
tail(Huntingdons_Disease_DE$meta.analysis$AW.weight)
tail(Huntingdons_Disease_DE$meta.analysis$pval)

#Major_Depressive_Disorder
Major_Depressive_Disorder_DE<-MetaDE.rawdata(Major_Depressive_Disorder_QC,
                                             ind.method = rep("modt", length(names(Major_Depressive_Disorder_QC))),# moderate t-test - same as limma
                                             meta.method="AW.OC",
                                             nperm=100,
                                             miss.tol=0,
                                             ind.tail="high")

tail(Major_Depressive_Disorder_DE$meta.analysis$stat)
tail(Major_Depressive_Disorder_DE$meta.analysis$FDR)
tail(Major_Depressive_Disorder_DE$meta.analysis$AW.weight)
tail(Major_Depressive_Disorder_DE$meta.analysis$pval)

#Schizophrenia
Schizophrenia_DE<-MetaDE.rawdata(Schizophrenia_QC,
                                 ind.method = rep("modt", length(names(Schizophrenia_QC))),# moderate t-test - same as limma
                                 meta.method="AW.OC",
                                 nperm=100,
                                 miss.tol=0,
                                 ind.tail="high")


tail(Schizophrenia_DE$meta.analysis$stat)
tail(Schizophrenia_DE$meta.analysis$FDR)
tail(Schizophrenia_DE$meta.analysis$AW.weight)
tail(Schizophrenia_DE$meta.analysis$pval)

##### HEATMAP "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC" #####

setwd(DE_Heamaps_dir)

pdf("Bipolar_Disorder_DE_heatmaps.pdf")
heatmap.sig.genes(Bipolar_Disorder_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe AW.OC DEG ", line = 1.5)
dev.off()

#parietal lobe
pdf("Huntingdons_Disease_DE_heatmaps.pdf")
heatmap.sig.genes(Huntingdons_Disease_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe AW.OC DEG ", line = 1.5)
dev.off()

# cerebellum
pdf("Major_Depressive_Disorder_DE_heatmaps.pdf")
heatmap.sig.genes(Major_Depressive_Disorder_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum AW.OC DEG ", line = 1.5)
dev.off()

#Temporal Lobe
pdf("Schizophrenia_DE_heatmaps.pdf")
heatmap.sig.genes(Schizophrenia_DE, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe AW.OC DEG ", line = 1.5)
dev.off()

##### COUNT DE NUMBERS #####

count.DEnumber(Bipolar_Disorder_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Huntingdons_Disease_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Schizophrenia_DE,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Major_Depressive_Disorder_DE,
               p.cut=0.05,
               q.cut=0.05)

##### EXTRACT DEG PER BRAIN REGION ####

# using  ind.tail="high" and ind.tail="low" gives same number of DEG using all meta-analysis methods - using only up-genes object

Bipolar_Disorder_DE_AW.OC<-(as.data.frame(subset(Bipolar_Disorder_DE$meta.analysis$FDR, Bipolar_Disorder_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Huntingdons_Disease_DE_AW.OC<-(as.data.frame(subset(Huntingdons_Disease_DE$meta.analysis$FDR, Huntingdons_Disease_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Major_Depressive_Disorder_DE_AW.OC<-(as.data.frame(subset(Major_Depressive_Disorder_DE$meta.analysis$FDR, Major_Depressive_Disorder_DE$meta.analysis$FDR[,1]<=0.05)))[1]

Schizophrenia_DE_AW.OC<-(as.data.frame(subset(Schizophrenia_DE$meta.analysis$FDR, Schizophrenia_DE$meta.analysis$FDR[,1]<=0.05)))[1]

dim(Schizophrenia_DE_AW.OC)
head(Schizophrenia_DE_AW.OC)

dim(Bipolar_Disorder_DE_AW.OC)
head(Bipolar_Disorder_DE_AW.OC)

dim(Huntingdons_Disease_DE_AW.OC)
head(Huntingdons_Disease_DE_AW.OC)

dim(Major_Depressive_Disorder_DE_AW.OC)
head(Major_Depressive_Disorder_DE_AW.OC)

#apply function to:

# 1.extract weights
# 2.subset limma results to match studies in meta analysis
# 3.according to AW weights - make any study where weight=0 to 0 in individual study stats
# 4.keep only values that have contributed to AW.OC meta-analysis.
# 5.make new column to count number of negative logFC values across studies per entrez gene id
# 6.make new column to count number of positive logFC values across studies per entrez gene id 
# 7.decide if genes are up/down regulated

# apply function
Bipolar_Disorder_DEG_regulation<-calculate_AW_gene_regulation(Bipolar_Disorder_DE,log_FC_results)
Bipolar_Disorder_DEG_regulation[rownames(subset(Bipolar_Disorder_DEG_regulation, outcome=="problem")),]

Huntingdons_Disease_DEG_regulation<-calculate_AW_gene_regulation(Huntingdons_Disease_DE,log_FC_results)
Huntingdons_Disease_DEG_regulation[rownames(subset(Huntingdons_Disease_DEG_regulation, outcome=="problem")),]

Schizophrenia_DEG_regulation<-calculate_AW_gene_regulation(Schizophrenia_DE,log_FC_results)
Schizophrenia_DEG_regulation[rownames(subset(Schizophrenia_DEG_regulation, outcome=="problem")),]

Major_Depressive_Disorder_DEG_regulation<-calculate_AW_gene_regulation(Major_Depressive_Disorder_DE,log_FC_results)
Major_Depressive_Disorder_DEG_regulation[rownames(subset(Major_Depressive_Disorder_DEG_regulation, outcome=="problem")),]

# keep genes which are not problematic

Bipolar_Disorder_DEG_regulation_good_genes<-subset(Bipolar_Disorder_DEG_regulation, outcome!="problem")
Huntingdons_Disease_DEG_regulation_good_genes<-subset(Huntingdons_Disease_DEG_regulation, outcome!="problem")
Schizophrenia_DEG_regulation_good_genes<-subset(Schizophrenia_DEG_regulation, outcome!="problem")
Major_Depressive_Disorder_DEG_regulation_good_genes<-subset(Major_Depressive_Disorder_DEG_regulation, outcome!="problem")

dim(Bipolar_Disorder_DEG_regulation)
dim(Bipolar_Disorder_DEG_regulation_good_genes)

dim(Huntingdons_Disease_DEG_regulation)
dim(Huntingdons_Disease_DEG_regulation_good_genes)

dim(Schizophrenia_DEG_regulation)
dim(Schizophrenia_DEG_regulation_good_genes)

dim(Major_Depressive_Disorder_DEG_regulation)
dim(Major_Depressive_Disorder_DEG_regulation_good_genes)

# merge meta AW p value - FDR adjusted

Bipolar_Disorder_DEG_regulation_good_genes<-merge(Bipolar_Disorder_DEG_regulation_good_genes, Bipolar_Disorder_DE_AW.OC, by="row.names")
rownames(Bipolar_Disorder_DEG_regulation_good_genes)<-Bipolar_Disorder_DEG_regulation_good_genes$Row.names
Bipolar_Disorder_DEG_regulation_good_genes$Row.names<-NULL
dim(Bipolar_Disorder_DEG_regulation_good_genes)

Huntingdons_Disease_DEG_regulation_good_genes<-merge(Huntingdons_Disease_DEG_regulation_good_genes, Huntingdons_Disease_DE_AW.OC, by="row.names")
rownames(Huntingdons_Disease_DEG_regulation_good_genes)<-Huntingdons_Disease_DEG_regulation_good_genes$Row.names
Huntingdons_Disease_DEG_regulation_good_genes$Row.names<-NULL
dim(Huntingdons_Disease_DEG_regulation_good_genes)

Schizophrenia_DEG_regulation_good_genes<-merge(Schizophrenia_DEG_regulation_good_genes, Schizophrenia_DE_AW.OC, by="row.names")
rownames(Schizophrenia_DEG_regulation_good_genes)<-Schizophrenia_DEG_regulation_good_genes$Row.names
Schizophrenia_DEG_regulation_good_genes$Row.names<-NULL
dim(Schizophrenia_DEG_regulation_good_genes)

Major_Depressive_Disorder_DEG_regulation_good_genes<-merge(Major_Depressive_Disorder_DEG_regulation_good_genes, Major_Depressive_Disorder_DE_AW.OC, by="row.names")
rownames(Major_Depressive_Disorder_DEG_regulation_good_genes)<-Major_Depressive_Disorder_DEG_regulation_good_genes$Row.names
Major_Depressive_Disorder_DEG_regulation_good_genes$Row.names<-NULL
dim(Major_Depressive_Disorder_DEG_regulation_good_genes)


setwd(AW_results_dir)

# BIPOLAR DISORDER 

#extract AW results and outcome
Bipolar_Disorder_AW <- Bipolar_Disorder_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Bipolar_Disorder_DEG_regulation_good_genes))]
#move rownames to entrez id column
Bipolar_Disorder_AW$Entrez_Gene_ID<-rownames(Bipolar_Disorder_AW)
#rearrange columns
Bipolar_Disorder_AW<-Bipolar_Disorder_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Bipolar_Disorder_AW)[3]<-"Gene_Regulation"
#order by p val
Bipolar_Disorder_AW<-Bipolar_Disorder_AW[order(Bipolar_Disorder_AW$AW.OC),]
#check
head(Bipolar_Disorder_AW)
nrow(Bipolar_Disorder_AW)

#write
write.table(Bipolar_Disorder_AW, file="Bipolar_Disorder_AW.OC", row.names=F, sep="\t", quote=F)

# HUNTINGDON'S DISEASE

#extract AW results and outcome
Huntingdons_Disease_AW <- Huntingdons_Disease_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Huntingdons_Disease_DEG_regulation_good_genes))]
#move rownames to entrez id column
Huntingdons_Disease_AW$Entrez_Gene_ID<-rownames(Huntingdons_Disease_AW)
#rearrange columns
Huntingdons_Disease_AW<-Huntingdons_Disease_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Huntingdons_Disease_AW)[3]<-"Gene_Regulation"
#order by p val
Huntingdons_Disease_AW<-Huntingdons_Disease_AW[order(Huntingdons_Disease_AW$AW.OC),]
#check
head(Huntingdons_Disease_AW)
nrow(Huntingdons_Disease_AW)

#write
write.table(Huntingdons_Disease_AW, file="Huntingdons_Disease_AW.OC", row.names=F, sep="\t", quote=F)

# SCHIZOPHRENIA

#extract AW results and outcome
Schizophrenia_AW <- Schizophrenia_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Schizophrenia_DEG_regulation_good_genes))]
#move rownames to entrez id column
Schizophrenia_AW$Entrez_Gene_ID<-rownames(Schizophrenia_AW)
#rearrange columns
Schizophrenia_AW<-Schizophrenia_AW[,c(3,2,1)]
#change colname to "Gene_regulation"
colnames(Schizophrenia_AW)[3]<-"Gene_Regulation"
#order by p val
Schizophrenia_AW<-Schizophrenia_AW[order(Schizophrenia_AW$AW.OC),]
#check
head(Schizophrenia_AW)
nrow(Schizophrenia_AW)

#write
write.table(Schizophrenia_AW, file="Schizophrenia_AW.OC", row.names=F, sep="\t", quote=F)

# MAJOR DEPRESSIVE DISORDER

# #extract AW results and outcome
# Major_Depressive_Disorder_AW <- Major_Depressive_Disorder_DEG_regulation_good_genes[,grep("AW.OC|outcome", colnames(Major_Depressive_Disorder_DEG_regulation_good_genes))]
# #move rownames to entrez id column
# Major_Depressive_Disorder_AW$Entrez_Gene_ID<-rownames(Major_Depressive_Disorder_AW)
# #rearrange columns
# Major_Depressive_Disorder_AW<-Major_Depressive_Disorder_AW[,c(3,2,1)]
# #change colname to "Gene_regulation"
# colnames(Major_Depressive_Disorder_AW)[3]<-"Gene_Regulation"
# #order by p val
# Major_Depressive_Disorder_AW<-Major_Depressive_Disorder_AW[order(Major_Depressive_Disorder_AW$AW.OC),]
# #check
# head(Major_Depressive_Disorder_AW)
# nrow(Major_Depressive_Disorder_AW)
# 
# #write
# write.table(Major_Depressive_Disorder_AW, file="Major_Depressive_Disorder_AW.OC", row.names=F, sep="\t", quote=F)
# 

##### SAVE IMAGE #####

setwd(work_dir)

save.image("NON-AD_Meta_Analysis.Rdata")
