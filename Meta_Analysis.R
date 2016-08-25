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
##
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/Clean_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/6.Meta_Analysis/"

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

# create entrez gene id to gene symbol lookup

# convert entrez gene id to gene symbol
entrez_gene_symbol <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(entrez_gene_symbol)
# Convert to a list
gene_symbol_lookup_table <- as.data.frame(entrez_gene_symbol[mapped_genes])

#check
head(gene_symbol_lookup_table)

#check for unique entrez gene_id
anyDuplicated(gene_symbol_lookup_table$gene_id)
anyDuplicated(gene_symbol_lookup_table$symbol)

# convert data to "MetaOmics" format - keep only diagnosis- create function for conversion of data

convert_data_to_MetaOmics_format<-function(dataset){
  # remove unwanted pheno - keep only diagnosis + expression (merge by rownames)
  dataset<-merge(dataset[5], dataset[10:dim(dataset)[2]], by="row.names")
  rownames(dataset)<-dataset$Row.names
  dataset$Row.names<-NULL
  # recode diagnosis column - case=1, control=0
  dataset[dataset$Diagnosis=="Control",1]<-"0"
  dataset[dataset$Diagnosis=="AD",1]<-"1"
  # relabel Diagnosis to Label
  colnames(dataset)[1]<-"label"
  # chane dataframe to numeric and transpose
  dataset<-as.data.frame(t(data.matrix(dataset)))
  # add gene symbol
  dataset$SYMBOL<-NA
  # loop through each entrez gene id, and extract gene symbol - NA if no gene symbol
  for (x in 1:dim(dataset)[1]){
    #if gene entrez id in gene_symbol_lookup_table extract corresponding gene symbol
    if(rownames(dataset)[x] %in% gene_symbol_lookup_table$gene_id == T){
      dataset$SYMBOL[x]<-gene_symbol_lookup_table[gene_symbol_lookup_table$gene_id==rownames(dataset)[x],2]
    }
    else{
      dataset$SYMBOL[x]<-NA
    }
  }
  # move SYMBOL to 1st column
  dataset<-dataset %>% select(SYMBOL, everything())  
  print(head(dataset)[1:5])
  return(dataset)
}

# apply function to all datasets - leave out occipital lobe

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
                       matched=F) # matched format false - i.e need to match entrez id to gene symbol

Parietal_Lobe<-MetaDE.Read(c("AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format", 
                            "AMP_MSBB_U133B_Superior_Parietal_Lobule_MetaOmics_format", 
                            "E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
                            "E_GEOD_5281_Posterior_Singulate_MetaOmics_format"), # 4 datasets
                          via="txt", # txt format
                          skip=rep(1,4), # skip 1 row - 4 datasets
                          log=T, # data is log transformed
                          matched=F) # matched format false - i.e need to match entrez id to gene symbol

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
                           matched=F) # matched format false - i.e need to match entrez id to gene symbol

Cerebellum<-MetaDE.Read(c("inhouse_data_Cerebellum_MetaOmics_format", 
                          "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"), # 2 datasets
                           via="txt", # txt format
                           skip=rep(1,2), # skip 1 row - 2 datasets
                           log=T, # data is log transformed
                           matched=F) # matched format false - i.e need to match entrez id to gene symbol

setwd(work_dir)

##### MAP ENTREZ TO GENE SYMBOL #####

setwd(work_dir)

# map gene symbol to entrez id - frontal lobe

Frontal_Lobe_GS<-MetaDE.match(Frontal_Lobe, pool.replicate="IQR")

Parietal_Lobe_GS<-MetaDE.match(Parietal_Lobe, pool.replicate="IQR")

Temporal_Lobe_GS<-MetaDE.match(Temporal_Lobe, pool.replicate="IQR")

Cerebellum_GS<-MetaDE.match(Cerebellum, pool.replicate="IQR")

##### CREATE FORK - EXTRACT COMMON PROBES #####

# common probes

Frontal_Lobe_GS_merged<-MetaDE.merge(Frontal_Lobe_GS, MVperc=0)
Parietal_Lobe_GS_merged<-MetaDE.merge(Parietal_Lobe_GS, MVperc=0)
Temporal_Lobe_GS_merged<-MetaDE.merge(Temporal_Lobe_GS, MVperc=0)
Cerebellum_GS_merged<-MetaDE.merge(Cerebellum_GS, MVperc=0)

dim(Frontal_Lobe_GS_merged[[1]][[1]])
dim(Parietal_Lobe_GS_merged[[1]][[1]])
dim(Temporal_Lobe_GS_merged[[1]][[1]])
dim(Cerebellum_GS_merged[[1]][[1]])

# no common probes - full data merge - minimum 20% common required due to knn.impute function

Frontal_Lobe_GS_full_merged<-MetaDE.merge(Frontal_Lobe_GS, MVperc=0.8)
Parietal_Lobe_GS_full_merged<-MetaDE.merge(Parietal_Lobe_GS, MVperc=0.8)
Temporal_Lobe_GS_full_merged<-MetaDE.merge(Temporal_Lobe_GS, MVperc=0.8)
Cerebellum_GS_full_merged<-MetaDE.merge(Cerebellum_GS, MVperc=0.8)

dim(Frontal_Lobe_GS_full_merged[[1]][[1]])
dim(Parietal_Lobe_GS_full_merged[[1]][[1]])
dim(Temporal_Lobe_GS_full_merged[[1]][[1]])
dim(Cerebellum_GS_full_merged[[1]][[1]])

##### FILTER GENES #####

#filter out 30% un-expressed genes and then 30% non-informative genes across datasets

# filter genes - on common probes

Frontal_Lobe_GS_merged_Filtered<-MetaDE.filter(Frontal_Lobe_GS_merged,c(0.3,0.3)) 
Parietal_Lobe_GS_merged_Filtered<-MetaDE.filter(Parietal_Lobe_GS_merged,c(0.3,0.3)) 
Temporal_Lobe_GS_merged_Filtered<-MetaDE.filter(Temporal_Lobe_GS_merged,c(0.3,0.3)) 
Cerebellum_GS_merged_Filtered<-MetaDE.filter(Cerebellum_GS_merged,c(0.3,0.3)) 

dim(Frontal_Lobe_GS_merged_Filtered[[1]][[1]])
dim(Parietal_Lobe_GS_merged_Filtered[[1]][[1]])
dim(Temporal_Lobe_GS_merged_Filtered[[1]][[1]])
dim(Cerebellum_GS_merged_Filtered[[1]][[1]])

# fiter genes on full data

Frontal_Lobe_GS_full_merged_Filtered<-MetaDE.filter(Frontal_Lobe_GS_full_merged,c(0.3,0.3)) 
Parietal_Lobe_GS_full_merged_Filtered<-MetaDE.filter(Parietal_Lobe_GS_full_merged,c(0.3,0.3)) 
Temporal_Lobe_GS_full_merged_Filtered<-MetaDE.filter(Temporal_Lobe_GS_full_merged,c(0.3,0.3)) 
Cerebellum_GS_full_merged_Filtered<-MetaDE.filter(Cerebellum_GS_full_merged,c(0.3,0.3)) 

dim(Frontal_Lobe_GS_full_merged_Filtered[[1]][[1]])
dim(Parietal_Lobe_GS_full_merged_Filtered[[1]][[1]])
dim(Temporal_Lobe_GS_full_merged_Filtered[[1]][[1]])
dim(Cerebellum_GS_full_merged_Filtered[[1]][[1]])

##### META QC - CREATE QC OBJECT #####

# meta QC

# create QC object with the following code:
 
Frontal_Lobe_common_QC<-list()
Parietal_Lobe_common_QC<-list()
Temporal_Lobe_common_QC<-list()
Cerebellum_common_QC<-list()

# moves diagnosis label to column name - Frontal_Lobe_GS_merged_Filtered 
for(i in 1:14){ 
  colnames(Frontal_Lobe_GS_merged_Filtered[[i]][[1]])<-Frontal_Lobe_GS_merged_Filtered[[i]][[2]]
  Frontal_Lobe_common_QC[[i]]<-impute.knn(Frontal_Lobe_GS_merged_Filtered[[i]][[1]])$data
}

names(Frontal_Lobe_common_QC)<-names(Frontal_Lobe_GS_merged_Filtered)

head(Frontal_Lobe_common_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)[,1:5]
dim(Frontal_Lobe_common_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)

# moves diagnosis label to column name - Parietal_Lobe_GS_merged_Filtered
for(i in 1:4){ 
  colnames(Parietal_Lobe_GS_merged_Filtered[[i]][[1]])<-Parietal_Lobe_GS_merged_Filtered[[i]][[2]]
  Parietal_Lobe_common_QC[[i]]<-impute.knn(Parietal_Lobe_GS_merged_Filtered[[i]][[1]])$data
}

names(Parietal_Lobe_common_QC)<-names(Parietal_Lobe_GS_merged_Filtered)

head(Parietal_Lobe_common_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)[,1:5]
dim(Parietal_Lobe_common_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)

# moves diagnosis label to column name - Temporal_Lobe_GS_merged_Filtered
for(i in 1:26){ 
  colnames(Temporal_Lobe_GS_merged_Filtered[[i]][[1]])<-Temporal_Lobe_GS_merged_Filtered[[i]][[2]]
  Temporal_Lobe_common_QC[[i]]<-impute.knn(Temporal_Lobe_GS_merged_Filtered[[i]][[1]])$data
}

names(Temporal_Lobe_common_QC)<-names(Temporal_Lobe_GS_merged_Filtered)

head(Temporal_Lobe_common_QC$E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format)[,1:5]
dim(Temporal_Lobe_common_QC$E_GEOD_48350_Entorhinal_Cortex_MetaOmics_format)

# moves diagnosis label to column name - Cerebellum_GS_merged_Filtered
for(i in 1:2){ 
  colnames(Cerebellum_GS_merged_Filtered[[i]][[1]])<-Cerebellum_GS_merged_Filtered[[i]][[2]]
  Cerebellum_common_QC[[i]]<-impute.knn(Cerebellum_GS_merged_Filtered[[i]][[1]])$data
}

names(Cerebellum_common_QC)<-names(Cerebellum_GS_merged_Filtered)

head(Cerebellum_common_QC$inhouse_data_Cerebellum_MetaOmics_format)[,1:5]
dim(Cerebellum_common_QC$inhouse_data_Cerebellum_MetaOmics_format)

# same for full data

Frontal_Lobe_full_QC<-list()
Parietal_Lobe_full_QC<-list()
Temporal_Lobe_full_QC<-list()
Cerebellum_full_QC<-list()

# moves diagnosis label to column name - Frontal_Lobe_GS_full_merged_Filtered - impute missing - 
for(i in 1:14){ 
  colnames(Frontal_Lobe_GS_full_merged_Filtered[[i]][[1]])<-Frontal_Lobe_GS_full_merged_Filtered[[i]][[2]]
  Frontal_Lobe_full_QC[[i]]<-impute.knn(Frontal_Lobe_GS_full_merged_Filtered[[i]][[1]])$data
}

names(Frontal_Lobe_full_QC)<-names(Frontal_Lobe_GS_full_merged_Filtered)

head(Frontal_Lobe_full_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)[,1:5]
dim(Frontal_Lobe_full_QC$E_GEOD_36980_Frontal_Cortex_MetaOmics_format)

# moves diagnosis label to column name - Parietal_Lobe_GS_full_merged_Filtered - impute missing - 
for(i in 1:4){ 
  colnames(Parietal_Lobe_GS_full_merged_Filtered[[i]][[1]])<-Parietal_Lobe_GS_full_merged_Filtered[[i]][[2]]
  Parietal_Lobe_full_QC[[i]]<-impute.knn(Parietal_Lobe_GS_full_merged_Filtered[[i]][[1]])$data
}

names(Parietal_Lobe_full_QC)<-names(Parietal_Lobe_GS_full_merged_Filtered)

head(Parietal_Lobe_full_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)[,1:5]
dim(Parietal_Lobe_full_QC$AMP_MSBB_U133A_Superior_Parietal_Lobule_MetaOmics_format)

# moves diagnosis label to column name - Temporal_Lobe_GS_full_merged_Filtered - impute missing - 
for(i in 1:26){ 
  colnames(Temporal_Lobe_GS_full_merged_Filtered[[i]][[1]])<-Temporal_Lobe_GS_full_merged_Filtered[[i]][[2]]
  Temporal_Lobe_full_QC[[i]]<-impute.knn(Temporal_Lobe_GS_full_merged_Filtered[[i]][[1]])$data
}

names(Temporal_Lobe_full_QC)<-names(Temporal_Lobe_GS_full_merged_Filtered)

head(Temporal_Lobe_full_QC$E_GEOD_28146_Hippocampus_CA1_MetaOmics_format)[,1:5]
dim(Temporal_Lobe_full_QC$E_GEOD_28146_Hippocampus_CA1_MetaOmics_format)

# moves diagnosis label to column name - Cerebellum_GS_full_merged_Filtered - impute missing - 
for(i in 1:2){ 
  colnames(Cerebellum_GS_full_merged_Filtered[[i]][[1]])<-Cerebellum_GS_full_merged_Filtered[[i]][[2]]
  Cerebellum_full_QC[[i]]<-impute.knn(Cerebellum_GS_full_merged_Filtered[[i]][[1]])$data
}

names(Cerebellum_full_QC)<-names(Cerebellum_GS_full_merged_Filtered)

head(Cerebellum_full_QC$inhouse_data_Cerebellum_MetaOmics_format)[,1:5]
dim(Cerebellum_full_QC$inhouse_data_Cerebellum_MetaOmics_format)

##### RUN QC CHECK #####

# create function to run quantitative quality control measures

run_meta_qc<-function(dataset) {
  MetaQC(dataset, # data
  #"c2.cp.biocarta.v3.0.symbols.gmt", #pathway file from biocarta
  "c2.cp.biocarta.v5.1.symbols.gmt", # latest version downloaded 15/08/2016 
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

Frontal_Lobe_common_QC_check<-run_meta_qc(Frontal_Lobe_common_QC)
Parietal_Lobe_common_QC_check<-run_meta_qc(Parietal_Lobe_common_QC)
Temporal_Lobe_common_QC_check<-run_meta_qc(Temporal_Lobe_common_QC)
Cerebellum_common_QC_check<-run_meta_qc(Cerebellum_common_QC)
Frontal_Lobe_full_QC_check<-run_meta_qc(Frontal_Lobe_full_QC)
Parietal_Lobe_full_QC_check<-run_meta_qc(Parietal_Lobe_full_QC)
Temporal_Lobe_full_QC_check<-run_meta_qc(Temporal_Lobe_full_QC)
Cerebellum_full_QC_check<-run_meta_qc(Cerebellum_full_QC)

Frontal_Lobe_common_QC_check
Parietal_Lobe_common_QC_check
Temporal_Lobe_common_QC_check
Cerebellum_common_QC_check
Frontal_Lobe_full_QC_check
Parietal_Lobe_full_QC_check
Temporal_Lobe_full_QC_check
Cerebellum_full_QC_check

# create function to run 2nd part of QC

# runQC- **script freezes when 1st downloading fileForCQCp. downloaded manually and inserted into working folder http://software.broadinstitute.org/gsea/downloads.jsp

run_meta_qc2<-function(dataset){
  runQC(dataset, 
        nPath=NULL, # 
        B=1e4, # The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
        pvalCut=.05, #P-value threshold used for AQC calculation.
        #pvalAdjust=T, #whether to apply p-value adjustment due to multiple testing (B-H procedure is used).
        #fileForCQCp="c2.all.v3.0.symbols.gmt")
        fileForCQCp="c2.all.v5.1.symbols.gmt") # latest version of c2.all.v5.1.symbols.gmt downloaded 15/08/2016
}
 
# ** ERROR IN QC - ? LOW NUMBER OF COMMON GENES

#run_meta_qc2(Frontal_Lobe_common_QC_check)
#run_meta_qc2(Parietal_Lobe_common_QC_check)
#run_meta_qc2(Temporal_Lobe_common_QC_check)
#run_meta_qc2(Cerebellum_common_QC_check)

run_meta_qc2(Frontal_Lobe_full_QC_check)
run_meta_qc2(Parietal_Lobe_full_QC_check)
run_meta_qc2(Temporal_Lobe_full_QC_check)
#run_meta_qc2(Cerebellum_full_QC_check) - "Error in { : task 1 failed - "not enough 'y' observations"

#Frontal_Lobe_common_QC_check
#Parietal_Lobe_common_QC_check
#Temporal_Lobe_common_QC_check
#Cerebellum_common_QC_check

Frontal_Lobe_full_QC_check
Parietal_Lobe_full_QC_check
Temporal_Lobe_full_QC_check
#Cerebellum_full_QC_check

#plot PCA

# plot(Frontal_Lobe_common_QC_check)
# plot(Parietal_Lobe_common_QC_check)
# plot(Temporal_Lobe_common_QC_check)
# plot(Cerebellum_common_QC_check)
plot(Frontal_Lobe_full_QC_check)
plot(Parietal_Lobe_full_QC_check)
plot(Temporal_Lobe_full_QC_check)
# plot(Cerebellum_full_QC_check)

##### PLOT QC PLOTS TO PDF #####

setwd(MetaQC_plots_dir)

pdf("Frontal_lobe_MetaQC.pdf")
plot(Frontal_Lobe_full_QC_check)
mtext("Frontal Lobe Datasets", cex=2, line = 1.5)
dev.off()

pdf("Parietal_lobe_MetaQC.pdf")
plot(Parietal_Lobe_full_QC_check)
mtext("Parietal Lobe Datasets", cex=2, line = 1.5)
dev.off()

pdf("Temporal_lobe_MetaQC.pdf")
plot(Temporal_Lobe_full_QC_check)
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
                        matched=F) # matched format false - i.e need to match entrez id to gene symbol
  
  #setwd to work dir
  setwd(work_dir)
  # map gene symbol to entrez id - frontal lobe
  raw_data_GS<-MetaDE.match(raw_data, pool.replicate="IQR")
  
  # merge by common probes, 0=100% common, 0.8=20% common
  raw_data_GS_merged<-MetaDE.merge(raw_data_GS, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_GS_merged[[1]][[1]])[1]), sep=" ")
  
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_GS_merged_Filtered<-MetaDE.filter(raw_data_GS_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_GS_merged_Filtered[[1]][[1]])[1]), sep=" ")
  
  # create QC object 
  Data_QC<-list()
  
  # moves diagnosis label to column name - Frontal_Lobe_GS_merged_Filtered 
  for(i in 1:number_of_datasets){ 
    colnames(raw_data_GS_merged_Filtered[[i]][[1]])<-raw_data_GS_merged_Filtered[[i]][[2]]
    Data_QC[[i]]<-impute.knn(raw_data_GS_merged_Filtered[[i]][[1]])$data
  }
  names(Data_QC)<-names(raw_data_GS_merged_Filtered)
  
  # run quantitative quality control measures
  Data_QC_measures<-MetaQC(Data_QC, # data
                           #"c2.cp.biocarta.v3.0.symbols.gmt", #pathway file from biocarta
                           "c2.cp.biocarta.v5.1.symbols.gmt", # latest version downloaded 15/08/2016 
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
  #run 2nd part of QC
  runQC(Data_QC_measures, 
        nPath=NULL, # 
        B=1e4, # The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
        pvalCut=.05, #P-value threshold used for AQC calculation.
        #pvalAdjust=T, #whether to apply p-value adjustment due to multiple testing (B-H procedure is used).
        #fileForCQCp="c2.all.v3.0.symbols.gmt")
        fileForCQCp="c2.all.v5.1.symbols.gmt") # latest version of c2.all.v5.1.symbols.gmt downloaded 15/08/2016
  plot(Data_QC_measures)
  return(Data_QC_measures)
}

# apply function to lobe regions without AMP MSBB datasets

Frontal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_36980_Frontal_Cortex_MetaOmics_format", 
                                              "E_GEOD_48350_Superior_Frontal_Gyrus_MetaOmics_format", 
                                              "E_GEOD_5281_Superior_Frontal_Gyrus_MetaOmics_format", 
                                              "inhouse_data_Frontal_Cortex_MetaOmics_format"),
                                   MVperc=0.2,
                                   gene_filter=c(0.3,0.3))

# Parietal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_48350_Postcentral_Gyrus_MetaOmics_format", 
#                                                "E_GEOD_5281_Posterior_Singulate_MetaOmics_format"),
#                                     MVperc=0.2,
#                                     gene_filter=c(0.3,0.3))

                                    
Temporal_lobe_no_AMP_QC<-run_MetaQC(datasets=c("E_GEOD_28146_Hippocampus_CA1_MetaOmics_format", 
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

# Cerebellum_no_AMP_QC<-run_MetaQC(datasets=c("inhouse_data_Cerebellum_MetaOmics_format", 
#                                   "AMP_MAYOeGWAS_Cerebellum_MetaOmics_format"),
#                        MVperc=0.2,
#                        gene_filter=c(0.3,0.3))

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
  # map gene symbol to entrez id - frontal lobe
  raw_data_GS<-MetaDE.match(raw_data, pool.replicate="IQR")
  
  # merge by common probes, 0=100% common, 0.8=20% common
  raw_data_GS_merged<-MetaDE.merge(raw_data_GS, MVperc=MVperc)
  print(paste("number of genes after merging: ", dim(raw_data_GS_merged[[1]][[1]])[1]), sep=" ")
  
  # gene_filter = i.e c(0.3, 0.3) filter out 30% un-expressed genes and then 30% non-informative genes across datasets, 
  raw_data_GS_merged_Filtered<-MetaDE.filter(raw_data_GS_merged, gene_filter) 
  print(paste("number of genes after filtering un-expressed and non-informative genes: ",dim(raw_data_GS_merged_Filtered[[1]][[1]])[1]), sep=" ")
  return(raw_data_GS_merged_Filtered)
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

##### META - DE ANALYSIS UPREGULATED #####

# Frontal Lobe
Frontal_lobe_DE_up<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
                                   ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                   meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                   nperm=300,
                                   miss.tol=0,
                                   ind.tail="high")

tail(Frontal_lobe_DE_up$meta.analysis$stat)
tail(Frontal_lobe_DE_up$meta.analysis$FDR)
tail(Frontal_lobe_DE_up$meta.analysis$AW.weight)
tail(Frontal_lobe_DE_up$meta.analysis$pval)

#Parietal Lobe
Parietal_lobe_DE_up<-MetaDE.rawdata(Parietal_lobe_no_AMP_QC_filtered,
                                    ind.method = rep("modt", length(names(Parietal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                    meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                    nperm=300,
                                    miss.tol=0,
                                    ind.tail="high")

tail(Parietal_lobe_DE_up$meta.analysis$stat)
tail(Parietal_lobe_DE_up$meta.analysis$FDR)
tail(Parietal_lobe_DE_up$meta.analysis$AW.weight)
tail(Parietal_lobe_DE_up$meta.analysis$pval)

#Cerebellum
Cerebellum_DE_up<-MetaDE.rawdata(Cerebellum_no_AMP_QC_filtered,
                                 ind.method = rep("modt", length(names(Cerebellum_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                 meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                 nperm=300,
                                 miss.tol=0,
                                 ind.tail="high")

tail(Cerebellum_DE_up$meta.analysis$stat)
tail(Cerebellum_DE_up$meta.analysis$FDR)
tail(Cerebellum_DE_up$meta.analysis$AW.weight)
tail(Cerebellum_DE_up$meta.analysis$pval)

# #Temporal - requires 251.1G
# Temporal_lobe_DE_up<-MetaDE.rawdata(Temporal_lobe_no_AMP_QC_filtered,
#                                     ind.method = rep("modt", length(names(Temporal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                     meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
#                                     nperm=20,
#                                     miss.tol=0,
#                                     ind.tail="high")

# # temporal lobe fails to aquire 251 ram reqiured - save and run on cluster
# # setwd(work_dir)
# # save(Temporal_lobe_no_AMP_QC_filtered, file="Temporal_lobe_no_AMP_QC_filtered.Rdata")

#read in  temporal DE run on cluster

load("Meta_DE/Temporal_lobe_run_on_cluster/Temporal_lobe_DE_up.Rdata")

tail(Temporal_lobe_DE_up$meta.analysis$stat)
tail(Temporal_lobe_DE_up$meta.analysis$FDR)
tail(Temporal_lobe_DE_up$meta.analysis$AW.weight)
tail(Temporal_lobe_DE_up$meta.analysis$pval)

##### META - DE ANALYSIS DOWN-REGULATED #####

# Frontal Lobe
Frontal_lobe_DE_down<-MetaDE.rawdata(Frontal_lobe_no_AMP_QC_filtered,
                                     ind.method = rep("modt", length(names(Frontal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                     meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                     nperm=300,
                                     miss.tol=0,
                                     ind.tail="low")

tail(Frontal_lobe_DE_down$meta.analysis$stat)
tail(Frontal_lobe_DE_down$meta.analysis$FDR)
tail(Frontal_lobe_DE_down$meta.analysis$AW.weight)
tail(Frontal_lobe_DE_down$meta.analysis$pval)

#Parietal Lobe
Parietal_lobe_DE_down<-MetaDE.rawdata(Parietal_lobe_no_AMP_QC_filtered,
                                      ind.method = rep("modt", length(names(Parietal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                      meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                      nperm=300,
                                      miss.tol=0,
                                      ind.tail="low")

tail(Parietal_lobe_DE_down$meta.analysis$stat)
tail(Parietal_lobe_DE_down$meta.analysis$FDR)
tail(Parietal_lobe_DE_down$meta.analysis$AW.weight)
tail(Parietal_lobe_DE_down$meta.analysis$pval)

#Temporal - requires 251.1G
# Temporal_lobe_DE_down<-MetaDE.rawdata(Temporal_lobe_no_AMP_QC_filtered,
#                                       ind.method = rep("modt", length(names(Temporal_lobe_no_AMP_QC_filtered))),# moderate t-test - same as limma
#                                       meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
#                                       nperm=300,
#                                       miss.tol=0,
#                                       ind.tail="low")
# 
# head(Temporal_lobe_DE_down$meta.analysis$stat)
# head(Temporal_lobe_DE_down$meta.analysis$FDR)
# head(Temporal_lobe_DE_down$meta.analysis$AW.weight)
# head(Temporal_lobe_DE_down$meta.analysis$pval)

#Cerebellum
Cerebellum_DE_down<-MetaDE.rawdata(Cerebellum_no_AMP_QC_filtered,
                                   ind.method = rep("modt", length(names(Cerebellum_no_AMP_QC_filtered))),# moderate t-test - same as limma
                                   meta.method=c("AW.OC", "maxP.OC", "minP.OC", "Fisher.OC"),
                                   nperm=300,
                                   miss.tol=0,
                                   ind.tail="low")

tail(Cerebellum_DE_down$meta.analysis$stat)
tail(Cerebellum_DE_down$meta.analysis$FDR)
tail(Cerebellum_DE_down$meta.analysis$AW.weight)
tail(Cerebellum_DE_down$meta.analysis$pval)

##### CALCULATE EFFECT SIZE - META #####

#This functions is used to calculate the effect size, standardized mean difference
#the p-values are calculated using parametric method by assupming the z-scores following a standard normal distribution

#Frontal Lobe - calculate for individual an then meta - nested
Frontal_lobe_DE_ES<-MetaDE.ES(ind.cal.ES(Frontal_lobe_no_AMP_QC_filtered,
                                         paired=rep(F, length(names(Frontal_lobe_no_AMP_QC_filtered))),
                                         nperm=300,
                                         miss.tol=0),
                              meta.method="REM")

head(Frontal_lobe_DE_ES$pval)
tail(Frontal_lobe_DE_ES$pval)


#Parietal Lobe
Parietal_lobe_DE_ES<-MetaDE.ES(ind.cal.ES(Parietal_lobe_no_AMP_QC_filtered,
                                          paired=rep(F, length(names(Parietal_lobe_no_AMP_QC_filtered))),
                                          nperm=300,
                                          miss.tol=0),
                               meta.method="REM")

head(Parietal_lobe_DE_ES$pval)
tail(Parietal_lobe_DE_ES$pval)

#Temporal Lobe
Temporal_lobe_DE_ES<-MetaDE.ES(ind.cal.ES(Temporal_lobe_no_AMP_QC_filtered,
                                          paired=rep(F, length(names(Temporal_lobe_no_AMP_QC_filtered))),
                                          nperm=300,
                                          miss.tol=0),
                               meta.method="REM")

head(Temporal_lobe_DE_ES$pval)
tail(Temporal_lobe_DE_ES$pval)

#Cerebellum
Cerebellum_DE_ES<-MetaDE.ES(ind.cal.ES(Cerebellum_no_AMP_QC_filtered,
                                       paired=rep(F, length(names(Cerebellum_no_AMP_QC_filtered))),
                                       nperm=300,
                                       miss.tol=0),
                            meta.method="REM")

head(Cerebellum_DE_ES$pval)
tail(Cerebellum_DE_ES$pval)

##### HEATMAP #####

# meta DE analysis done using "AW.OC", "maxP.OC", "minP.OC", "Fisher.OC". plot heat map for each method
setwd(DE_Heamaps_dir)

pdf("Frontal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Frontal_lobe_DE_up, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe AW.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_up, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe maxP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_up, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe minP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_up, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe Fisher.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_down, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe AW.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_down, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe maxP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_down, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe minP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Frontal_lobe_DE_down, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Frontal Lobe Fisher.OC DEG DOWN", line = 1.5)
dev.off()

#parietal lobe
pdf("Parietal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Parietal_lobe_DE_up, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe AW.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_up, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe maxP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_up, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe minP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_up, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe Fisher.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_down, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe AW.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_down, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe maxP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_down, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe minP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Parietal_lobe_DE_down, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Parietal Lobe Fisher.OC DEG DOWN", line = 1.5)
dev.off()

# cerebellum
pdf("Cerebellum_DE_heatmaps.pdf")
heatmap.sig.genes(Cerebellum_DE_up, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum AW.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_up, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum maxP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_up, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum minP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_up, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum Fisher.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_down, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum AW.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_down, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum maxP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_down, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum minP.OC DEG DOWN", line = 1.5)

heatmap.sig.genes(Cerebellum_DE_down, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Cerebellum Fisher.OC DEG DOWN", line = 1.5)
dev.off()

#Temporal Lobe
pdf("Temporal_Lobe_DE_heatmaps.pdf")
heatmap.sig.genes(Temporal_lobe_DE_up, meta.method="AW.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe AW.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE_up, meta.method="maxP.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe maxP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE_up, meta.method="minP.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe minP.OC DEG UP ", line = 1.5)

heatmap.sig.genes(Temporal_lobe_DE_up, meta.method="Fisher.OC", fdr.cut=0.05, color="GR")
mtext("Temporal Lobe Fisher.OC DEG UP ", line = 1.5)
dev.off()

##### COUNT DE NUMBERS #####
# up regulated DE results

count.DEnumber(Frontal_lobe_DE_up,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Parietal_lobe_DE_up,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Temporal_lobe_DE_up,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Cerebellum_DE_up,
               p.cut=0.05,
               q.cut=0.05)

# down regulated DE results

count.DEnumber(Frontal_lobe_DE_down,
               p.cut=0.05,
               q.cut=0.05)

count.DEnumber(Parietal_lobe_DE_down,
               p.cut=0.05,
               q.cut=0.05)

# count.DEnumber(Temporal_lobe_DE_down,
#                p.cut=0.05,
#                q.cut=0.05)

count.DEnumber(Cerebellum_DE_down,
               p.cut=0.05,
               q.cut=0.05)

##### DE NUMBERS PER META_METHOD - PLOT #####

setwd(Meta_DE_dir)

pdf("DE_numbers_for_each_method.pdf")

draw.DEnumber(Frontal_lobe_DE_up,
              0.05,
              FDR=T)
mtext("Frontal Lobe UP DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Frontal_lobe_DE_down,
              0.05,
              FDR=T)
mtext("Frontal Lobe DOWN DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Parietal_lobe_DE_up,
              0.05,
              FDR=T)
mtext("Parietal Lobe UP DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Parietal_lobe_DE_down,
              0.05,
              FDR=T)
mtext("Parietal Lobe DOWN DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Cerebellum_DE_up,
              0.05,
              FDR=T)
mtext("Cerebellum UP DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Cerebellum_DE_down,
              0.05,
              FDR=T)
mtext("Cerebellum DOWN DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

draw.DEnumber(Temporal_lobe_DE_up,
              0.05,
              FDR=T)
mtext("Temporal Lobe UP DEG numbers in indivdual and Meta-analysis", line = 2, cex=1.2)

dev.off()

##### EXTRACT DEG PER BRAIN REGION ####

# using  ind.tail="high" and ind.tail="low" gives same number of DEG using all meta-analysis methods - using only up-genes object

Frontal_lobe_DE_AW.OC<-(as.data.frame(subset(Frontal_lobe_DE_up$meta.analysis$FDR, Frontal_lobe_DE_up$meta.analysis$FDR[,1]<=0.05)))[1]

Parietal_lobe_DE_AW.OC<-(as.data.frame(subset(Parietal_lobe_DE_up$meta.analysis$FDR, Parietal_lobe_DE_up$meta.analysis$FDR[,1]<=0.05)))[1]

Cerebellum_DE_AW.OC<-(as.data.frame(subset(Cerebellum_DE_up$meta.analysis$FDR, Cerebellum_DE_up$meta.analysis$FDR[,1]<=0.05)))[1]

Temporal_lobe_DE_AW.OC<-(as.data.frame(subset(Temporal_lobe_DE_up$meta.analysis$FDR, Temporal_lobe_DE_up$meta.analysis$FDR[,1]<=0.05)))[1]

dim(Temporal_lobe_DE_AW.OC)
head(Temporal_lobe_DE_AW.OC)

dim(Frontal_lobe_DE_AW.OC)
head(Frontal_lobe_DE_AW.OC)

dim(Parietal_lobe_DE_AW.OC)
head(Parietal_lobe_DE_AW.OC)

dim(Cerebellum_DE_AW.OC)
head(Cerebellum_DE_AW.OC)

##### META PATHWAY ANALYSIS #####

setwd(Meta_Pathway_dir)

# pathway database - file downloaded manually from website ("http://software.broadinstitute.org/gsea/downloads.jsp")
pathway_database<-getGmt("c5.all.v5.1.symbols.gmt") 

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
subset(Frontal_Lobe_pathway_MAPE$qvalue, MAPE_I<1.1)

#heatmap of pathways
plotMAPE(Frontal_Lobe_pathway_MAPE, 
         cutoff=1.1, 
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

save.image("Meta_analysis.Rdata")


