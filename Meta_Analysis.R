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
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/Clean_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/6.Meta_Analysis/"

setwd(work_dir)

# create dir DE and patway analysis

dir.create(paste(work_dir,"Meta_DE", sep="/"))

Meta_DE_dir=paste(work_dir,"Meta_DE", sep="/")

dir.create(paste(work_dir,"Meta_Pathway", sep="/"))

Meta_Pathway_dir=paste(work_dir,"Meta_Pathway", sep="/")

dir.create(paste(work_dir,"MetaQC_plots", sep="/"))

MetaQC_plots_dir=paste(work_dir,"MetaQC_plots", sep="/")

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

#create directory
dir.create(paste(work_dir, "Data_in_metaomics_format", sep="/"))

Data_in_metaomics_format=paste(work_dir, "Data_in_metaomics_format", sep="/")

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

##### COMPILE REPORT #####

#rmarkdown::render("Meta_Analysis.R") # make sure this is hashed out and file is saved. run separetly on command line

##### SAVE IMAGE #####

setwd(work_dir)

save.image("Meta_analysis.Rdata")


