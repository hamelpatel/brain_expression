
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                               INDIVIDUAL DE ANALYSIS                                                             #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# DESCRIPTION
# Remove AD where BRAAK<=3
# Remove Control where BRAAK >=3
# Differential expression of AD and control 
# heamap
#

rm(list=ls())

##### SET PARAMETERS #####

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

work_dir="/media/hamel/1TB/Projects/Brain_expression/3.DE_analysis/"

setwd(work_dir)

##### LOAD LIBRARIES #####

library(limma)
library(gplots)
library(RColorBrewer)
library(pamr)

##### LOAD DATASETS - 11 tables #####

E_GEOD_28146<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-28146/Pre-Processing/clean_data/E-GEOD-28146_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_29378<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-29378/Pre-Processing/clean_data/E-GEOD-29378_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_36980<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-36980/clean_data/E-GEOD-36980_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_48350<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-48350/Pre-Processing/clean_data/E-GEOD-48350_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_1297<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-1297/Pre-Processing/clean_data/E-GEOD-1297_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
E_GEOD_5281<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-5281/Pre-Processing/clean_data/E-GEOD-5281_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
inhouse_data<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/inhouse_data/clean_data/inhouse_brain_exprs_and_pheno2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MSBB_U133A<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MSBB/Pre-Processing/AffymetrixHG-U133A/clean_data/AMP_MSBB_U133A_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MSBB_U133B<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MSBB/Pre-Processing/AffymetrixHG-U133B/Pre-Processing/clean_data/AMP_MSBB_U133B_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MAYOeGWAS_Temporal_Cortex<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MAyoeGWAS/Pre-processing/clean_data/AMP_MayoeGWAS_Temporal_Cortex_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)
AMP_MAYOeGWAS_Cerebellum<-read.table("/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/AMP_MAyoeGWAS/Pre-processing/clean_data/AMP_MayoeGWAS_Cerebellum_pre-processed_data2.txt", sep="\t", head=T, check.names = F, stringsAsFactors=FALSE)

##### CHECK DATA FRAMES #####

#check tables read in OK

head(E_GEOD_28146[1:10])
head(E_GEOD_29378[1:10])
head(E_GEOD_36980[1:10])
head(E_GEOD_48350[1:10])
head(E_GEOD_1297[1:10])
head(E_GEOD_5281[1:10])
head(inhouse_data[1:10])
head(AMP_MSBB_U133A[1:10])
head(AMP_MSBB_U133B[1:10])
head(AMP_MAYOeGWAS_Temporal_Cortex[1:10])
head(AMP_MAYOeGWAS_Cerebellum[1:10])

# check colnames same in 1st 9 columns

colnames(E_GEOD_28146[1:9])
colnames(E_GEOD_29378[1:9])
colnames(E_GEOD_36980[1:9])
colnames(E_GEOD_48350[1:9])
colnames(E_GEOD_1297[1:9])
colnames(E_GEOD_5281[1:9])
colnames(inhouse_data[1:9])
colnames(AMP_MSBB_U133A[1:9])
colnames(AMP_MSBB_U133B[1:9])
colnames(AMP_MAYOeGWAS_Temporal_Cortex[1:9])
colnames(AMP_MAYOeGWAS_Cerebellum[1:9])

# check number of genes in each dataset

dim(E_GEOD_28146)
dim(E_GEOD_29378)
dim(E_GEOD_36980)
dim(E_GEOD_48350)
dim(E_GEOD_1297)
dim(E_GEOD_5281)
dim(inhouse_data)
dim(AMP_MSBB_U133A)
dim(AMP_MSBB_U133B)
dim(AMP_MAYOeGWAS_Temporal_Cortex)
dim(AMP_MAYOeGWAS_Cerebellum)

# check if any na exist 

any(is.na(E_GEOD_28146))=="TRUE"
any(is.na(E_GEOD_29378))=="TRUE"
any(is.na(E_GEOD_36980))=="TRUE"
any(is.na(E_GEOD_48350))=="TRUE"
any(is.na(E_GEOD_1297))=="TRUE"
any(is.na(E_GEOD_5281))=="TRUE"
any(is.na(inhouse_data))=="TRUE"
any(is.na(AMP_MSBB_U133A))=="TRUE"
any(is.na(AMP_MSBB_U133B))=="TRUE"
any(is.na(AMP_MAYOeGWAS_Temporal_Cortex))=="TRUE"
any(is.na(AMP_MAYOeGWAS_Cerebellum))=="TRUE"

# check diagnosis same

table(E_GEOD_28146$Diagnosis, exclude=NULL)
table(E_GEOD_29378$Diagnosis, exclude=NULL)
table(E_GEOD_36980$Diagnosis, exclude=NULL)
table(E_GEOD_48350$Diagnosis, exclude=NULL)
table(E_GEOD_1297$Diagnosis, exclude=NULL)
table(E_GEOD_5281$Diagnosis, exclude=NULL)
table(inhouse_data$Diagnosis, exclude=NULL)
table(AMP_MSBB_U133A$Diagnosis, exclude=NULL)
table(AMP_MSBB_U133B$Diagnosis, exclude=NULL)
table(AMP_MAYOeGWAS_Temporal_Cortex$Diagnosis, exclude=NULL)
table(AMP_MAYOeGWAS_Cerebellum$Diagnosis, exclude=NULL)

#check Tissue

table(E_GEOD_28146$Tissue, exclude=NULL)
table(E_GEOD_29378$Tissue, exclude=NULL)
table(E_GEOD_36980$Tissue, exclude=NULL)
table(E_GEOD_48350$Tissue, exclude=NULL)
table(E_GEOD_1297$Tissue, exclude=NULL)
table(E_GEOD_5281$Tissue, exclude=NULL)
table(inhouse_data$Tissue, exclude=NULL)
table(AMP_MSBB_U133A$Tissue, exclude=NULL)
table(AMP_MSBB_U133B$Tissue, exclude=NULL)
table(AMP_MAYOeGWAS_Temporal_Cortex$Tissue, exclude=NULL)
table(AMP_MAYOeGWAS_Cerebellum$Tissue, exclude=NULL)

# change E_GEOD_36980 Tissue
E_GEOD_36980[E_GEOD_36980$Tissue=="frontal cortex",1]<-"Frontal_Cortex"
E_GEOD_36980[E_GEOD_36980$Tissue=="hippocampus",1]<-"Hippocampus"
E_GEOD_36980[E_GEOD_36980$Tissue=="Temporal cortex",1]<-"Temporal_Cortex"

table(E_GEOD_36980$Tissue, exclude=NULL)

# table of diagnosis and BRAAK score -
# Control - 0-1
# AD - 3,4,5,6
# remove AD + BRAAK<=3

names(E_GEOD_28146[1:10])

table(E_GEOD_28146[c(5,7)], exclude=NULL)
table(E_GEOD_36980[c(5,7)], exclude=NULL)
table(E_GEOD_1297[c(5,7)], exclude=NULL)
table(E_GEOD_5281[c(5,7)], exclude=NULL)
table(AMP_MAYOeGWAS_Temporal_Cortex[c(5,7)], exclude=NULL)
table(AMP_MAYOeGWAS_Cerebellum[c(5,7)], exclude=NULL)
table(E_GEOD_29378[c(5,7)], exclude=NULL)


# remove AD + BRAAK<=3
table(E_GEOD_48350[c(5,7)], exclude=NULL)

E_GEOD_48350_clean<-subset(E_GEOD_48350, BRAAK!="1" & BRAAK!="2" & BRAAK!="3")

table(E_GEOD_48350_clean[c(5,7)], exclude=NULL)


table(E_GEOD_48350[c(5,7)], exclude=NULL)
table(E_GEOD_48350[c(5,5)], exclude=NULL)

# remove AD + BRAAK<3 - too mandy AD samples with BRAAK as 0 - possibly unknown
table(inhouse_data[c(5,7)], exclude=NULL)

# remove controls where BRAAK is 3 and 5 | AD where BRAAK=3 - leave AD + BRAAK=0 - posibly unknown

inhouse_to_remove<-subset(inhouse_data, Diagnosis=="Control" & BRAAK=="5" | Diagnosis=="Control" & BRAAK=="3" | Diagnosis=="AD" & BRAAK=="3")[1:10]
inhouse_to_remove

inhouse_data_clean<-inhouse_data[!(rownames(inhouse_data))%in%(rownames(inhouse_to_remove)),]
table(inhouse_data_clean[c(5,7)], exclude=NULL)

# AMP_MSBB_U133A - decide what braak score for control, AD and incipient 
# 0-1-2 control
# 3 incipient
# 4 moderate
# 5-6 severe

table(AMP_MSBB_U133A[c(5,7)], exclude=NULL)
table(AMP_MSBB_U133A[c(5,8)], exclude=NULL)
table(AMP_MSBB_U133A[c(5,4)], exclude=NULL)

AMP_MSBB_U133A_clean<-subset(AMP_MSBB_U133A, Diagnosis=="AD" & BRAAK %in% c("4", "5", "6") | Diagnosis=="Control" & BRAAK %in% c("0", "1", "2"))

table(AMP_MSBB_U133A_clean[c(5,7)], exclude=NULL)
table(AMP_MSBB_U133A_clean[c(5,8)], exclude=NULL)
table(AMP_MSBB_U133A_clean[c(5,4)], exclude=NULL)


dim(AMP_MSBB_U133A_clean)

# AMP_MSBB_U133A - decide what braak score for control, AD and incipient 
# 0-1-2 control
# 3 incipient
# 4 moderate
# 5-6 severe

table(AMP_MSBB_U133B[c(5,7)], exclude=NULL)
table(AMP_MSBB_U133B[c(5,8)], exclude=NULL)
table(AMP_MSBB_U133B[c(5,4)], exclude=NULL)

AMP_MSBB_U133B_clean<-subset(AMP_MSBB_U133B, Diagnosis=="AD" & BRAAK %in% c("4", "5", "6") | Diagnosis=="Control" & BRAAK %in% c("0", "1", "2"))
table(AMP_MSBB_U133B_clean[c(5,7)], exclude=NULL)
table(AMP_MSBB_U133B_clean[c(5,8)], exclude=NULL)
table(AMP_MSBB_U133B_clean[c(5,4)], exclude=NULL)

# final clean data sets - adjust name - add "_clean to name"

E_GEOD_28146_clean<-E_GEOD_28146
E_GEOD_36980_clean<-E_GEOD_36980
E_GEOD_1297_clean<-E_GEOD_1297
E_GEOD_5281_clean<-E_GEOD_5281
AMP_MAYOeGWAS_Temporal_Cortex_clean<-AMP_MAYOeGWAS_Temporal_Cortex
AMP_MAYOeGWAS_Cerebellum_clean<-AMP_MAYOeGWAS_Cerebellum
E_GEOD_29378_clean<-E_GEOD_29378

# check diag sub cat

table(E_GEOD_28146_clean$Diagnosis_sub_cat, exclude=NULL)
table(E_GEOD_29378_clean$Diagnosis_sub_cat, exclude=NULL)
table(E_GEOD_36980_clean$Diagnosis_sub_cat, exclude=NULL)
table(E_GEOD_48350_clean$Diagnosis_sub_cat, exclude=NULL)
table(E_GEOD_1297_clean$Diagnosis_sub_cat, exclude=NULL)
table(E_GEOD_5281_clean$Diagnosis_sub_cat, exclude=NULL)
table(inhouse_data_clean$Diagnosis_sub_cat, exclude=NULL)
table(AMP_MSBB_U133A_clean$Diagnosis_sub_cat, exclude=NULL)
table(AMP_MSBB_U133B_clean$Diagnosis_sub_cat, exclude=NULL)
table(AMP_MAYOeGWAS_Temporal_Cortex_clean$Diagnosis_sub_cat, exclude=NULL)
table(AMP_MAYOeGWAS_Cerebellum_clean$Diagnosis_sub_cat, exclude=NULL)

##### CHECK BEFORE AFTER SAMPLE REMOVAL #####

dim(E_GEOD_28146)
dim(E_GEOD_29378)
dim(E_GEOD_36980)
dim(E_GEOD_48350)
dim(E_GEOD_1297)
dim(E_GEOD_5281)
dim(inhouse_data)
dim(AMP_MSBB_U133A)
dim(AMP_MSBB_U133B)
dim(AMP_MAYOeGWAS_Temporal_Cortex)
dim(AMP_MAYOeGWAS_Cerebellum)

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

##### DIFFERENTIAL EXPRESSION ON EACH DATASET BY CASE/CONTROL PER TISSUE FUNCTION - #####

#create function to run differential expression on each dataset

run_diff_exprs_analysis<-function(dataset, x){
  #sort dataset by Diagnosis column - want AD samples 1st
  dataset_sorted<-dataset[order(dataset$Diagnosis),]
  #split dataset into expres and diagnosis
  dataset_exprs<-dataset_sorted[10:dim(dataset_sorted)[2]]
  dataset_pheno<-dataset_sorted[,5]
  #replace sample name to numbers
  rownames(dataset_exprs)<-c(1:dim(dataset)[1])
  # setup experimental desgin
  design <- model.matrix(~0 + dataset_pheno)
  # change colnames to AD and Control
  colnames(design)<-c("AD", "Control")
  # transpose dataset, convert to numeric 
  transposed_dataset_exprs<-t(dataset_exprs)
  #run diff expression
  dataset_exprs_fit <- lmFit(transposed_dataset_exprs, design, method="robust", maxit=x)
  dataset_exprs_contrast_matrix<- makeContrasts(AD-Control, levels=design)
  dataset_exprs_contrast_fit <-contrasts.fit(dataset_exprs_fit, dataset_exprs_contrast_matrix)
  dataset_exprs_contrast_ebayes <- eBayes(dataset_exprs_contrast_fit, robust=T)
  dataset_exprs_top_genes <- topTable(dataset_exprs_contrast_ebayes, number=(dim(transposed_dataset_exprs)[1]), coef=1, adjust.method="fdr", confint=TRUE) 
  return(dataset_exprs_top_genes)
}

# apply differential expression on each dataset - need to subset each dataset by tissue - i.e have only same tissue type in each dataset

#####  E_GEOD_28146 DE #####

table(E_GEOD_28146_clean$Tissue, exclude=NULL)
table(E_GEOD_28146_clean$Diagnosis, exclude=NULL)

E_GEOD_28146_top_genes_Hippocampus_CA1<-run_diff_exprs_analysis(E_GEOD_28146_clean, 100)
E_GEOD_28146_sig_genes_Hippocampus_CA1<-subset(E_GEOD_28146_top_genes_Hippocampus_CA1, adj.P.Val <= 0.05)
head(E_GEOD_28146_sig_genes_Hippocampus_CA1)
dim(E_GEOD_28146_clean)
dim(E_GEOD_28146_sig_genes_Hippocampus_CA1)

#####  E_GEOD_29378 DE ####

table(E_GEOD_29378_clean$Tissue, exclude=NULL)
table(E_GEOD_29378_clean$Diagnosis, exclude=NULL)

#split dataset by tissue

head(E_GEOD_29378_clean)[1:5]

E_GEOD_29378_clean_CA1<-E_GEOD_29378_clean[E_GEOD_29378_clean$Tissue=="Hippocampus_CA1",]
E_GEOD_29378_clean_CA3<-E_GEOD_29378_clean[E_GEOD_29378_clean$Tissue=="Hippocampus_CA3",]

table(E_GEOD_29378_clean_CA1$Diagnosis)
table(E_GEOD_29378_clean_CA3$Diagnosis)
dim(E_GEOD_29378_clean_CA3)

E_GEOD_29378_top_genes_CA1<-run_diff_exprs_analysis(E_GEOD_29378_clean_CA1, 100)
E_GEOD_29378_sig_genes_CA1<-subset(E_GEOD_29378_top_genes_CA1, adj.P.Val <= 0.05)
head(E_GEOD_29378_sig_genes_CA1)
dim(E_GEOD_29378_sig_genes_CA1)

E_GEOD_29378_top_genes_CA3<-run_diff_exprs_analysis(E_GEOD_29378_clean_CA3, 100)
E_GEOD_29378_sig_genes_CA3<-subset(E_GEOD_29378_top_genes_CA3, adj.P.Val <= 0.05)
head(E_GEOD_29378_sig_genes_CA3)
dim(E_GEOD_29378_sig_genes_CA3)

#####  E_GEOD_36980 DE #####

table(E_GEOD_36980_clean$Tissue, exclude=NULL)
table(E_GEOD_36980_clean$Diagnosis, exclude=NULL)

#split dataset by tissue
dim(E_GEOD_36980_clean)
head(E_GEOD_36980_clean)[1:5]

E_GEOD_36980_clean_Frontal_Cortex<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Frontal_Cortex",]
E_GEOD_36980_clean_Hippocampus<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Hippocampus",]
E_GEOD_36980_clean_Temporal_Cortex<-E_GEOD_36980_clean[E_GEOD_36980_clean$Tissue=="Temporal_Cortex",]

table(E_GEOD_36980_clean_Frontal_Cortex$Tissue, exclude=NULL)
table(E_GEOD_36980_clean_Frontal_Cortex$Diagnosis, exclude=NULL)
table(E_GEOD_36980_clean_Hippocampus$Tissue, exclude=NULL)
table(E_GEOD_36980_clean_Hippocampus$Diagnosis, exclude=NULL)
table(E_GEOD_36980_clean_Temporal_Cortex$Tissue, exclude=NULL)
table(E_GEOD_36980_clean_Temporal_Cortex$Diagnosis, exclude=NULL)

E_GEOD_36980_top_genes_Frontal_Cortex<-run_diff_exprs_analysis(E_GEOD_36980_clean_Frontal_Cortex, 100)
E_GEOD_36980_sig_genes_Frontal_Cortex<-subset(E_GEOD_36980_top_genes_Frontal_Cortex, adj.P.Val <= 0.05)
head(E_GEOD_36980_sig_genes_Frontal_Cortex)
dim(E_GEOD_36980_sig_genes_Frontal_Cortex)

E_GEOD_36980_top_genes_Hippocampus<-run_diff_exprs_analysis(E_GEOD_36980_clean_Hippocampus, 1000)
E_GEOD_36980_sig_genes_Hippocampus<-subset(E_GEOD_36980_top_genes_Hippocampus, adj.P.Val <= 0.05)
head(E_GEOD_36980_sig_genes_Hippocampus)
dim(E_GEOD_36980_sig_genes_Hippocampus)

E_GEOD_36980_top_genes_Temporal_Cortex<-run_diff_exprs_analysis(E_GEOD_36980_clean_Temporal_Cortex, 100)
E_GEOD_36980_sig_genes_Temporal_Cortex<-subset(E_GEOD_36980_top_genes_Temporal_Cortex, adj.P.Val <= 0.05)
head(E_GEOD_36980_sig_genes_Temporal_Cortex)
dim(E_GEOD_36980_sig_genes_Temporal_Cortex)

#####  E_GEOD_48350 DE #####

table(E_GEOD_48350_clean$Tissue, exclude=NULL)
table(E_GEOD_48350_clean$Diagnosis, exclude=NULL)
dim(E_GEOD_48350)

#split dataset by tissue

head(E_GEOD_48350_clean)[1:5]

E_GEOD_48350_clean_Entorhinal_Cortex<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Entorhinal_Cortex",]
E_GEOD_48350_clean_Hippocampus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Hippocampus",]
E_GEOD_48350_clean_Postcentral_Gyrus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Postcentral_Gyrus",]
E_GEOD_48350_clean_Superior_Frontal_Gyrus<-E_GEOD_48350_clean[E_GEOD_48350_clean$Tissue=="Superior_Frontal_Gyrus",]

table(E_GEOD_48350_clean_Entorhinal_Cortex$Tissue, exclude=NULL)
table(E_GEOD_48350_clean_Entorhinal_Cortex$Diagnosis, exclude=NULL)
table(E_GEOD_48350_clean_Hippocampus$Tissue, exclude=NULL)
table(E_GEOD_48350_clean_Hippocampus$Diagnosis, exclude=NULL)
table(E_GEOD_48350_clean_Postcentral_Gyrus$Tissue, exclude=NULL)
table(E_GEOD_48350_clean_Postcentral_Gyrus$Diagnosis, exclude=NULL)
table(E_GEOD_48350_clean_Superior_Frontal_Gyrus$Tissue, exclude=NULL)
table(E_GEOD_48350_clean_Superior_Frontal_Gyrus$Diagnosis, exclude=NULL)

E_GEOD_48350_top_genes_Entorhinal_Cortex<-run_diff_exprs_analysis(E_GEOD_48350_clean_Entorhinal_Cortex, 100)
E_GEOD_48350_sig_genes_Entorhinal_Cortex<-subset(E_GEOD_48350_top_genes_Entorhinal_Cortex, adj.P.Val <= 0.05)
head(E_GEOD_48350_sig_genes_Entorhinal_Cortex)
dim(E_GEOD_48350_sig_genes_Entorhinal_Cortex)

E_GEOD_48350_top_genes_Hippocampus<-run_diff_exprs_analysis(E_GEOD_48350_clean_Hippocampus, 100)
E_GEOD_48350_sig_genes_Hippocampus<-subset(E_GEOD_48350_top_genes_Hippocampus, adj.P.Val <= 0.05)
head(E_GEOD_48350_sig_genes_Hippocampus)
dim(E_GEOD_48350_sig_genes_Hippocampus)

E_GEOD_48350_top_genes_Postcentral_Gyrus<-run_diff_exprs_analysis(E_GEOD_48350_clean_Postcentral_Gyrus, 100)
E_GEOD_48350_sig_genes_Postcentral_Gyrus<-subset(E_GEOD_48350_top_genes_Postcentral_Gyrus, adj.P.Val <= 0.05)
head(E_GEOD_48350_sig_genes_Postcentral_Gyrus)
dim(E_GEOD_48350_sig_genes_Postcentral_Gyrus)

E_GEOD_48350_top_genes_Superior_Frontal_Gyrus<-run_diff_exprs_analysis(E_GEOD_48350_clean_Superior_Frontal_Gyrus, 100)
E_GEOD_48350_sig_genes_Superior_Frontal_Gyrus<-subset(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus, adj.P.Val <= 0.05)
head(E_GEOD_48350_sig_genes_Superior_Frontal_Gyrus)
dim(E_GEOD_48350_sig_genes_Superior_Frontal_Gyrus)

#####  E_GEOD_1297 DE #####

table(E_GEOD_1297_clean$Tissue, exclude=NULL)
table(E_GEOD_1297_clean$Diagnosis, exclude=NULL)

E_GEOD_1297_top_genes_Hippocampus<-run_diff_exprs_analysis(E_GEOD_1297_clean, 1000)
E_GEOD_1297_sig_genes_Hippocampus<-subset(E_GEOD_1297_top_genes_Hippocampus, adj.P.Val <= 0.05)
head(E_GEOD_1297_sig_genes_Hippocampus)
dim(E_GEOD_1297_sig_genes_Hippocampus)

#####  E_GEOD_5281 DE #####

table(E_GEOD_5281_clean$Tissue, exclude=NULL)
table(E_GEOD_5281_clean$Diagnosis, exclude=NULL)

#split dataset by tissue

head(E_GEOD_5281_clean)[1:5]

E_GEOD_5281_clean_Entorhinal_Cortex<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Entorhinal_Cortex",]
E_GEOD_5281_clean_Hippocampus_CA1<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Hippocampus_CA1",]
E_GEOD_5281_clean_Medial_Temporal_Gyrus<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Medial_Temporal_Gyrus",]
E_GEOD_5281_clean_Superior_Frontal_Gyrus<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Superior_Frontal_Gyrus",]
E_GEOD_5281_clean_Posterior_Singulate<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Posterior_Singulate",]
E_GEOD_5281_clean_Primary_Visual_Cortex<-E_GEOD_5281_clean[E_GEOD_5281_clean$Tissue=="Primary_Visual_Cortex",]


table(E_GEOD_5281_clean_Entorhinal_Cortex$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Entorhinal_Cortex$Diagnosis, exclude=NULL)
table(E_GEOD_5281_clean_Hippocampus_CA1$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Hippocampus_CA1$Diagnosis, exclude=NULL)
table(E_GEOD_5281_clean_Medial_Temporal_Gyrus$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Medial_Temporal_Gyrus$Diagnosis, exclude=NULL)
table(E_GEOD_5281_clean_Superior_Frontal_Gyrus$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Superior_Frontal_Gyrus$Diagnosis, exclude=NULL)
table(E_GEOD_5281_clean_Posterior_Singulate$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Posterior_Singulate$Diagnosis, exclude=NULL)
table(E_GEOD_5281_clean_Primary_Visual_Cortex$Tissue, exclude=NULL)
table(E_GEOD_5281_clean_Primary_Visual_Cortex$Diagnosis, exclude=NULL)


E_GEOD_5281_top_genes_Entorhinal_Cortex<-run_diff_exprs_analysis(E_GEOD_5281_clean_Entorhinal_Cortex, 100)
E_GEOD_5281_sig_genes_Entorhinal_Cortex<-subset(E_GEOD_5281_top_genes_Entorhinal_Cortex, adj.P.Val <= 0.05)
head(E_GEOD_5281_sig_genes_Entorhinal_Cortex)
dim(E_GEOD_5281_sig_genes_Entorhinal_Cortex)

E_GEOD_5281_top_genes_Hippocampus_CA1<-run_diff_exprs_analysis(E_GEOD_5281_clean_Hippocampus_CA1, 100)
E_GEOD_5281_sig_genes_Hippocampus_CA1<-subset(E_GEOD_5281_top_genes_Hippocampus_CA1, adj.P.Val <= 0.05)
head(E_GEOD_5281_sig_genes_Hippocampus_CA1)
dim(E_GEOD_5281_sig_genes_Hippocampus_CA1)

E_GEOD_5281_top_genes_Medial_Temporal_Gyrus<-run_diff_exprs_analysis(E_GEOD_5281_clean_Medial_Temporal_Gyrus, 1000)
E_GEOD_5281_sig_genes_Medial_Temporal_Gyrus<-subset(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus, adj.P.Val <= 0.05)
head(E_GEOD_5281_sig_genes_Medial_Temporal_Gyrus)
dim(E_GEOD_5281_sig_genes_Medial_Temporal_Gyrus)

E_GEOD_5281_top_genes_Superior_Frontal_Gyrus<-run_diff_exprs_analysis(E_GEOD_5281_clean_Superior_Frontal_Gyrus, 100)
E_GEOD_5281_sig_genes_Superior_Frontal_Gyrus<-subset(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus, adj.P.Val <= 0.05)
head(E_GEOD_5281_sig_genes_Superior_Frontal_Gyrus)
dim(E_GEOD_5281_sig_genes_Superior_Frontal_Gyrus)

E_GEOD_5281_top_genes_Posterior_Singulate<-run_diff_exprs_analysis(E_GEOD_5281_clean_Posterior_Singulate, 1000)
E_GEOD_5281_sig_genes_Posterior_Singulate<-subset(E_GEOD_5281_top_genes_Posterior_Singulate, adj.P.Val <= 0.05)
head(E_GEOD_5281_sig_genes_Posterior_Singulate)
dim(E_GEOD_5281_sig_genes_Posterior_Singulate)

E_GEOD_5281_top_genes_Primary_Visual_Cortex<-run_diff_exprs_analysis(E_GEOD_5281_clean_Primary_Visual_Cortex, 100)
E_GEOD_5281_sig_genes_Primary_Visual_Cortex<-subset(E_GEOD_5281_top_genes_Primary_Visual_Cortex, adj.P.Val <= 0.05)
dim(E_GEOD_5281_sig_genes_Primary_Visual_Cortex)
head(E_GEOD_5281_sig_genes_Primary_Visual_Cortex)

##### inhouse data DE ####

# remove CO_AD diagnosis

inhouse_data_clean<-subset(inhouse_data_clean, Diagnosis!="CO_AD")

table(inhouse_data_clean$Tissue, exclude=NULL)
table(inhouse_data_clean$Diagnosis, exclude=NULL)

#split dataset by tissue

head(inhouse_data_clean)[1:5]

inhouse_data_clean_Entorhinal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Entorhinal_Cortex",]
inhouse_data_clean_Cerebellum<-inhouse_data_clean[inhouse_data_clean$Tissue=="Cerebellum",]
inhouse_data_clean_Frontal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Frontal_Cortex",]
inhouse_data_clean_Temporal_Cortex<-inhouse_data_clean[inhouse_data_clean$Tissue=="Temporal_Cortex",]


table(inhouse_data_clean_Entorhinal_Cortex$Tissue, exclude=NULL)
table(inhouse_data_clean_Entorhinal_Cortex$Diagnosis, exclude=NULL)
table(inhouse_data_clean_Cerebellum$Tissue, exclude=NULL)
table(inhouse_data_clean_Cerebellum$Diagnosis, exclude=NULL)
table(inhouse_data_clean_Frontal_Cortex$Tissue, exclude=NULL)
table(inhouse_data_clean_Frontal_Cortex$Diagnosis, exclude=NULL)
table(inhouse_data_clean_Temporal_Cortex$Tissue, exclude=NULL)
table(inhouse_data_clean_Temporal_Cortex$Diagnosis, exclude=NULL)


inhouse_data_clean_top_genes_Entorhinal_Cortex<-run_diff_exprs_analysis(inhouse_data_clean_Entorhinal_Cortex, 100)
inhouse_data_clean_sig_genes_Entorhinal_Cortex<-subset(inhouse_data_clean_top_genes_Entorhinal_Cortex, adj.P.Val <= 0.05)
head(inhouse_data_clean_sig_genes_Entorhinal_Cortex)
dim(inhouse_data_clean_sig_genes_Entorhinal_Cortex)

inhouse_data_clean_top_genes_Cerebellum<-run_diff_exprs_analysis(inhouse_data_clean_Cerebellum, 100)
inhouse_data_clean_sig_genes_Cerebellum<-subset(inhouse_data_clean_top_genes_Cerebellum, adj.P.Val <= 0.05)
head(inhouse_data_clean_sig_genes_Cerebellum)
dim(inhouse_data_clean_sig_genes_Cerebellum)

inhouse_data_clean_top_genes_Frontal_Cortex<-run_diff_exprs_analysis(inhouse_data_clean_Frontal_Cortex, 100)
inhouse_data_clean_sig_genes_Frontal_Cortex<-subset(inhouse_data_clean_top_genes_Frontal_Cortex, adj.P.Val <= 0.05)
head(inhouse_data_clean_sig_genes_Frontal_Cortex)
dim(inhouse_data_clean_sig_genes_Frontal_Cortex)

inhouse_data_clean_top_genes_Temporal_Cortex<-run_diff_exprs_analysis(inhouse_data_clean_Temporal_Cortex, 100)
inhouse_data_clean_sig_genes_Temporal_Cortex<-subset(inhouse_data_clean_top_genes_Temporal_Cortex, adj.P.Val <= 0.05)
head(inhouse_data_clean_sig_genes_Temporal_Cortex)
dim(inhouse_data_clean_sig_genes_Temporal_Cortex)

##### AMP_MAYOeGWAS_Temporal_Cortex DE #####

table(AMP_MAYOeGWAS_Temporal_Cortex_clean$Tissue, exclude=NULL)
table(AMP_MAYOeGWAS_Temporal_Cortex_clean$Diagnosis, exclude=NULL)

AMP_MAYOeGWAS_top_genes_Temporal_Cortex<-run_diff_exprs_analysis(AMP_MAYOeGWAS_Temporal_Cortex_clean, 1000)
AMP_MAYOeGWAS_sig_genes_Temporal_Cortex<-subset(AMP_MAYOeGWAS_top_genes_Temporal_Cortex, adj.P.Val <= 0.05)
head(AMP_MAYOeGWAS_sig_genes_Temporal_Cortex)
dim(AMP_MAYOeGWAS_sig_genes_Temporal_Cortex)

##### AMP_MAYOeGWAS_Cerebellum DE #####

table(AMP_MAYOeGWAS_Cerebellum_clean$Tissue, exclude=NULL)
table(AMP_MAYOeGWAS_Cerebellum_clean$Diagnosis, exclude=NULL)

AMP_MAYOeGWAS_top_genes_Cerebellum<-run_diff_exprs_analysis(AMP_MAYOeGWAS_Cerebellum_clean, 1000)
AMP_MAYOeGWAS_sig_genes_Cerebellum<-subset(AMP_MAYOeGWAS_top_genes_Cerebellum, adj.P.Val <= 0.05)
head(AMP_MAYOeGWAS_sig_genes_Cerebellum)
dim(AMP_MAYOeGWAS_sig_genes_Cerebellum)

##### AMP_MSBB_U133A DE #####

table(AMP_MSBB_U133A_clean$Tissue, exclude=NULL)
table(AMP_MSBB_U133A_clean$Diagnosis, exclude=NULL)

#split dataset by tissue

head(AMP_MSBB_U133A_clean)[1:5]

AMP_MSBB_U133A_clean_Frontal_Pole<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Frontal_Pole",]
AMP_MSBB_U133A_clean_Precentral_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Precentral_Gyrus",]
AMP_MSBB_U133A_clean_Inferior_Frontal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Inferior_Frontal_Gyrus",]
AMP_MSBB_U133A_clean_Dorsolateral_Prefrontal_Cortex<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Dorsolateral_Prefrontal_Cortex",]
AMP_MSBB_U133A_clean_Superior_Parietal_Lobule<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Superior_Parietal_Lobule",]
AMP_MSBB_U133A_clean_Prefrontal_Cortex<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Prefrontal_Cortex",]
AMP_MSBB_U133A_clean_Parahippocampal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Parahippocampal_Gyrus",]
AMP_MSBB_U133A_clean_Hippocampus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Hippocampus",]
AMP_MSBB_U133A_clean_Inferior_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Inferior_Temporal_Gyrus",]
AMP_MSBB_U133A_clean_Middle_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Middle_Temporal_Gyrus",]
AMP_MSBB_U133A_clean_Superior_Temporal_Gyrus<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Superior_Temporal_Gyrus",]
AMP_MSBB_U133A_clean_Temporal_Pole<-AMP_MSBB_U133A_clean[AMP_MSBB_U133A_clean$Tissue=="Temporal_Pole",]

table(AMP_MSBB_U133A_clean_Frontal_Pole$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Precentral_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Inferior_Frontal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Dorsolateral_Prefrontal_Cortex$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Superior_Parietal_Lobule$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Prefrontal_Cortex$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Parahippocampal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Hippocampus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Inferior_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Middle_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Superior_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133A_clean_Temporal_Pole$Tissue , exclude=NULL)

table(AMP_MSBB_U133A_clean_Frontal_Pole$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Precentral_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Inferior_Frontal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Dorsolateral_Prefrontal_Cortex$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Superior_Parietal_Lobule$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Prefrontal_Cortex$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Parahippocampal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Hippocampus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Inferior_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Middle_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Superior_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133A_clean_Temporal_Pole$Diagnosis , exclude=NULL)

#run diff expr

AMP_MSBB_U133A_top_genes_Frontal_Pole<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Frontal_Pole, 100)
AMP_MSBB_U133A_sig_genes_Frontal_Pole<-subset(AMP_MSBB_U133A_top_genes_Frontal_Pole, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Frontal_Pole)
dim(AMP_MSBB_U133A_sig_genes_Frontal_Pole)

AMP_MSBB_U133A_top_genes_Precentral_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Precentral_Gyrus, 100)
AMP_MSBB_U133A_sig_genes_Precentral_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Precentral_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Precentral_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Precentral_Gyrus)

AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Inferior_Frontal_Gyrus, 100)
AMP_MSBB_U133A_sig_genes_Inferior_Frontal_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Inferior_Frontal_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Inferior_Frontal_Gyrus)

AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Dorsolateral_Prefrontal_Cortex, 100)
AMP_MSBB_U133A_sig_genes_Dorsolateral_Prefrontal_Cortex<-subset(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Dorsolateral_Prefrontal_Cortex)
dim(AMP_MSBB_U133A_sig_genes_Dorsolateral_Prefrontal_Cortex)

AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Superior_Parietal_Lobule, 1000)
AMP_MSBB_U133A_sig_genes_Superior_Parietal_Lobule<-subset(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Superior_Parietal_Lobule)
dim(AMP_MSBB_U133A_sig_genes_Superior_Parietal_Lobule)

AMP_MSBB_U133A_top_genes_Prefrontal_Cortex<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Prefrontal_Cortex, 1000)
AMP_MSBB_U133A_sig_genes_Prefrontal_Cortex<-subset(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Prefrontal_Cortex)
dim(AMP_MSBB_U133A_sig_genes_Prefrontal_Cortex)

AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Parahippocampal_Gyrus, 100)
AMP_MSBB_U133A_sig_genes_Parahippocampal_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Parahippocampal_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Parahippocampal_Gyrus)

AMP_MSBB_U133A_top_genes_Hippocampus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Hippocampus, 100)
AMP_MSBB_U133A_sig_genes_Hippocampus<-subset(AMP_MSBB_U133A_top_genes_Hippocampus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Hippocampus)
dim(AMP_MSBB_U133A_sig_genes_Hippocampus)

AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Inferior_Temporal_Gyrus, 100)
AMP_MSBB_U133A_sig_genes_Inferior_Temporal_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Inferior_Temporal_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Inferior_Temporal_Gyrus)

AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Middle_Temporal_Gyrus, 500)
AMP_MSBB_U133A_sig_genes_Middle_Temporal_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Middle_Temporal_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Middle_Temporal_Gyrus)

AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Superior_Temporal_Gyrus, 100)
AMP_MSBB_U133A_sig_genes_Superior_Temporal_Gyrus<-subset(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Superior_Temporal_Gyrus)
dim(AMP_MSBB_U133A_sig_genes_Superior_Temporal_Gyrus)

AMP_MSBB_U133A_top_genes_Temporal_Pole<-run_diff_exprs_analysis(AMP_MSBB_U133A_clean_Temporal_Pole, 100)
AMP_MSBB_U133A_sig_genes_Temporal_Pole<-subset(AMP_MSBB_U133A_top_genes_Temporal_Pole, adj.P.Val <= 0.05)
head(AMP_MSBB_U133A_sig_genes_Temporal_Pole)
dim(AMP_MSBB_U133A_sig_genes_Temporal_Pole)

##### AMP_MSBB_U133B DE #####

table(AMP_MSBB_U133B_clean$Tissue, exclude=NULL)

table(AMP_MSBB_U133B_clean$Diagnosis, exclude=NULL)

#split dataset by tissue

head(AMP_MSBB_U133B_clean)[1:5]

AMP_MSBB_U133B_clean_Frontal_Pole<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Frontal_Pole",]
AMP_MSBB_U133B_clean_Precentral_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Precentral_Gyrus",]
AMP_MSBB_U133B_clean_Inferior_Frontal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Inferior_Frontal_Gyrus",]
AMP_MSBB_U133B_clean_Dorsolateral_Prefrontal_Cortex<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Dorsolateral_Prefrontal_Cortex",]
AMP_MSBB_U133B_clean_Superior_Parietal_Lobule<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Superior_Parietal_Lobule",]
AMP_MSBB_U133B_clean_Prefrontal_Cortex<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Prefrontal_Cortex",]
AMP_MSBB_U133B_clean_Parahippocampal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Parahippocampal_Gyrus",]
AMP_MSBB_U133B_clean_Hippocampus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Hippocampus",]
AMP_MSBB_U133B_clean_Inferior_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Inferior_Temporal_Gyrus",]
AMP_MSBB_U133B_clean_Middle_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Middle_Temporal_Gyrus",]
AMP_MSBB_U133B_clean_Superior_Temporal_Gyrus<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Superior_Temporal_Gyrus",]
AMP_MSBB_U133B_clean_Temporal_Pole<-AMP_MSBB_U133B_clean[AMP_MSBB_U133B_clean$Tissue=="Temporal_Pole",]

table(AMP_MSBB_U133B_clean_Frontal_Pole$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Precentral_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Inferior_Frontal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Dorsolateral_Prefrontal_Cortex$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Superior_Parietal_Lobule$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Prefrontal_Cortex$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Parahippocampal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Hippocampus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Inferior_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Middle_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Superior_Temporal_Gyrus$Tissue , exclude=NULL)
table(AMP_MSBB_U133B_clean_Temporal_Pole$Tissue , exclude=NULL)

table(AMP_MSBB_U133B_clean_Frontal_Pole$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Precentral_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Inferior_Frontal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Dorsolateral_Prefrontal_Cortex$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Superior_Parietal_Lobule$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Prefrontal_Cortex$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Parahippocampal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Hippocampus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Inferior_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Middle_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Superior_Temporal_Gyrus$Diagnosis , exclude=NULL)
table(AMP_MSBB_U133B_clean_Temporal_Pole$Diagnosis , exclude=NULL)

#run diff expr

AMP_MSBB_U133B_top_genes_Frontal_Pole<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Frontal_Pole, 100)
AMP_MSBB_U133B_sig_genes_Frontal_Pole<-subset(AMP_MSBB_U133B_top_genes_Frontal_Pole, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Frontal_Pole)
dim(AMP_MSBB_U133B_sig_genes_Frontal_Pole)

AMP_MSBB_U133B_top_genes_Precentral_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Precentral_Gyrus, 500)
AMP_MSBB_U133B_sig_genes_Precentral_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Precentral_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Precentral_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Precentral_Gyrus)

AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Inferior_Frontal_Gyrus, 100)
AMP_MSBB_U133B_sig_genes_Inferior_Frontal_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Inferior_Frontal_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Inferior_Frontal_Gyrus)

AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Dorsolateral_Prefrontal_Cortex, 100)
AMP_MSBB_U133B_sig_genes_Dorsolateral_Prefrontal_Cortex<-subset(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Dorsolateral_Prefrontal_Cortex)
dim(AMP_MSBB_U133B_sig_genes_Dorsolateral_Prefrontal_Cortex)

AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Superior_Parietal_Lobule, 1000)
AMP_MSBB_U133B_sig_genes_Superior_Parietal_Lobule<-subset(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Superior_Parietal_Lobule)
dim(AMP_MSBB_U133B_sig_genes_Superior_Parietal_Lobule)

AMP_MSBB_U133B_top_genes_Prefrontal_Cortex<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Prefrontal_Cortex, 1000)
AMP_MSBB_U133B_sig_genes_Prefrontal_Cortex<-subset(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Prefrontal_Cortex)
dim(AMP_MSBB_U133B_sig_genes_Prefrontal_Cortex)

AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Parahippocampal_Gyrus, 100)
AMP_MSBB_U133B_sig_genes_Parahippocampal_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Parahippocampal_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Parahippocampal_Gyrus)

AMP_MSBB_U133B_top_genes_Hippocampus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Hippocampus, 100)
AMP_MSBB_U133B_sig_genes_Hippocampus<-subset(AMP_MSBB_U133B_top_genes_Hippocampus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Hippocampus)
dim(AMP_MSBB_U133B_sig_genes_Hippocampus)

AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Inferior_Temporal_Gyrus, 100)
AMP_MSBB_U133B_sig_genes_Inferior_Temporal_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Inferior_Temporal_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Inferior_Temporal_Gyrus)

AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Middle_Temporal_Gyrus, 100)
AMP_MSBB_U133B_sig_genes_Middle_Temporal_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Middle_Temporal_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Middle_Temporal_Gyrus)

AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Superior_Temporal_Gyrus, 100)
AMP_MSBB_U133B_sig_genes_Superior_Temporal_Gyrus<-subset(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Superior_Temporal_Gyrus)
dim(AMP_MSBB_U133B_sig_genes_Superior_Temporal_Gyrus)

AMP_MSBB_U133B_top_genes_Temporal_Pole<-run_diff_exprs_analysis(AMP_MSBB_U133B_clean_Temporal_Pole, 1000)
AMP_MSBB_U133B_sig_genes_Temporal_Pole<-subset(AMP_MSBB_U133B_top_genes_Temporal_Pole, adj.P.Val <= 0.05)
head(AMP_MSBB_U133B_sig_genes_Temporal_Pole)
dim(AMP_MSBB_U133B_sig_genes_Temporal_Pole)

##### MERGE ALL DIFF RESULTS INTO 1 TABLE #####

# list diff tables - 47 dataframes

head(E_GEOD_28146_top_genes_Hippocampus_CA1)
head(E_GEOD_29378_top_genes_CA1)
head(E_GEOD_29378_top_genes_CA3)
head(E_GEOD_36980_top_genes_Frontal_Cortex)
head(E_GEOD_36980_top_genes_Hippocampus)
head(E_GEOD_36980_top_genes_Temporal_Cortex)
head(E_GEOD_48350_top_genes_Entorhinal_Cortex)
head(E_GEOD_48350_top_genes_Hippocampus)
head(E_GEOD_48350_top_genes_Postcentral_Gyrus)
head(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus)
head(E_GEOD_1297_top_genes_Hippocampus)
head(E_GEOD_5281_top_genes_Entorhinal_Cortex)
head(E_GEOD_5281_top_genes_Hippocampus_CA1)
head(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus)
head(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus)
head(E_GEOD_5281_top_genes_Posterior_Singulate)
head(E_GEOD_5281_top_genes_Primary_Visual_Cortex)
head(inhouse_data_clean_top_genes_Entorhinal_Cortex)
head(inhouse_data_clean_top_genes_Cerebellum)
head(inhouse_data_clean_top_genes_Frontal_Cortex)
head(inhouse_data_clean_top_genes_Temporal_Cortex)
head(AMP_MAYOeGWAS_top_genes_Temporal_Cortex)
head(AMP_MAYOeGWAS_top_genes_Cerebellum)
head(AMP_MSBB_U133A_top_genes_Frontal_Pole)
head(AMP_MSBB_U133A_top_genes_Precentral_Gyrus)
head(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus)
head(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex)
head(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule)
head(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex)
head(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus)
head(AMP_MSBB_U133A_top_genes_Hippocampus)
head(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus)
head(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus)
head(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus)
head(AMP_MSBB_U133A_top_genes_Temporal_Pole)
head(AMP_MSBB_U133B_top_genes_Frontal_Pole)
head(AMP_MSBB_U133B_top_genes_Precentral_Gyrus)
head(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus)
head(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex)
head(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule)
head(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex)
head(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus)
head(AMP_MSBB_U133B_top_genes_Hippocampus)
head(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus)
head(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus)
head(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus)
head(AMP_MSBB_U133B_top_genes_Temporal_Pole)

# add dataset ID and Tissue region into colnames

colnames(E_GEOD_28146_top_genes_Hippocampus_CA1)<-c("E_GEOD_28146_Hippocampus_CA1_logFC", "E_GEOD_28146_Hippocampus_CA1_CI.L", "E_GEOD_28146_Hippocampus_CA1_CI.R", "E_GEOD_28146_Hippocampus_CA1_AveExpr", "E_GEOD_28146_Hippocampus_CA1_t", "E_GEOD_28146_Hippocampus_CA1_P.Value",  "E_GEOD_28146_Hippocampus_CA1_adj.P.Value", "E_GEOD_28146_Hippocampus_CA1_B")
colnames(E_GEOD_29378_top_genes_CA1)<-c("E_GEOD_29378_Hippocampus_CA1_logFC", "E_GEOD_29378_Hippocampus_CA1_CI.L", "E_GEOD_29378_Hippocampus_CA1_CI.R", "E_GEOD_29378_Hippocampus_CA1_AveExpr", "E_GEOD_29378_Hippocampus_CA1_t", "E_GEOD_29378_Hippocampus_CA1_P.Value",  "E_GEOD_29378_Hippocampus_CA1_adj.P.Value", "E_GEOD_29378_Hippocampus_CA1_B")
colnames(E_GEOD_29378_top_genes_CA3)<-c("E_GEOD_29378_Hippocampus_CA3_logFC", "E_GEOD_29378_Hippocampus_CA3_CI.L", "E_GEOD_29378_Hippocampus_CA3_CI.R", "E_GEOD_29378_Hippocampus_CA3_AveExpr", "E_GEOD_29378_Hippocampus_CA3_t", "E_GEOD_29378_Hippocampus_CA3_P.Value",  "E_GEOD_29378_Hippocampus_CA3_adj.P.Value", "E_GEOD_29378_Hippocampus_CA3_B")
colnames(E_GEOD_36980_top_genes_Frontal_Cortex)<-c("E_GEOD_36980_Frontal_Cortex_logFC", "E_GEOD_36980_Frontal_Cortex_CI.L", "E_GEOD_36980_Frontal_Cortex_CI.R", "E_GEOD_36980_Frontal_Cortex_AveExpr", "E_GEOD_36980_Frontal_Cortex_t", "E_GEOD_36980_Frontal_Cortex_P.Value",  "E_GEOD_36980_Frontal_Cortex_adj.P.Value", "E_GEOD_36980_Frontal_Cortex_B")
colnames(E_GEOD_36980_top_genes_Hippocampus)<-c("E_GEOD_36980_Hippocampus_logFC", "E_GEOD_36980_Hippocampus_CI.L", "E_GEOD_36980_Hippocampus_CI.R", "E_GEOD_36980_Hippocampus_AveExpr", "E_GEOD_36980_Hippocampus_t", "E_GEOD_36980_Hippocampus_P.Value",  "E_GEOD_36980_Hippocampus_adj.P.Value", "E_GEOD_36980_Hippocampus_B")
colnames(E_GEOD_36980_top_genes_Temporal_Cortex)<-c("E_GEOD_36980_Temporal_Cortex_logFC", "E_GEOD_36980_Temporal_Cortex_CI.L", "E_GEOD_36980_Temporal_Cortex_CI.R", "E_GEOD_36980_Temporal_Cortex_AveExpr", "E_GEOD_36980_Temporal_Cortex_t", "E_GEOD_36980_Temporal_Cortex_P.Value",  "E_GEOD_36980_Temporal_Cortex_adj.P.Value", "E_GEOD_36980_Temporal_Cortex_B")
colnames(E_GEOD_48350_top_genes_Entorhinal_Cortex)<-c("E_GEOD_48350_Entorhinal_Cortex_logFC", "E_GEOD_48350_Entorhinal_Cortex_CI.L", "E_GEOD_48350_Entorhinal_Cortex_CI.R", "E_GEOD_48350_Entorhinal_Cortex_AveExpr", "E_GEOD_48350_Entorhinal_Cortex_t", "E_GEOD_48350_Entorhinal_Cortex_P.Value",  "E_GEOD_48350_Entorhinal_Cortex_adj.P.Value", "E_GEOD_48350_Entorhinal_Cortex_B")
colnames(E_GEOD_48350_top_genes_Hippocampus)<-c("E_GEOD_48350_Hippocampus_logFC", "E_GEOD_48350_Hippocampus_CI.L", "E_GEOD_48350_Hippocampus_CI.R", "E_GEOD_48350_Hippocampus_AveExpr", "E_GEOD_48350_Hippocampus_t", "E_GEOD_48350_Hippocampus_P.Value",  "E_GEOD_48350_Hippocampus_adj.P.Value", "E_GEOD_48350_Hippocampus_B")
colnames(E_GEOD_48350_top_genes_Postcentral_Gyrus)<-c("E_GEOD_48350_Postcentral_Gyrus_logFC", "E_GEOD_48350_Postcentral_Gyrus_CI.L", "E_GEOD_48350_Postcentral_Gyrus_CI.R", "E_GEOD_48350_Postcentral_Gyrus_AveExpr", "E_GEOD_48350_Postcentral_Gyrus_t", "E_GEOD_48350_Postcentral_Gyrus_P.Value",  "E_GEOD_48350_Postcentral_Gyrus_adj.P.Value", "E_GEOD_48350_Postcentral_Gyrus_B")
colnames(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus)<-c("E_GEOD_48350_Superior_Frontal_Gyrus_logFC", "E_GEOD_48350_Superior_Frontal_Gyrus_CI.L", "E_GEOD_48350_Superior_Frontal_Gyrus_CI.R", "E_GEOD_48350_Superior_Frontal_Gyrus_AveExpr", "E_GEOD_48350_Superior_Frontal_Gyrus_t", "E_GEOD_48350_Superior_Frontal_Gyrus_P.Value",  "E_GEOD_48350_Superior_Frontal_Gyrus_adj.P.Value", "E_GEOD_48350_Superior_Frontal_Gyrus_B")
colnames(E_GEOD_1297_top_genes_Hippocampus)<-c("E_GEOD_1297_Hippocampus_logFC", "E_GEOD_1297_Hippocampus_CI.L", "E_GEOD_1297_Hippocampus_CI.R", "E_GEOD_1297_Hippocampus_AveExpr", "E_GEOD_1297_Hippocampus_t", "E_GEOD_1297_Hippocampus_P.Value",  "E_GEOD_1297_Hippocampus_adj.P.Value", "E_GEOD_1297_Hippocampus_B")
colnames(E_GEOD_5281_top_genes_Entorhinal_Cortex)<-c("E_GEOD_5281_Entorhinal_Cortex_logFC", "E_GEOD_5281_Entorhinal_Cortex_CI.L", "E_GEOD_5281_Entorhinal_Cortex_CI.R", "E_GEOD_5281_Entorhinal_Cortex_AveExpr", "E_GEOD_5281_Entorhinal_Cortex_t", "E_GEOD_5281_Entorhinal_Cortex_P.Value",  "E_GEOD_5281_Entorhinal_Cortex_adj.P.Value", "E_GEOD_5281_Entorhinal_Cortex_B")
colnames(E_GEOD_5281_top_genes_Hippocampus_CA1)<-c("E_GEOD_5281_Hippocampus_CA1_logFC", "E_GEOD_5281_Hippocampus_CA1_CI.L", "E_GEOD_5281_Hippocampus_CA1_CI.R", "E_GEOD_5281_Hippocampus_CA1_AveExpr", "E_GEOD_5281_Hippocampus_CA1_t", "E_GEOD_5281_Hippocampus_CA1_P.Value",  "E_GEOD_5281_Hippocampus_CA1_adj.P.Value", "E_GEOD_5281_Hippocampus_CA1_B")
colnames(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus)<-c("E_GEOD_5281_Medial_Temporal_Gyrus_logFC", "E_GEOD_5281_Medial_Temporal_Gyrus_CI.L", "E_GEOD_5281_Medial_Temporal_Gyrus_CI.R", "E_GEOD_5281_Medial_Temporal_Gyrus_AveExpr", "E_GEOD_5281_Medial_Temporal_Gyrus_t", "E_GEOD_5281_Medial_Temporal_Gyrus_P.Value",  "E_GEOD_5281_Medial_Temporal_Gyrus_adj.P.Value", "E_GEOD_5281_Medial_Temporal_Gyrus_B")
colnames(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus)<-c("E_GEOD_5281_Superior_Frontal_Gyrus_logFC", "E_GEOD_5281_Superior_Frontal_Gyrus_CI.L", "E_GEOD_5281_Superior_Frontal_Gyrus_CI.R", "E_GEOD_5281_Superior_Frontal_Gyrus_AveExpr", "E_GEOD_5281_Superior_Frontal_Gyrus_t", "E_GEOD_5281_Superior_Frontal_Gyrus_P.Value",  "E_GEOD_5281_Superior_Frontal_Gyrus_adj.P.Value", "E_GEOD_5281_Superior_Frontal_Gyrus_B")
colnames(E_GEOD_5281_top_genes_Posterior_Singulate)<-c("E_GEOD_5281_Posterior_Singulate_logFC", "E_GEOD_5281_Posterior_Singulate_CI.L", "E_GEOD_5281_Posterior_Singulate_CI.R", "E_GEOD_5281_Posterior_Singulate_AveExpr", "E_GEOD_5281_Posterior_Singulate_t", "E_GEOD_5281_Posterior_Singulate_P.Value",  "E_GEOD_5281_Posterior_Singulate_adj.P.Value", "E_GEOD_5281_Posterior_Singulate_B")
colnames(E_GEOD_5281_top_genes_Primary_Visual_Cortex)<-c("E_GEOD_5281_Primary_Visual_Cortex_logFC", "E_GEOD_5281_Primary_Visual_Cortex_CI.L", "E_GEOD_5281_Primary_Visual_Cortex_CI.R", "E_GEOD_5281_Primary_Visual_Cortex_AveExpr", "E_GEOD_5281_Primary_Visual_Cortex_t", "E_GEOD_5281_Primary_Visual_Cortex_P.Value",  "E_GEOD_5281_Primary_Visual_Cortex_adj.P.Value", "E_GEOD_5281_Primary_Visual_Cortex_B")
colnames(inhouse_data_clean_top_genes_Entorhinal_Cortex)<-c("inhouse_Entorhinal_Cortex_logFC", "inhouse_Entorhinal_Cortex_CI.L", "inhouse_Entorhinal_Cortex_CI.R", "inhouse_Entorhinal_Cortex_AveExpr", "inhouse_Entorhinal_Cortex_t", "inhouse_Entorhinal_Cortex_P.Value",  "inhouse_Entorhinal_Cortex_adj.P.Value", "inhouse_Entorhinal_Cortex_B")
colnames(inhouse_data_clean_top_genes_Cerebellum)<-c("inhouse_Cerebellum_logFC", "inhouse_Cerebellum_CI.L", "inhouse_Cerebellum_CI.R", "inhouse_Cerebellum_AveExpr", "inhouse_Cerebellum_t", "inhouse_Cerebellum_P.Value",  "inhouse_Cerebellum_adj.P.Value", "inhouse_Cerebellum_B")
colnames(inhouse_data_clean_top_genes_Frontal_Cortex)<-c("inhouse_Frontal_Cortex_logFC", "inhouse_Frontal_Cortex_CI.L", "inhouse_Frontal_Cortex_CI.R", "inhouse_Frontal_Cortex_AveExpr", "inhouse_Frontal_Cortex_t", "inhouse_Frontal_Cortex_P.Value",  "inhouse_Frontal_Cortex_adj.P.Value", "inhouse_Frontal_Cortex_B")
colnames(inhouse_data_clean_top_genes_Temporal_Cortex)<-c("inhouse_Temporal_Cortex_logFC", "inhouse_Temporal_Cortex_CI.L", "inhouse_Temporal_Cortex_CI.R", "inhouse_Temporal_Cortex_AveExpr", "inhouse_Temporal_Cortex_t", "inhouse_Temporal_Cortex_P.Value",  "inhouse_Temporal_Cortex_adj.P.Value", "inhouse_Temporal_Cortex_B")
colnames(AMP_MAYOeGWAS_top_genes_Temporal_Cortex)<-c("AMP_MAYOeGWAS_Temporal_Cortex_logFC", "AMP_MAYOeGWAS_Temporal_Cortex_CI.L", "AMP_MAYOeGWAS_Temporal_Cortex_CI.R", "AMP_MAYOeGWAS_Temporal_Cortex_AveExpr", "AMP_MAYOeGWAS_Temporal_Cortex_t", "AMP_MAYOeGWAS_Temporal_Cortex_P.Value", "AMP_MAYOeGWAS_Temporal_Cortex_adj.P.Value", "AMP_MAYOeGWAS_Temporal_Cortex_B")
colnames(AMP_MAYOeGWAS_top_genes_Cerebellum)<-c("AMP_MAYOeGWAS_Cerebellum_logFC", "AMP_MAYOeGWAS_Cerebellum_CI.L", "AMP_MAYOeGWAS_Cerebellum_CI.R", "AMP_MAYOeGWAS_Cerebellum_AveExpr", "AMP_MAYOeGWAS_Cerebellum_t", "AMP_MAYOeGWAS_Cerebellum_P.Value", "AMP_MAYOeGWAS_Cerebellum_adj.P.Value", "AMP_MAYOeGWAS_Cerebellum_B")
colnames(AMP_MSBB_U133A_top_genes_Frontal_Pole)<-c("AMP_MSBB_U133A_Frontal_Pole_logFC", "AMP_MSBB_U133A_Frontal_Pole_CI.L", "AMP_MSBB_U133A_Frontal_Pole_CI.R", "AMP_MSBB_U133A_Frontal_Pole_AveExpr", "AMP_MSBB_U133A_Frontal_Pole_t", "AMP_MSBB_U133A_Frontal_Pole_P.Value", "AMP_MSBB_U133A_Frontal_Pole_adj.P.Value", "AMP_MSBB_U133A_Frontal_Pole_B")
colnames(AMP_MSBB_U133A_top_genes_Precentral_Gyrus)<-c("AMP_MSBB_U133A_Precentral_Gyrus_logFC", "AMP_MSBB_U133A_Precentral_Gyrus_CI.L", "AMP_MSBB_U133A_Precentral_Gyrus_CI.R", "AMP_MSBB_U133A_Precentral_Gyrus_AveExpr", "AMP_MSBB_U133A_Precentral_Gyrus_t", "AMP_MSBB_U133A_Precentral_Gyrus_P.Value", "AMP_MSBB_U133A_Precentral_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Precentral_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus)<-c("AMP_MSBB_U133A_Inferior_Frontal_Gyrus_logFC", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_CI.L", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_CI.R", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_AveExpr", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_t", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_P.Value", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Inferior_Frontal_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex)<-c("AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_logFC", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_CI.L", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_CI.R", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_AveExpr", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_t", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_P.Value", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_adj.P.Value", "AMP_MSBB_U133A_Dorsolateral_Prefrontal_Cortex_B")
colnames(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule)<-c("AMP_MSBB_U133A_Superior_Parietal_Lobule_logFC", "AMP_MSBB_U133A_Superior_Parietal_Lobule_CI.L", "AMP_MSBB_U133A_Superior_Parietal_Lobule_CI.R", "AMP_MSBB_U133A_Superior_Parietal_Lobule_AveExpr", "AMP_MSBB_U133A_Superior_Parietal_Lobule_t", "AMP_MSBB_U133A_Superior_Parietal_Lobule_P.Value", "AMP_MSBB_U133A_Superior_Parietal_Lobule_adj.P.Value", "AMP_MSBB_U133A_Superior_Parietal_Lobule_B")
colnames(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex)<-c("AMP_MSBB_U133A_Prefrontal_Cortex_logFC", "AMP_MSBB_U133A_Prefrontal_Cortex_CI.L", "AMP_MSBB_U133A_Prefrontal_Cortex_CI.R", "AMP_MSBB_U133A_Prefrontal_Cortex_AveExpr", "AMP_MSBB_U133A_Prefrontal_Cortex_t", "AMP_MSBB_U133A_Prefrontal_Cortex_P.Value", "AMP_MSBB_U133A_Prefrontal_Cortex_adj.P.Value", "AMP_MSBB_U133A_Prefrontal_Cortex_B")
colnames(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus)<-c("AMP_MSBB_U133A_Parahippocampal_Gyrus_logFC", "AMP_MSBB_U133A_Parahippocampal_Gyrus_CI.L", "AMP_MSBB_U133A_Parahippocampal_Gyrus_CI.R", "AMP_MSBB_U133A_Parahippocampal_Gyrus_AveExpr", "AMP_MSBB_U133A_Parahippocampal_Gyrus_t", "AMP_MSBB_U133A_Parahippocampal_Gyrus_P.Value", "AMP_MSBB_U133A_Parahippocampal_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Parahippocampal_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Hippocampus)<-c("AMP_MSBB_U133A_Hippocampus_logFC", "AMP_MSBB_U133A_Hippocampus_CI.L", "AMP_MSBB_U133A_Hippocampus_CI.R", "AMP_MSBB_U133A_Hippocampus_AveExpr", "AMP_MSBB_U133A_Hippocampus_t", "AMP_MSBB_U133A_Hippocampus_P.Value", "AMP_MSBB_U133A_Hippocampus_adj.P.Value", "AMP_MSBB_U133A_Hippocampus_B")
colnames(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus)<-c("AMP_MSBB_U133A_Inferior_Temporal_Gyrus_logFC", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_CI.L", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_CI.R", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_t", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_P.Value", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Inferior_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus)<-c("AMP_MSBB_U133A_Middle_Temporal_Gyrus_logFC", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_CI.L", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_CI.R", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_t", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_P.Value", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Middle_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus)<-c("AMP_MSBB_U133A_Superior_Temporal_Gyrus_logFC", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_CI.L", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_CI.R", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_t", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_P.Value", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133A_Superior_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133A_top_genes_Temporal_Pole)<-c("AMP_MSBB_U133A_Temporal_Pole_logFC", "AMP_MSBB_U133A_Temporal_Pole_CI.L", "AMP_MSBB_U133A_Temporal_Pole_CI.R", "AMP_MSBB_U133A_Temporal_Pole_AveExpr", "AMP_MSBB_U133A_Temporal_Pole_t", "AMP_MSBB_U133A_Temporal_Pole_P.Value", "AMP_MSBB_U133A_Temporal_Pole_adj.P.Value", "AMP_MSBB_U133A_Temporal_Pole_B")
colnames(AMP_MSBB_U133B_top_genes_Frontal_Pole)<-c("AMP_MSBB_U133B_Frontal_Pole_logFC", "AMP_MSBB_U133B_Frontal_Pole_CI.L", "AMP_MSBB_U133B_Frontal_Pole_CI.R", "AMP_MSBB_U133B_Frontal_Pole_AveExpr", "AMP_MSBB_U133B_Frontal_Pole_t", "AMP_MSBB_U133B_Frontal_Pole_P.Value", "AMP_MSBB_U133B_Frontal_Pole_adj.P.Value", "AMP_MSBB_U133B_Frontal_Pole_B")
colnames(AMP_MSBB_U133B_top_genes_Precentral_Gyrus)<-c("AMP_MSBB_U133B_Precentral_Gyrus_logFC", "AMP_MSBB_U133B_Precentral_Gyrus_CI.L", "AMP_MSBB_U133B_Precentral_Gyrus_CI.R", "AMP_MSBB_U133B_Precentral_Gyrus_AveExpr", "AMP_MSBB_U133B_Precentral_Gyrus_t", "AMP_MSBB_U133B_Precentral_Gyrus_P.Value", "AMP_MSBB_U133B_Precentral_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Precentral_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus)<-c("AMP_MSBB_U133B_Inferior_Frontal_Gyrus_logFC", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_CI.L", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_CI.R", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_AveExpr", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_t", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_P.Value", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Inferior_Frontal_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex)<-c("AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_logFC", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_CI.L", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_CI.R", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_AveExpr", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_t", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_P.Value", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_adj.P.Value", "AMP_MSBB_U133B_Dorsolateral_Prefrontal_Cortex_B")
colnames(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule)<-c("AMP_MSBB_U133B_Superior_Parietal_Lobule_logFC", "AMP_MSBB_U133B_Superior_Parietal_Lobule_CI.L", "AMP_MSBB_U133B_Superior_Parietal_Lobule_CI.R", "AMP_MSBB_U133B_Superior_Parietal_Lobule_AveExpr", "AMP_MSBB_U133B_Superior_Parietal_Lobule_t", "AMP_MSBB_U133B_Superior_Parietal_Lobule_P.Value", "AMP_MSBB_U133B_Superior_Parietal_Lobule_adj.P.Value", "AMP_MSBB_U133B_Superior_Parietal_Lobule_B")
colnames(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex)<-c("AMP_MSBB_U133B_Prefrontal_Cortex_logFC", "AMP_MSBB_U133B_Prefrontal_Cortex_CI.L", "AMP_MSBB_U133B_Prefrontal_Cortex_CI.R", "AMP_MSBB_U133B_Prefrontal_Cortex_AveExpr", "AMP_MSBB_U133B_Prefrontal_Cortex_t", "AMP_MSBB_U133B_Prefrontal_Cortex_P.Value", "AMP_MSBB_U133B_Prefrontal_Cortex_adj.P.Value", "AMP_MSBB_U133B_Prefrontal_Cortex_B")
colnames(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus)<-c("AMP_MSBB_U133B_Parahippocampal_Gyrus_logFC", "AMP_MSBB_U133B_Parahippocampal_Gyrus_CI.L", "AMP_MSBB_U133B_Parahippocampal_Gyrus_CI.R", "AMP_MSBB_U133B_Parahippocampal_Gyrus_AveExpr", "AMP_MSBB_U133B_Parahippocampal_Gyrus_t", "AMP_MSBB_U133B_Parahippocampal_Gyrus_P.Value", "AMP_MSBB_U133B_Parahippocampal_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Parahippocampal_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Hippocampus)<-c("AMP_MSBB_U133B_Hippocampus_logFC", "AMP_MSBB_U133B_Hippocampus_CI.L", "AMP_MSBB_U133B_Hippocampus_CI.R", "AMP_MSBB_U133B_Hippocampus_AveExpr", "AMP_MSBB_U133B_Hippocampus_t", "AMP_MSBB_U133B_Hippocampus_P.Value", "AMP_MSBB_U133B_Hippocampus_adj.P.Value", "AMP_MSBB_U133B_Hippocampus_B")
colnames(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus)<-c("AMP_MSBB_U133B_Inferior_Temporal_Gyrus_logFC", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_CI.L", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_CI.R", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_t", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_P.Value", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Inferior_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus)<-c("AMP_MSBB_U133B_Middle_Temporal_Gyrus_logFC", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_CI.L", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_CI.R", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_t", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_P.Value", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Middle_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus)<-c("AMP_MSBB_U133B_Superior_Temporal_Gyrus_logFC", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_CI.L", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_CI.R", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_AveExpr", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_t", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_P.Value", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_adj.P.Value", "AMP_MSBB_U133B_Superior_Temporal_Gyrus_B")
colnames(AMP_MSBB_U133B_top_genes_Temporal_Pole)<-c("AMP_MSBB_U133B_Temporal_Pole_logFC", "AMP_MSBB_U133B_Temporal_Pole_CI.L", "AMP_MSBB_U133B_Temporal_Pole_CI.R", "AMP_MSBB_U133B_Temporal_Pole_AveExpr", "AMP_MSBB_U133B_Temporal_Pole_t", "AMP_MSBB_U133B_Temporal_Pole_P.Value", "AMP_MSBB_U133B_Temporal_Pole_adj.P.Value", "AMP_MSBB_U133B_Temporal_Pole_B")

# merge

MyMerge       <- function(x, y){
  df<- merge(x, y, by= "row.names", all=T)
  rownames(df)<- df$Row.names
  df$Row.names<- NULL
  return(df)
}

merged_diff_exprs<-Reduce(MyMerge, list(E_GEOD_28146_top_genes_Hippocampus_CA1,
                                        E_GEOD_29378_top_genes_CA1,
                                        E_GEOD_29378_top_genes_CA3,
                                        E_GEOD_36980_top_genes_Frontal_Cortex,
                                        E_GEOD_36980_top_genes_Hippocampus,
                                        E_GEOD_36980_top_genes_Temporal_Cortex,
                                        E_GEOD_48350_top_genes_Entorhinal_Cortex,
                                        E_GEOD_48350_top_genes_Hippocampus,
                                        E_GEOD_48350_top_genes_Postcentral_Gyrus,
                                        E_GEOD_48350_top_genes_Superior_Frontal_Gyrus,
                                        E_GEOD_1297_top_genes_Hippocampus,
                                        E_GEOD_5281_top_genes_Entorhinal_Cortex,
                                        E_GEOD_5281_top_genes_Hippocampus_CA1,
                                        E_GEOD_5281_top_genes_Medial_Temporal_Gyrus,
                                        E_GEOD_5281_top_genes_Superior_Frontal_Gyrus,
                                        E_GEOD_5281_top_genes_Posterior_Singulate,
                                        E_GEOD_5281_top_genes_Primary_Visual_Cortex,
                                        inhouse_data_clean_top_genes_Entorhinal_Cortex,
                                        inhouse_data_clean_top_genes_Cerebellum,
                                        inhouse_data_clean_top_genes_Frontal_Cortex,
                                        inhouse_data_clean_top_genes_Temporal_Cortex,
                                        AMP_MAYOeGWAS_top_genes_Temporal_Cortex,
                                        AMP_MAYOeGWAS_top_genes_Cerebellum,
                                        AMP_MSBB_U133A_top_genes_Frontal_Pole,
                                        AMP_MSBB_U133A_top_genes_Precentral_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex,
                                        AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule,
                                        AMP_MSBB_U133A_top_genes_Prefrontal_Cortex,
                                        AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Hippocampus,
                                        AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus,
                                        AMP_MSBB_U133A_top_genes_Temporal_Pole,
                                        AMP_MSBB_U133B_top_genes_Frontal_Pole,
                                        AMP_MSBB_U133B_top_genes_Precentral_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex,
                                        AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule,
                                        AMP_MSBB_U133B_top_genes_Prefrontal_Cortex,
                                        AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Hippocampus,
                                        AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus,
                                        AMP_MSBB_U133B_top_genes_Temporal_Pole))


dim(merged_diff_exprs)
head(merged_diff_exprs)[100:110]


# separate logFC and adj.P

logFC_values<-grep("logFC", colnames(merged_diff_exprs))
logFC_data<-as.matrix(subset(merged_diff_exprs[c(logFC_values)]))
head(logFC_data)
dim(logFC_data)

adj_P_values<-grep("adj.P.Val", colnames(merged_diff_exprs))
adj_P_data<-as.matrix(subset(merged_diff_exprs[c(adj_P_values)]))
head(adj_P_data)
dim(adj_P_data)

# change colnames

colnames(adj_P_data)<-gsub("_adj.P.Val", "", colnames(adj_P_data))
head(adj_P_data)

colnames(logFC_data)<-gsub("_logFC", "", colnames(logFC_data))
head(logFC_data)

##### HEATMAP - adj.P - all probes #####

#"error in hclustfun(distr) : NA/NaN/Inf in foreign function call (arg 11)" - remove rows where Na in more than 60% (minimum number needed for heatmap through trial and error)

# remove row (probe) if more than 40% missing
adj_P_data_min_na<-adj_P_data[-which(rowMeans(is.na(adj_P_data))>0.4),]
#dim(adj_P_data_min_na)

dim(adj_P_data)
dim(adj_P_data_min_na)

# heatmap

mycolour<-colorRampPalette(c("red", "yellow", "black", "purple", "blue"))(n = 100)

heatmap.2(log2(t(adj_P_data_min_na)),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Heatmap of adj.P.Value",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### HEATMAP - adj.P - common probes #####

adj_P_common_genes<-na.omit(adj_P_data)
dim(adj_P_data)
dim(adj_P_common_genes)

# heatmap

heatmap.2(log2(t(adj_P_common_genes)),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "row",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.6,
          las=2,
          main="Heatmap of adj.P.Value using common probes only",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### HEATMAP - logFC - all probes #####

#"error in hclustfun(distr) : NA/NaN/Inf in foreign function call (arg 11)" - remove rows where Na in more than 60% (minimum number needed for heatmap through trial and error)

# remove row (probe) if more than 60% missing
logFC_data_min_na<-logFC_data[-which(rowMeans(is.na(logFC_data))>0.4),]
dim(logFC_data_min_na)

# heatmap

heatmap.2(t(logFC_data_min_na),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Heatmap of LogFC",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### HEATMAP - logFC - common probes #####

# common probes only
logFC_data_common_probes<-na.omit(logFC_data)
dim(logFC_data_common_probes)

# heatmap

mycolour<-colorRampPalette(c("red", "orange", "yellow", "black", "purple", "blue", "green"))(n = 200)

heatmap.2(t(logFC_data_common_probes),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Heatmap of LogFC using common probes only",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### HEATMAP - LogFC - only 1/0 #####

#convert pos logFC to 1 and neg logFC to -1.
logFC_data_common_probes
dim(logFC_data_common_probes)
logFC_data_common_probes_binary<-apply(logFC_data_common_probes, 2, as.numeric)
head(logFC_data_common_probes_binary)[,1:5]
logFC_data_common_probes_binary[logFC_data_common_probes_binary<0]<-"-1"
head(logFC_data_common_probes_binary)[,1:5]
logFC_data_common_probes_binary<-apply(logFC_data_common_probes_binary, 2, as.numeric)
logFC_data_common_probes_binary[logFC_data_common_probes_binary>0]<-"1"
logFC_data_common_probes_binary<-apply(logFC_data_common_probes_binary, 2, as.numeric)
head(logFC_data_common_probes_binary)[,1:5]

dim(logFC_data_common_probes_binary)
str(logFC_data_common_probes_binary)
class(logFC_data_common_probes_binary)

heatmap.2(t(logFC_data_common_probes_binary),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Heatmap of LogFC using common probes only",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### PLOT TO PDF #####

# create folder within work_dir

setwd(work_dir)

dir.create(paste(work_dir,"Heatmaps", sep="/"))

pdf_plot_dir=paste(work_dir,"Heatmaps", sep="/")

setwd(pdf_plot_dir)

#plot DE adj P value all probes

pdf("Differential_Expression_Heatmap_adj_P_value_all_probes.pdf")
heatmap.2(log2(t(adj_P_data_min_na)),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Log2 adj.P.Value all probes",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")
dev.off()

#plot DE adj P value common probes

pdf("Differential_Expression_Heatmap_adj_P_value_common_probes.pdf")
heatmap.2(log2(t(adj_P_common_genes)),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "row",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.6,
          las=2,
          main="Log2 adj.P.Value common probes",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")
dev.off()

pdf("Differential_Expression_Heatmap_logFC_all_probes.pdf")
heatmap.2(t(logFC_data_min_na),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="LogFC all probes",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")
dev.off()

pdf("Differential_Expression_Heatmap_logFC_common_probes.pdf")
heatmap.2(t(logFC_data_common_probes),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="LogFC common probes",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")
dev.off()

pdf("Differential_Expression_Heatmap_logFC_common_probes_0_and_1.pdf")
heatmap.2(t(logFC_data_common_probes_binary),
          Rowv=TRUE,
          Colv = TRUE,
          scale="none",
          col=mycolour,
          #col=redgreen(100),
          #dendrogram = "column",
          dendrogram = "both",
          notecex=1,
          notecol="black",
          na.color=par("bg"),
          margins = c(5, 15),
          sepwidth=c(0.01,0.01),
          key=TRUE,
          symkey=TRUE,
          trace="none",
          cexRow=0.6,
          cexCol=0.8,
          las=2,
          main="Heatmap of LogFC using common probes only",
          xlab="Entrez Gene ID", 
          ylab="Dataset Tissue")

##### WITHOUT AMP MSBB #####




##### SAVE IMAGE #####

setwd(work_dir)

save.image("DE ANALYSIS.Rdata")
