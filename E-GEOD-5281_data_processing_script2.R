
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                     AD - E-GEOD-5281  - PIEPLINE V02                                                   #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Affymetrix
# EXPRESSION CHIP - U133 plus 2.0
# NUMBER OF SAMPLES - 161
# BRAIN REGIONS - Entorhinal cortex, hippocampus, posterior cingulate cortex, superior frontal gyrus, middle temporal gyrus, vusual cortex
# 
#
# REPROCESSING DATA USING NEW PIPELINE

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

raw_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/E-GEOD-5281/Raw_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/2.Re-process_with_new_pipeline/AD/E-GEOD-5281/Pre-Processing"

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
library(hgu133plus2.db)

##### DOWNLOAD RAW DATA #####

setwd(raw_dir)

brain_data_raw=getAE("E-GEOD-5281 ", type = "full")

##### CREATE R EXPRESSION OBJECT FROM RAW DATA #####

#convert MAGE-TAB files into expresssion set - using processed data

test<-ReadAffy()

brain_data<-test

##### PLOTS OF RAW DATA ####

setwd(work_dir)

boxplot(log2(exprs(brain_data)))
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

#normalise

brain_data_normalised<-rsn(log2(exprs(brain_data_background_corrected)))

brain_data_normalised_as_data_frame<-as.data.frame(brain_data_normalised)

# change colnames - remove ".CEL" from each column name

sample_id_cel<-colnames(brain_data_normalised_as_data_frame)

sample_id_corrected<-gsub(".CEL", "", sample_id_cel)
head(sample_id_corrected)
head(sample_id_cel)

colnames(brain_data_normalised_as_data_frame)<-sample_id_corrected


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

#manually read in pheno data

setwd(raw_dir)

phenotype_data<-read.table("E-GEOD-5281.sdrf.txt", sep="\t", header=T, comment.char="", stringsAsFactor=F, fill=T, quote="")

head(phenotype_data)

#move sample id to rownames

phenotype_data$sample_id<-gsub(" 1", "", phenotype_data$Source.Name)

rownames(phenotype_data)<-phenotype_data$sample_id

head(phenotype_data)

table(phenotype_data$Comment..Sample_source_name.)

table(phenotype_data$Characteristics.disease.state.)

table(phenotype_data$Characteristics.organ.region.)

table(phenotype_data$FactorValue..SEX..1)

table(phenotype_data$FactorValue..ETHNICITY.)

##### SAVE PROBE IDs #####

probe_ids<-rownames(brain_data_normalised)

head(probe_ids)
length(probe_ids)

setwd(work_dir)

write(file = "E-GEOD-5281_probe_IDs.txt", probe_ids, sep="\t")

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(phenotype_data))

##### GENDER CHECK ##### 

# XIST:
# 1427262_at 
# 1427263_at 
# 1436936_s_at
# 1442137_at 
# 1458435_at 

# PRKY

# individual gender uknown - total 16 males and 19 females. estimatung gender based on XIST gene

##### PROBE ID DETECTION #####

# separate case control - Characteristics.disease.state. column 

sample_case_control_ID<-phenotype_data$Source.Name

sample_case_control_ID

case_ID<-rownames(phenotype_data[phenotype_data$Characteristics.disease.state.=="Alzheimer's Disease",])

control_ID<-rownames(phenotype_data[phenotype_data$Characteristics.disease.state.=="normal",])

case_ID

control_ID

case_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%case_ID]
control_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%control_ID]

head(case_exprs)
dim(case_exprs)
head(control_exprs)
dim(control_exprs)

# separate by brain region - Characteristics.organ.region.

Entorhinal_Cortex<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="Entorhinal Cortex",])
hippocampus<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="hippocampus",])
Medial_Temporal_Gyrus<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="Medial Temporal Gyrus",])
Posterior_Singulate<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="Posterior Singulate",])
Primary_Visual_Cortex<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="Primary Visual Cortex",])
Superior_Frontal_Gyrus<-rownames(phenotype_data[phenotype_data$Characteristics.organ.region.=="Superior Frontal Gyrus",])

Entorhinal_Cortex
hippocampus
Medial_Temporal_Gyrus
Posterior_Singulate
Primary_Visual_Cortex
Superior_Frontal_Gyrus

# separate by brain region and disorder

Entorhinal_Cortex_case_exprs<-case_exprs[,colnames(case_exprs)%in%Entorhinal_Cortex]
hippocampus_case_exprs<-case_exprs[,colnames(case_exprs)%in%hippocampus]
Medial_Temporal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Medial_Temporal_Gyrus]
Posterior_Singulate_case_exprs<-case_exprs[,colnames(case_exprs)%in%Posterior_Singulate]
Primary_Visual_Cortex_case_exprs<-case_exprs[,colnames(case_exprs)%in%Primary_Visual_Cortex]
Superior_Frontal_Gyrus_case_exprs<-case_exprs[,colnames(case_exprs)%in%Superior_Frontal_Gyrus]

Entorhinal_Cortex_control_exprs<-control_exprs[,colnames(control_exprs)%in%Entorhinal_Cortex]
hippocampus_control_exprs<-control_exprs[,colnames(control_exprs)%in%hippocampus]
Medial_Temporal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Medial_Temporal_Gyrus]
Posterior_Singulate_control_exprs<-control_exprs[,colnames(control_exprs)%in%Posterior_Singulate]
Primary_Visual_Cortex_control_exprs<-control_exprs[,colnames(control_exprs)%in%Primary_Visual_Cortex]
Superior_Frontal_Gyrus_control_exprs<-control_exprs[,colnames(control_exprs)%in%Superior_Frontal_Gyrus]

# check dataframe

dim(Entorhinal_Cortex_case_exprs)
dim(hippocampus_case_exprs)
dim(Medial_Temporal_Gyrus_case_exprs)
dim(Posterior_Singulate_case_exprs)
dim(Primary_Visual_Cortex_case_exprs)
dim(Superior_Frontal_Gyrus_case_exprs)

head(Entorhinal_Cortex_case_exprs)[1:5]
head(hippocampus_case_exprs)[1:5]
head(Medial_Temporal_Gyrus_case_exprs)[1:5]
head(Posterior_Singulate_case_exprs)[1:5]
head(Primary_Visual_Cortex_case_exprs)[1:5]
head(Superior_Frontal_Gyrus_case_exprs)[1:5]

dim(Entorhinal_Cortex_control_exprs)
dim(hippocampus_control_exprs)
dim(Medial_Temporal_Gyrus_control_exprs)
dim(Posterior_Singulate_control_exprs)
dim(Primary_Visual_Cortex_control_exprs)
dim(Superior_Frontal_Gyrus_control_exprs)

head(Entorhinal_Cortex_control_exprs)[1:5]
head(hippocampus_control_exprs)[1:5]
head(Medial_Temporal_Gyrus_control_exprs)[1:5]
head(Posterior_Singulate_control_exprs)[1:5]
head(Primary_Visual_Cortex_control_exprs)[1:5]
head(Superior_Frontal_Gyrus_control_exprs)[1:5]

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
  # return good probes
  return(good_probes)
}

#

# apply function - case

Entorhinal_Cortex_case_expressed_probes_list<-extract_good_probe_list(Entorhinal_Cortex_case_exprs, 0.9, 0.8)
head(Entorhinal_Cortex_case_expressed_probes_list)
length(Entorhinal_Cortex_case_expressed_probes_list)
nrow(Entorhinal_Cortex_case_exprs)

hippocampus_case_expressed_probes_list<-extract_good_probe_list(hippocampus_case_exprs, 0.9, 0.8)
head(hippocampus_case_expressed_probes_list)
length(hippocampus_case_expressed_probes_list)
nrow(hippocampus_case_exprs)

Medial_Temporal_Gyrus_case_expressed_probes_list<-extract_good_probe_list(Medial_Temporal_Gyrus_case_exprs, 0.9, 0.8)
head(Medial_Temporal_Gyrus_case_expressed_probes_list)
length(Medial_Temporal_Gyrus_case_expressed_probes_list)
nrow(Medial_Temporal_Gyrus_case_exprs)

Posterior_Singulate_case_expressed_probes_list<-extract_good_probe_list(Posterior_Singulate_case_exprs, 0.9, 0.8)
head(Posterior_Singulate_case_expressed_probes_list)
length(Posterior_Singulate_case_expressed_probes_list)
nrow(Posterior_Singulate_case_exprs)

Primary_Visual_Cortex_case_expressed_probes_list<-extract_good_probe_list(Primary_Visual_Cortex_case_exprs, 0.9, 0.8)
head(Primary_Visual_Cortex_case_expressed_probes_list)
length(Primary_Visual_Cortex_case_expressed_probes_list)
nrow(Primary_Visual_Cortex_case_exprs)

Superior_Frontal_Gyrus_case_expressed_probes_list<-extract_good_probe_list(Superior_Frontal_Gyrus_case_exprs, 0.9, 0.8)
head(Superior_Frontal_Gyrus_case_expressed_probes_list)
length(Superior_Frontal_Gyrus_case_expressed_probes_list)
nrow(Superior_Frontal_Gyrus_case_exprs)

#apply function to control

Entorhinal_Cortex_control_expressed_probes_list<-extract_good_probe_list(Entorhinal_Cortex_control_exprs, 0.9, 0.8)
head(Entorhinal_Cortex_control_expressed_probes_list)
length(Entorhinal_Cortex_control_expressed_probes_list)
nrow(Entorhinal_Cortex_control_exprs)

hippocampus_control_expressed_probes_list<-extract_good_probe_list(hippocampus_control_exprs, 0.9, 0.8)
head(hippocampus_control_expressed_probes_list)
length(hippocampus_control_expressed_probes_list)
nrow(hippocampus_control_exprs)

Medial_Temporal_Gyrus_control_expressed_probes_list<-extract_good_probe_list(Medial_Temporal_Gyrus_control_exprs, 0.9, 0.8)
head(Medial_Temporal_Gyrus_control_expressed_probes_list)
length(Medial_Temporal_Gyrus_control_expressed_probes_list)
nrow(Medial_Temporal_Gyrus_control_exprs)

Posterior_Singulate_control_expressed_probes_list<-extract_good_probe_list(Posterior_Singulate_control_exprs, 0.9, 0.8)
head(Posterior_Singulate_control_expressed_probes_list)
length(Posterior_Singulate_control_expressed_probes_list)
nrow(Posterior_Singulate_control_exprs)

Primary_Visual_Cortex_control_expressed_probes_list<-extract_good_probe_list(Primary_Visual_Cortex_control_exprs, 0.9, 0.8)
head(Primary_Visual_Cortex_control_expressed_probes_list)
length(Primary_Visual_Cortex_control_expressed_probes_list)
nrow(Primary_Visual_Cortex_control_exprs)

Superior_Frontal_Gyrus_control_expressed_probes_list<-extract_good_probe_list(Superior_Frontal_Gyrus_control_exprs, 0.9, 0.8)
head(Superior_Frontal_Gyrus_control_expressed_probes_list)
length(Superior_Frontal_Gyrus_control_expressed_probes_list)
nrow(Superior_Frontal_Gyrus_control_exprs)

# merge list of good probes from both case + control + and all brain regions, sort and keep unique values

good_probe_list<-unique(sort(c(Entorhinal_Cortex_case_expressed_probes_list,
                               hippocampus_case_expressed_probes_list,
                               Medial_Temporal_Gyrus_case_expressed_probes_list,
                               Posterior_Singulate_case_expressed_probes_list,
                               Primary_Visual_Cortex_case_expressed_probes_list,
                               Superior_Frontal_Gyrus_case_expressed_probes_list,
                               Entorhinal_Cortex_control_expressed_probes_list,
                               hippocampus_control_expressed_probes_list,
                               Medial_Temporal_Gyrus_control_expressed_probes_list,
                               Posterior_Singulate_control_expressed_probes_list,
                               Primary_Visual_Cortex_control_expressed_probes_list,
                               Superior_Frontal_Gyrus_control_expressed_probes_list)))

length(good_probe_list)


# extract good probes from dataset

brain_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]

Entorhinal_Cortex_case_exprs_good_probes<-Entorhinal_Cortex_case_exprs[rownames(Entorhinal_Cortex_case_exprs)%in%good_probe_list,]

hippocampus_case_exprs_good_probes<-hippocampus_case_exprs[rownames(hippocampus_case_exprs)%in%good_probe_list,]

Medial_Temporal_Gyrus_case_exprs_good_probes<-Medial_Temporal_Gyrus_case_exprs[rownames(Medial_Temporal_Gyrus_case_exprs)%in%good_probe_list,]

Posterior_Singulate_case_exprs_good_probes<-Posterior_Singulate_case_exprs[rownames(Posterior_Singulate_case_exprs)%in%good_probe_list,]

Primary_Visual_Cortex_case_exprs_good_probes<-Primary_Visual_Cortex_case_exprs[rownames(Primary_Visual_Cortex_case_exprs)%in%good_probe_list,]

Superior_Frontal_Gyrus_case_exprs_good_probes<-Superior_Frontal_Gyrus_case_exprs[rownames(Superior_Frontal_Gyrus_case_exprs)%in%good_probe_list,]

Entorhinal_Cortex_control_exprs_good_probes<-Entorhinal_Cortex_control_exprs[rownames(Entorhinal_Cortex_control_exprs)%in%good_probe_list,]

hippocampus_control_exprs_good_probes<-hippocampus_control_exprs[rownames(hippocampus_control_exprs)%in%good_probe_list,]

Medial_Temporal_Gyrus_control_exprs_good_probes<-Medial_Temporal_Gyrus_control_exprs[rownames(Medial_Temporal_Gyrus_control_exprs)%in%good_probe_list,]

Posterior_Singulate_control_exprs_good_probes<-Posterior_Singulate_control_exprs[rownames(Posterior_Singulate_control_exprs)%in%good_probe_list,]

Primary_Visual_Cortex_control_exprs_good_probes<-Primary_Visual_Cortex_control_exprs[rownames(Primary_Visual_Cortex_control_exprs)%in%good_probe_list,]

Superior_Frontal_Gyrus_control_exprs_good_probes<-Superior_Frontal_Gyrus_control_exprs[rownames(Superior_Frontal_Gyrus_control_exprs)%in%good_probe_list,]


# check dataframe - all should be 6658

head(Entorhinal_Cortex_case_exprs_good_probes)[1:5]
dim(Entorhinal_Cortex_case_exprs_good_probes)

head(hippocampus_case_exprs_good_probes)[1:5]
dim(hippocampus_case_exprs_good_probes)

head(Medial_Temporal_Gyrus_case_exprs_good_probes)[1:5]
dim(Medial_Temporal_Gyrus_case_exprs_good_probes)

head(Posterior_Singulate_case_exprs_good_probes)[1:5]
dim(Posterior_Singulate_case_exprs_good_probes)

head(Primary_Visual_Cortex_case_exprs_good_probes)[1:5]
dim(Primary_Visual_Cortex_case_exprs_good_probes)

head(Superior_Frontal_Gyrus_case_exprs_good_probes)[1:5]
dim(Superior_Frontal_Gyrus_case_exprs_good_probes)

head(Entorhinal_Cortex_control_exprs_good_probes)[1:5]
dim(Entorhinal_Cortex_control_exprs_good_probes)

head(hippocampus_control_exprs_good_probes)[1:5]
dim(hippocampus_control_exprs_good_probes)

head(Medial_Temporal_Gyrus_control_exprs_good_probes)[1:5]
dim(Medial_Temporal_Gyrus_control_exprs_good_probes)

head(Posterior_Singulate_control_exprs_good_probes)[1:5]
dim(Posterior_Singulate_control_exprs_good_probes)

head(Primary_Visual_Cortex_control_exprs_good_probes)[1:5]
dim(Primary_Visual_Cortex_control_exprs_good_probes)

head(Superior_Frontal_Gyrus_control_exprs_good_probes)[1:5]
dim(Superior_Frontal_Gyrus_control_exprs_good_probes)

dim(brain_data_normalised)

# create dataframe with good probes - contain all brain regions, case and control

brain_data_normalised_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]
dim(brain_data_normalised_exprs_good_probes)
head(brain_data_normalised_exprs_good_probes)[1:5]

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

phenotype<-phenotype_data[10]
head(phenotype)

brain_region<-phenotype_data[13]
head(brain_region)

# phenotyp and sample iD
Diagnosis_lookup<-phenotype_data[10]
colnames(Diagnosis_lookup)<-"Diagnosis"
Diagnosis_lookup

# order of samples in expression data
Diagnosis_temp<-colnames(brain_data_normalised_exprs_good_probes)
Diagnosis_temp

# match order
Diagnosis_pca<-Diagnosis_lookup[match(Diagnosis_temp, rownames(Diagnosis_lookup)),]

# assign color to group
Diagnosis_pca_color<-labels2colors(as.character(Diagnosis_pca))

# pca plot - color by disease - case/control
plot(pca$rotation[,1:2], main=" PCA plot coloured by Disease Status",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomleft', unique(Diagnosis_pca), fill=unique(Diagnosis_pca_color))

#pca plot - color by brain region

# phenotype and sample iD
brain_region_lookup<-phenotype_data[13]
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
legend('bottomleft', unique(brain_region_pca), fill=unique(brain_region_pca_color))

##### SVA ####

# separate by brain region

Entorhinal_Cortex_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Entorhinal_Cortex]
Hippocampus_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%hippocampus]
Medial_Temporal_Gyrus_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Medial_Temporal_Gyrus]
Posterior_Singulate_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Posterior_Singulate]
Primary_Visual_Cortex_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Primary_Visual_Cortex]
Superior_Frontal_Gyrus_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Superior_Frontal_Gyrus]

dim(Entorhinal_Cortex_good_probes)
dim(Hippocampus_good_probes)
dim(Medial_Temporal_Gyrus_good_probes)
dim(Posterior_Singulate_good_probes)
dim(Primary_Visual_Cortex_good_probes)
dim(Superior_Frontal_Gyrus_good_probes)


## add in diagnosis column

Entorhinal_Cortex_good_probes<-merge(Diagnosis, t(Entorhinal_Cortex_good_probes), by="row.names")
rownames(Entorhinal_Cortex_good_probes)<-Entorhinal_Cortex_good_probes$Row.names
Entorhinal_Cortex_good_probes$Row.names<-NULL
dim(Entorhinal_Cortex_good_probes)
head(Entorhinal_Cortex_good_probes)[1:5]

Hippocampus_good_probes<-merge(Diagnosis, t(Hippocampus_good_probes), by="row.names")
rownames(Hippocampus_good_probes)<-Hippocampus_good_probes$Row.names
Hippocampus_good_probes$Row.names<-NULL
dim(Hippocampus_good_probes)
head(Hippocampus_good_probes)[1:5]

Medial_Temporal_Gyrus_good_probes<-merge(Diagnosis, t(Medial_Temporal_Gyrus_good_probes), by="row.names")
rownames(Medial_Temporal_Gyrus_good_probes)<-Medial_Temporal_Gyrus_good_probes$Row.names
Medial_Temporal_Gyrus_good_probes$Row.names<-NULL
dim(Medial_Temporal_Gyrus_good_probes)
head(Medial_Temporal_Gyrus_good_probes)[1:5]

Posterior_Singulate_good_probes<-merge(Diagnosis, t(Posterior_Singulate_good_probes), by="row.names")
rownames(Posterior_Singulate_good_probes)<-Posterior_Singulate_good_probes$Row.names
Posterior_Singulate_good_probes$Row.names<-NULL
dim(Posterior_Singulate_good_probes)
head(Posterior_Singulate_good_probes)[1:5]

Primary_Visual_Cortex_good_probes<-merge(Diagnosis, t(Primary_Visual_Cortex_good_probes), by="row.names")
rownames(Primary_Visual_Cortex_good_probes)<-Primary_Visual_Cortex_good_probes$Row.names
Primary_Visual_Cortex_good_probes$Row.names<-NULL
dim(Primary_Visual_Cortex_good_probes)
head(Primary_Visual_Cortex_good_probes)[1:5]

Superior_Frontal_Gyrus_good_probes<-merge(Diagnosis, t(Superior_Frontal_Gyrus_good_probes), by="row.names")
rownames(Superior_Frontal_Gyrus_good_probes)<-Superior_Frontal_Gyrus_good_probes$Row.names
Superior_Frontal_Gyrus_good_probes$Row.names<-NULL
dim(Superior_Frontal_Gyrus_good_probes)
head(Superior_Frontal_Gyrus_good_probes)[1:5]

# create sva function

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_diagnosis[1]
  dataset_exprs<-t(sorted_by_diagnosis[2:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis, data=dataset_pheno)
  # check number of SV in data
  num.sv(dataset_exprs, mod, method="leek")
}

# check sv

check_SV_in_data(Entorhinal_Cortex_good_probes)
check_SV_in_data(Hippocampus_good_probes)
check_SV_in_data(Medial_Temporal_Gyrus_good_probes)
check_SV_in_data(Posterior_Singulate_good_probes)
check_SV_in_data(Primary_Visual_Cortex_good_probes)
check_SV_in_data(Superior_Frontal_Gyrus_good_probes)

# create function to adjust for SV

adjust_for_sva<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_sva_pheno<-sorted_by_diagnosis[1]
  dataset_sva_exprs<-t(sorted_by_diagnosis[2:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis, data=dataset_sva_pheno)
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

Entorhinal_Cortex_good_probes_sva_adjusted<-adjust_for_sva(Entorhinal_Cortex_good_probes)
Hippocampus_good_probes_sva_adjusted<-adjust_for_sva(Hippocampus_good_probes)
Medial_Temporal_Gyrus_good_probes_sva_adjusted<-adjust_for_sva(Medial_Temporal_Gyrus_good_probes)
Posterior_Singulate_good_probes_sva_adjusted<-adjust_for_sva(Posterior_Singulate_good_probes)
Primary_Visual_Cortex_good_probes_sva_adjusted<-adjust_for_sva(Primary_Visual_Cortex_good_probes)
Superior_Frontal_Gyrus_good_probes_sva_adjusted<-adjust_for_sva(Superior_Frontal_Gyrus_good_probes)

# separate case and control

Entorhinal_Cortex_good_probes_sva_adjusted_case<-Entorhinal_Cortex_good_probes_sva_adjusted[Entorhinal_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]
Hippocampus_good_probes_sva_adjusted_case<-Hippocampus_good_probes_sva_adjusted[Hippocampus_good_probes_sva_adjusted$Diagnosis=="AD",]
Medial_Temporal_Gyrus_good_probes_sva_adjusted_case<-Medial_Temporal_Gyrus_good_probes_sva_adjusted[Medial_Temporal_Gyrus_good_probes_sva_adjusted$Diagnosis=="AD",]
Posterior_Singulate_good_probes_sva_adjusted_case<-Posterior_Singulate_good_probes_sva_adjusted[Posterior_Singulate_good_probes_sva_adjusted$Diagnosis=="AD",]
Primary_Visual_Cortex_good_probes_sva_adjusted_case<-Primary_Visual_Cortex_good_probes_sva_adjusted[Primary_Visual_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]
Superior_Frontal_Gyrus_good_probes_sva_adjusted_case<-Superior_Frontal_Gyrus_good_probes_sva_adjusted[Superior_Frontal_Gyrus_good_probes_sva_adjusted$Diagnosis=="AD",]

Entorhinal_Cortex_good_probes_sva_adjusted_control<-Entorhinal_Cortex_good_probes_sva_adjusted[Entorhinal_Cortex_good_probes_sva_adjusted$Diagnosis=="Control",]
Hippocampus_good_probes_sva_adjusted_control<-Hippocampus_good_probes_sva_adjusted[Hippocampus_good_probes_sva_adjusted$Diagnosis=="Control",]
Medial_Temporal_Gyrus_good_probes_sva_adjusted_control<-Medial_Temporal_Gyrus_good_probes_sva_adjusted[Medial_Temporal_Gyrus_good_probes_sva_adjusted$Diagnosis=="Control",]
Posterior_Singulate_good_probes_sva_adjusted_control<-Posterior_Singulate_good_probes_sva_adjusted[Posterior_Singulate_good_probes_sva_adjusted$Diagnosis=="Control",]
Primary_Visual_Cortex_good_probes_sva_adjusted_control<-Primary_Visual_Cortex_good_probes_sva_adjusted[Primary_Visual_Cortex_good_probes_sva_adjusted$Diagnosis=="Control",]
Superior_Frontal_Gyrus_good_probes_sva_adjusted_control<-Superior_Frontal_Gyrus_good_probes_sva_adjusted[Superior_Frontal_Gyrus_good_probes_sva_adjusted$Diagnosis=="Control",]

dim(Entorhinal_Cortex_good_probes_sva_adjusted_case)
dim(Hippocampus_good_probes_sva_adjusted_case)
dim(Medial_Temporal_Gyrus_good_probes_sva_adjusted_case)
dim(Posterior_Singulate_good_probes_sva_adjusted_case)
dim(Primary_Visual_Cortex_good_probes_sva_adjusted_case)
dim(Superior_Frontal_Gyrus_good_probes_sva_adjusted_case)


dim(Entorhinal_Cortex_good_probes_sva_adjusted_control)
dim(Hippocampus_good_probes_sva_adjusted_control)
dim(Medial_Temporal_Gyrus_good_probes_sva_adjusted_control)
dim(Posterior_Singulate_good_probes_sva_adjusted_control)
dim(Primary_Visual_Cortex_good_probes_sva_adjusted_control)
dim(Superior_Frontal_Gyrus_good_probes_sva_adjusted_control)

##### SAMPLE NETWORK PLOT ######

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

##### SAMPLE NETWORK PLOT ######

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

Entorhinal_Cortex_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_case, -3)
Entorhinal_Cortex_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control, -3)

# run sample network on Hippocampus

Hippocampus_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Hippocampus_good_probes_sva_adjusted_case, -3)
Hippocampus_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Hippocampus_good_probes_sva_adjusted_control, -3)

# run sample network on Medial_Temporal_Gyrus

Medial_Temporal_Gyrus_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Medial_Temporal_Gyrus_good_probes_sva_adjusted_case, -3)
Medial_Temporal_Gyrus_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Medial_Temporal_Gyrus_good_probes_sva_adjusted_control, -3)

# run sample network on Posterior_Singulate

Posterior_Singulate_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Posterior_Singulate_good_probes_sva_adjusted_case, -3)
Posterior_Singulate_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Posterior_Singulate_good_probes_sva_adjusted_control, -3)

# run sample network on Primary_Visual_Cortex

Primary_Visual_Cortex_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Primary_Visual_Cortex_good_probes_sva_adjusted_case, -3)
Primary_Visual_Cortex_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Primary_Visual_Cortex_good_probes_sva_adjusted_control, -3)

# run sample network on Superior_Frontal_Gyrus

Superior_Frontal_Gyrus_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Superior_Frontal_Gyrus_good_probes_sva_adjusted_case, -3)
Superior_Frontal_Gyrus_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Superior_Frontal_Gyrus_good_probes_sva_adjusted_control, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("Entorhinal_Cortex_case_sample_network_analysis.pdf")
run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Entorhinal_Cortex_control_sample_network_analysis.pdf")
run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Hippocampus_case_sample_network_analysis.pdf")
run_sample_network_plot(Hippocampus_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Hippocampus_control_sample_network_analysis.pdf")
run_sample_network_plot(Hippocampus_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Medial_Temporal_Gyrus_case_sample_network_analysis.pdf")
run_sample_network_plot(Medial_Temporal_Gyrus_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Medial_Temporal_Gyrus_control_sample_network_analysis.pdf")
run_sample_network_plot(Medial_Temporal_Gyrus_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Posterior_Singulate_case_sample_network_analysis.pdf")
run_sample_network_plot(Posterior_Singulate_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Posterior_Singulate_control_sample_network_analysis.pdf")
run_sample_network_plot(Posterior_Singulate_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Primary_Visual_Cortex_case_sample_network_analysis.pdf")
run_sample_network_plot(Primary_Visual_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Primary_Visual_Cortex_control_sample_network_analysis.pdf")
run_sample_network_plot(Primary_Visual_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Superior_Frontal_Gyrus_case_sample_network_analysis.pdf")
run_sample_network_plot(Superior_Frontal_Gyrus_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Superior_Frontal_Gyrus_control_sample_network_analysis.pdf")
run_sample_network_plot(Superior_Frontal_Gyrus_good_probes_sva_adjusted_control, -3)
dev.off()

setwd(work_dir)

##### CREATE QC'd DATASET #####
# extract sample ID's from QC'd sample network file

# check colnames same in all dataframes
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Hippocampus_good_probes_sva_adjusted_control_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Superior_Frontal_Gyrus_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Superior_Frontal_Gyrus_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Medial_Temporal_Gyrus_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Medial_Temporal_Gyrus_good_probes_sva_adjusted_control_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Entorhinal_Cortex_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Entorhinal_Cortex_good_probes_sva_adjusted_control_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Posterior_Singulate_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Posterior_Singulate_good_probes_sva_adjusted_control_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Primary_Visual_Cortex_good_probes_sva_adjusted_case_QC))==F
any(colnames(Hippocampus_good_probes_sva_adjusted_case_QC)==colnames(Primary_Visual_Cortex_good_probes_sva_adjusted_control_QC))==F


brain_data_QCd<-rbind(Hippocampus_good_probes_sva_adjusted_case_QC,
                      Hippocampus_good_probes_sva_adjusted_control_QC,
                      Superior_Frontal_Gyrus_good_probes_sva_adjusted_case_QC,
                      Superior_Frontal_Gyrus_good_probes_sva_adjusted_control_QC,
                      Medial_Temporal_Gyrus_good_probes_sva_adjusted_case_QC,
                      Medial_Temporal_Gyrus_good_probes_sva_adjusted_control_QC,
                      Posterior_Singulate_good_probes_sva_adjusted_case_QC,
                      Posterior_Singulate_good_probes_sva_adjusted_control_QC,
                      Primary_Visual_Cortex_good_probes_sva_adjusted_case_QC,
                      Primary_Visual_Cortex_good_probes_sva_adjusted_control_QC,
                      Entorhinal_Cortex_good_probes_sva_adjusted_case_QC,
                      Entorhinal_Cortex_good_probes_sva_adjusted_control_QC)

dim(brain_data_QCd)
head(brain_data_QCd)[1:5]

##### PCA ON CLEAN DATA ####

# sumary of pca
summary_pca_QCd<-summary(t(brain_data_QCd[2:dim(brain_data_QCd)[2]]))

# pca matrix plot

#color

# assign color to group
Diagnosis_pca_color_QCd<-labels2colors(as.character(brain_data_QCd$Diagnosis))

# pca plot - color by disease - case/control
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Disease Status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomleft', as.character(unique(brain_data_QCd$Diagnosis)), fill=unique(Diagnosis_pca_color_QCd))

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
legend('bottomleft', as.character(unique(brain_region_pca_QCd)), fill=unique(brain_region_pca_color_QCd))

#plot to pdf - before/after QC

setwd(pca_dir)

pdf("PCA_plot_prior_and_after_QC.pdf")
# pca plot by diagnosis
plot(pca$rotation[,1:2], main=" PCA plot coloured by disease status ",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomleft', as.character(unique(Diagnosis_pca)), fill=unique(Diagnosis_pca_color))
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by disease status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomleft', as.character(unique(brain_data_QCd$Diagnosis)), fill=unique(Diagnosis_pca_color_QCd))
#pca plot by brain region
plot(pca$rotation[,1:2], main=" PCA plot coloured by brain region ",col="black", pch=21,bg=brain_region_pca_color)
legend('bottomleft', as.character(unique(brain_region_pca)), fill=unique(brain_region_pca_color))
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by brain region after QC",col="black", pch=21,bg=brain_region_pca_color_QCd)
legend('bottomleft', as.character(unique(brain_region_pca_QCd)), fill=unique(brain_region_pca_color_QCd), cex=0.75)
dev.off()

setwd(work_dir)

##### CONVERT PROBE ID TO ENTREZ ID #####

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using hgu133a.db
mapped_probes <- mappedkeys(hgu133plus2ENTREZID)

# Convert to a list
hgu133plus2.db_mapping <- as.data.frame(hgu133plus2ENTREZID[mapped_probes])
# arrange order of column by entrezgene probe_id
hgu133plus2.db_mapping<-hgu133plus2.db_mapping[c(2,1)]
colnames(hgu133plus2.db_mapping)[1]<-"entrezgene"

head(hgu133plus2.db_mapping)
dim(hgu133plus2.db_mapping)

#check any duplicated probe IDs
anyDuplicated(hgu133plus2.db_mapping$probe_id)

#check any dupliacted entrezgene IDs
anyDuplicated(hgu133plus2.db_mapping$entrezgene)

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

# using hgu133plus2.db

brain_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(brain_data_QCd[2:dim(brain_data_QCd)[2]], hgu133plus2.db_mapping)
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

rownames(Diagnosis)==rownames(brain_region)

Diagnosis_and_tissue_type_lookup<-cbind(brain_region, Diagnosis)

colnames(Diagnosis_and_tissue_type_lookup)<-c("Tissue", "Diagnosis")

brain_data_QCd_entrez_id_unique_pheno<-merge(Diagnosis_and_tissue_type_lookup, brain_data_QCd_entrez_id_unique, by="row.names")

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

rownames(brain_data_QCd_entrez_id_unique_pheno)<-brain_data_QCd_entrez_id_unique_pheno$Row.names

brain_data_QCd_entrez_id_unique_pheno$Row.names<-NULL

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

dim(brain_data_QCd_entrez_id_unique_pheno)

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(brain_data_QCd_entrez_id_unique_pheno, file="E-GEOD-5281_pre-processed_data2.txt", sep="\t")

##### SAVE PHENOTYPE DATAFRAME #####

setwd(clean_data_dir)

write.table(phenotype_data, file="E-GEOD-5281_phenotype_data2.txt", sep="\t")

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="E-GEOD-5281_data_processing2.Rdata")
