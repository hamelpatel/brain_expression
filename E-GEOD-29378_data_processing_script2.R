
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                     AD - E-GEOD-29378  - PIEPLINE V02                                                  #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Illumina
# EXPRESSION CHIP - HT 12v3
# NUMBER OF SAMPLES - 63
# BRAIN REGIONS - Hippocampus CA1 + CA2
# 
#
# REPROCESSING DATA USING NEW PIPELINE
#

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

raw_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/E-GEOD-29378/Raw_Data"

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/E-GEOD-29378/Pre-Processing"

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
library(illuminaHumanv3.db)
library(ggplot2)
library(reshape)
library(massiR)

##### DOWNLOAD RAW DATA #####

setwd(raw_dir)

brain_data_raw=getAE("E-GEOD-29378", type = "full")

##### CREATE R EXPRESSION OBJECT FROM RAW DATA #####

#convert MAGE-TAB files into expresssion set - using processed data

cnames=getcolproc(brain_data_raw)

cnames

brain_data=procset(brain_data_raw, cnames[2])

brain_data

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

# need to log2 transform

brain_data_normalised<-log2(exprs(brain_data))

brain_data_normalised_as_data_frame<-as.data.frame(brain_data_normalised)

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

table(phenotype_data$Comment..Sample_description.)

table(phenotype_data$Comment..Sample_source_name.)

names(phenotype_data)

Diagnosis<-phenotype_data[2]

colnames(Diagnosis)<-"Diagnosis"

Diagnosis

#rename case control status

Diagnosis[grep("Control", Diagnosis$Diagnosis), ]<-"CONTROL"

Diagnosis[grep("AD", Diagnosis$Diagnosis), ]<-"AD"

Diagnosis

Diagnosis_case<-subset(Diagnosis, grepl("AD", Diagnosis))

Diagnosis_control<-subset(Diagnosis, grepl("CONTROL", Diagnosis))

Diagnosis_case

Diagnosis_control

##### SAVE PROBE IDs #####

probe_ids<-rownames(brain_data_normalised)

head(probe_ids)
length(probe_ids)

setwd(work_dir)

write(file = "E-GEOD-29378_probe_IDs.txt", probe_ids, sep="\t")

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(phenotype_data))

##### GENDER CHECK ##### 

head(phenotype_data)
names(phenotype_data)

# gender_information
gender_info<-phenotype_data[18]
colnames(gender_info)<-"Gender"
table(gender_info$Gender)
head(gender_info)

#separae male/female IDs
female_samples<-subset(gender_info, Gender=="female")
male_samples<-subset(gender_info, Gender=="male")

head(female_samples)
head(male_samples)

head(brain_data_normalised_as_data_frame)[1:5]

# get Y choromosome genes
data(y.probes)
names(y.probes)

y_chromo_probes <- data.frame(y.probes["illumina_humanht_12"])

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

# gender_info$Gender<-as.character(gender_info$Gender)
# gender_info[gender_info$Gender=="M",]<-"male"
# gender_info[gender_info$Gender=="F",]<-"female"

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

# separate case control 

case_ID<-rownames(Diagnosis_case)

control_ID<-rownames(Diagnosis_control)

case_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%case_ID]
control_exprs<-brain_data_normalised_as_data_frame[,colnames(brain_data_normalised_as_data_frame)%in%control_ID]

head(case_exprs[1:5])
dim(case_exprs)
head(control_exprs[1:5])
dim(control_exprs)

# separate by brain region - FactorValue..ORGANISM.PART. - Frontal, Cortex Hippocampus_CA2, Temporal Cortex

table(phenotype_data$Comment..Sample_source_name.)

Hippocampus_CA1<-rownames(phenotype_data[phenotype_data$Comment..Sample_source_name.=="hippocampus - CA1",])
Hippocampus_CA2<-rownames(phenotype_data[phenotype_data$Comment..Sample_source_name.=="hippocampus - CA3",])

Hippocampus_CA1
Hippocampus_CA2


# separate by brain region and disorder

Hippocampus_CA1_case_exprs<-case_exprs[,colnames(case_exprs)%in%Hippocampus_CA1]
Hippocampus_CA2_case_exprs<-case_exprs[,colnames(case_exprs)%in%Hippocampus_CA2]

Hippocampus_CA1_control_exprs<-control_exprs[,colnames(control_exprs)%in%Hippocampus_CA1]
Hippocampus_CA2_control_exprs<-control_exprs[,colnames(control_exprs)%in%Hippocampus_CA2]

# check dataframe

dim(Hippocampus_CA1_case_exprs)
dim(Hippocampus_CA2_case_exprs)

head(Hippocampus_CA1_case_exprs)[1:5]
head(Hippocampus_CA2_case_exprs)[1:5]

dim(Hippocampus_CA1_control_exprs)
dim(Hippocampus_CA2_control_exprs)

head(Hippocampus_CA1_control_exprs)[1:5]
head(Hippocampus_CA2_control_exprs)[1:5]

#separate by gender

Hippocampus_CA1_control_exprs_F<-Hippocampus_CA1_control_exprs[colnames(Hippocampus_CA1_control_exprs)%in%rownames(female_samples)]
Hippocampus_CA2_control_exprs_F<-Hippocampus_CA2_control_exprs[colnames(Hippocampus_CA2_control_exprs)%in%rownames(female_samples)]

Hippocampus_CA1_case_exprs_F<-Hippocampus_CA1_case_exprs[colnames(Hippocampus_CA1_case_exprs)%in%rownames(female_samples)]
Hippocampus_CA2_case_exprs_F<-Hippocampus_CA2_case_exprs[colnames(Hippocampus_CA2_case_exprs)%in%rownames(female_samples)]

Hippocampus_CA1_control_exprs_M<-Hippocampus_CA1_control_exprs[colnames(Hippocampus_CA1_control_exprs)%in%rownames(male_samples)]
Hippocampus_CA2_control_exprs_M<-Hippocampus_CA2_control_exprs[colnames(Hippocampus_CA2_control_exprs)%in%rownames(male_samples)]

Hippocampus_CA1_case_exprs_M<-Hippocampus_CA1_case_exprs[colnames(Hippocampus_CA1_case_exprs)%in%rownames(male_samples)]
Hippocampus_CA2_case_exprs_M<-Hippocampus_CA2_case_exprs[colnames(Hippocampus_CA2_case_exprs)%in%rownames(male_samples)]

# calculate 90th percentile for each sample in each group - dataset quantile normalised - alll will have same 90th percentile cut-off

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

Hippocampus_CA1_case_expressed_probes_list_F<-extract_good_probe_list(Hippocampus_CA1_case_exprs_F, 0.9, 0.8)
length(Hippocampus_CA1_case_expressed_probes_list_F)

Hippocampus_CA2_case_expressed_probes_list_F<-extract_good_probe_list(Hippocampus_CA2_case_exprs_F, 0.9, 0.8)
length(Hippocampus_CA2_case_expressed_probes_list_F)

Hippocampus_CA1_control_expressed_probes_list_F<-extract_good_probe_list(Hippocampus_CA1_control_exprs_F, 0.9, 0.8)
length(Hippocampus_CA1_control_expressed_probes_list_F)

Hippocampus_CA2_control_expressed_probes_list_F<-extract_good_probe_list(Hippocampus_CA2_control_exprs_F, 0.9, 0.8)
length(Hippocampus_CA2_control_expressed_probes_list_F)

Hippocampus_CA1_case_expressed_probes_list_M<-extract_good_probe_list(Hippocampus_CA1_case_exprs_M, 0.9, 0.8)
length(Hippocampus_CA1_case_expressed_probes_list_M)

Hippocampus_CA2_case_expressed_probes_list_M<-extract_good_probe_list(Hippocampus_CA2_case_exprs_M, 0.9, 0.8)
length(Hippocampus_CA2_case_expressed_probes_list_M)

Hippocampus_CA1_control_expressed_probes_list_M<-extract_good_probe_list(Hippocampus_CA1_control_exprs_M, 0.9, 0.8)
length(Hippocampus_CA1_control_expressed_probes_list_M)

Hippocampus_CA2_control_expressed_probes_list_M<-extract_good_probe_list(Hippocampus_CA2_control_exprs_M, 0.9, 0.8)
length(Hippocampus_CA2_control_expressed_probes_list_M)

# merge list of good probes from both case + control + and all brain regions, sort and keep unique values

good_probe_list<-unique(sort(c(Hippocampus_CA1_control_expressed_probes_list_F,
                               Hippocampus_CA2_control_expressed_probes_list_F,
                               Hippocampus_CA1_case_expressed_probes_list_F,
                               Hippocampus_CA2_case_expressed_probes_list_F,
                               Hippocampus_CA1_control_expressed_probes_list_M,
                               Hippocampus_CA2_control_expressed_probes_list_M,
                               Hippocampus_CA1_case_expressed_probes_list_M,
                               Hippocampus_CA2_case_expressed_probes_list_M)))

length(good_probe_list)

# extract good probes from dataset

brain_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]

Hippocampus_CA1_case_exprs_good_probes<-Hippocampus_CA1_case_exprs[rownames(Hippocampus_CA1_case_exprs)%in%good_probe_list,]

Hippocampus_CA2_case_exprs_good_probes<-Hippocampus_CA2_case_exprs[rownames(Hippocampus_CA2_case_exprs)%in%good_probe_list,]

Hippocampus_CA1_control_exprs_good_probes<-Hippocampus_CA1_control_exprs[rownames(Hippocampus_CA1_control_exprs)%in%good_probe_list,]

Hippocampus_CA2_control_exprs_good_probes<-Hippocampus_CA2_control_exprs[rownames(Hippocampus_CA2_control_exprs)%in%good_probe_list,]

# check dataframe - all should be 4948

head(Hippocampus_CA1_case_exprs_good_probes)[1:5]
dim(Hippocampus_CA1_case_exprs_good_probes)

head(Hippocampus_CA2_case_exprs_good_probes)[1:5]
dim(Hippocampus_CA2_case_exprs_good_probes)

head(Hippocampus_CA1_control_exprs_good_probes)[1:5]
dim(Hippocampus_CA1_control_exprs_good_probes)

head(Hippocampus_CA2_control_exprs_good_probes)[1:5]
dim(Hippocampus_CA2_control_exprs_good_probes)

dim(brain_data_normalised)

# create dataframe with good probes - contain all brain regions, case and control

brain_data_normalised_exprs_good_probes<-brain_data_normalised_as_data_frame[rownames(brain_data_normalised_as_data_frame)%in%good_probe_list,]
dim(brain_data_normalised_exprs_good_probes)
head(brain_data_normalised_exprs_good_probes)[1:5]

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
plot_gender_specific_genes(Hippocampus_CA1_control_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA1_control")
plot_gender_specific_genes(Hippocampus_CA2_control_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA2_control")
plot_gender_specific_genes(Hippocampus_CA1_case_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA1_case")
plot_gender_specific_genes(Hippocampus_CA2_case_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA2_case")
dev.off()

# plot gender missmatches

cbind(Gender_missmatch, Diagnosis[rownames(Gender_missmatch),])
Gender_missmatch

gender_missmatch_tissue_lookup<-phenotype_data[3]

gender_missmatch_tissue_lookup<-merge(Gender_missmatch, gender_missmatch_tissue_lookup, by="row.names")
rownames(gender_missmatch_tissue_lookup)<-gender_missmatch_tissue_lookup$Row.names
gender_missmatch_tissue_lookup$Row.names<-NULL

gender_missmatch_tissue_lookup<-merge(gender_missmatch_tissue_lookup, Diagnosis, by="row.names")
rownames(gender_missmatch_tissue_lookup)<-gender_missmatch_tissue_lookup$Row.names
gender_missmatch_tissue_lookup$Row.names<-NULL

head(gender_missmatch_tissue_lookup)
table(gender_missmatch_tissue_lookup[c(3,4)])

# plot gender missmatches

pdf("Gender_missmatch.pdf")
plot_gender_specific_genes(Hippocampus_CA1_case_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA1_case_clinical_gender")
plot_gender_specific_genes(Hippocampus_CA2_case_exprs, gender_comparison[1], gene_list, 0.9, "Hippocampus_CA2_case_clinical_gender")
plot_gender_specific_genes(Hippocampus_CA1_case_exprs, gender_comparison[2], gene_list, 0.9, "Hippocampus_CA1_case_predicted_gender")
plot_gender_specific_genes(Hippocampus_CA2_case_exprs, gender_comparison[2], gene_list, 0.9, "Hippocampus_CA2_case_predicted_gender")
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

head(Diagnosis)

head(phenotype_data)
names(phenotype_data) - Comment..Sample_source_name.

brain_region_lookup<-phenotype_data[3]
colnames(brain_region_lookup)<-"brain_region"
head(brain_region_lookup)

# order of samples in expression data
Diagnosis_temp<-colnames(brain_data_normalised_exprs_good_probes)
Diagnosis_temp

# match order
Diagnosis_pca<-Diagnosis[match(Diagnosis_temp, rownames(Diagnosis)),]

# assign color to group
Diagnosis_pca_color<-labels2colors(as.character(Diagnosis_pca))

# pca plot - color by disease - case/control
plot(pca$rotation[,1:2], main=" PCA plot coloured by Disease Status",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomright', unique(Diagnosis_pca), fill=unique(Diagnosis_pca_color))

#pca plot - color by brain region

# order of samples in expression data
brain_region_temp<-colnames(brain_data_normalised_exprs_good_probes)
brain_region_temp

# match order
brain_region_pca<-brain_region_lookup[match(brain_region_temp, rownames(brain_region_lookup)),]

# assign color to group
brain_region_pca_color<-labels2colors(as.character(brain_region_pca))

# pca plot - color by disease - case/control
plot(pca$rotation[,1:2], main=" PCA plot coloured by Brain Region",col="black", pch=21,bg=brain_region_pca_color)
legend('bottomright', unique(brain_region_pca), fill=unique(brain_region_pca_color))

##### SVA ####

#separate by brain region
head(brain_data_normalised_exprs_good_probes)[1:5]

Hippocampus_CA1_exprs_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Hippocampus_CA1]
Hippocampus_CA2_exprs_good_probes<-brain_data_normalised_exprs_good_probes[,colnames(brain_data_normalised_exprs_good_probes)%in%Hippocampus_CA2]

dim(Hippocampus_CA1_exprs_good_probes)
dim(Hippocampus_CA2_exprs_good_probes)

# add diagnosis in

Hippocampus_CA1_exprs_good_probes<-merge (Diagnosis, t(Hippocampus_CA1_exprs_good_probes), by="row.names", all.y=T)
rownames(Hippocampus_CA1_exprs_good_probes)<-Hippocampus_CA1_exprs_good_probes$Row.names
Hippocampus_CA1_exprs_good_probes$Row.names<-NULL
head(Hippocampus_CA1_exprs_good_probes)[1:5]

Hippocampus_CA2_exprs_good_probes<-merge (Diagnosis, t(Hippocampus_CA2_exprs_good_probes), by="row.names", all.y=T)
rownames(Hippocampus_CA2_exprs_good_probes)<-Hippocampus_CA2_exprs_good_probes$Row.names
Hippocampus_CA2_exprs_good_probes$Row.names<-NULL
head(Hippocampus_CA2_exprs_good_probes)[1:5]

# add gender in 

Hippocampus_CA1_exprs_good_probes<-merge(gender_comparison[1], Hippocampus_CA1_exprs_good_probes, by="row.names")
rownames(Hippocampus_CA1_exprs_good_probes)<-Hippocampus_CA1_exprs_good_probes$Row.names
Hippocampus_CA1_exprs_good_probes$Row.names<-NULL

Hippocampus_CA2_exprs_good_probes<-merge(gender_comparison[1], Hippocampus_CA2_exprs_good_probes, by="row.names")
rownames(Hippocampus_CA2_exprs_good_probes)<-Hippocampus_CA2_exprs_good_probes$Row.names
Hippocampus_CA2_exprs_good_probes$Row.names<-NULL

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
check_SV_in_data(Hippocampus_CA1_exprs_good_probes)
check_SV_in_data(Hippocampus_CA2_exprs_good_probes)

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

Hippocampus_CA1_exprs_good_probes_sva_adjusted<-adjust_for_sva(Hippocampus_CA1_exprs_good_probes)
Hippocampus_CA2_exprs_good_probes_sva_adjusted<-adjust_for_sva(Hippocampus_CA2_exprs_good_probes)

# remove gender

Hippocampus_CA1_exprs_good_probes_sva_adjusted$Clinical_Gender<-NULL
Hippocampus_CA2_exprs_good_probes_sva_adjusted$Clinical_Gender<-NULL

# separate case and control

Hippocampus_CA1_good_probes_sva_adjusted_case<-Hippocampus_CA1_exprs_good_probes_sva_adjusted[Hippocampus_CA1_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Hippocampus_CA1_good_probes_sva_adjusted_control<-Hippocampus_CA1_exprs_good_probes_sva_adjusted[Hippocampus_CA1_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]

dim(Hippocampus_CA1_good_probes_sva_adjusted_case)
dim(Hippocampus_CA1_good_probes_sva_adjusted_control)

Hippocampus_CA2_good_probes_sva_adjusted_case<-Hippocampus_CA2_exprs_good_probes_sva_adjusted[Hippocampus_CA2_exprs_good_probes_sva_adjusted$Diagnosis=="AD",]
Hippocampus_CA2_good_probes_sva_adjusted_control<-Hippocampus_CA2_exprs_good_probes_sva_adjusted[Hippocampus_CA2_exprs_good_probes_sva_adjusted$Diagnosis=="CONTROL",]

dim(Hippocampus_CA2_good_probes_sva_adjusted_case)
dim(Hippocampus_CA2_good_probes_sva_adjusted_control)

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

# run sample network on Hippocampus_CA1

Hippocampus_CA1_case_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA1_good_probes_sva_adjusted_case, -3)
Hippocampus_CA1_control_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA1_good_probes_sva_adjusted_control, -3)

# run sample network on Hippocampus_CA2

Hippocampus_CA2_case_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA2_good_probes_sva_adjusted_case, -3)
Hippocampus_CA2_control_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA2_good_probes_sva_adjusted_control, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("Hippocampus_CA1_case_sample_network_analysis.pdf")
Hippocampus_CA1_case_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA1_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Hippocampus_CA1_control_sample_network_analysis.pdf")
Hippocampus_CA1_control_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA1_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Hippocampus_CA2_case_sample_network_analysis.pdf")
Hippocampus_CA2_case_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA2_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Hippocampus_CA2_control_sample_network_analysis.pdf")
Hippocampus_CA2_control_exprs_good_probes_QC<-run_sample_network_plot(Hippocampus_CA2_good_probes_sva_adjusted_control, -3)
dev.off()

##### CREATE QC'd DATASET #####

# no samples removed

Hippocampus_CA1_QCd<-Hippocampus_CA1_exprs_good_probes_sva_adjusted
Hippocampus_CA2_QCd<-Hippocampus_CA2_exprs_good_probes_sva_adjusted

Hippocampus_CA1_QCd$Diagnosis<-NULL
Hippocampus_CA2_QCd$Diagnosis<-NULL

any(colnames(Hippocampus_CA1_QCd)==colnames(Hippocampus_CA2_QCd))==F

brain_data_QCd<-rbind(Hippocampus_CA1_QCd, Hippocampus_CA2_QCd)

##### PCA ON CLEAN DATA ####

# calculate pca

pca_QCd<-prcomp(t(brain_data_QCd))

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
Diagnosis_temp_QCd<-colnames(t(brain_data_QCd))
Diagnosis_temp_QCd

# match order
Diagnosis_pca_QCd<-Diagnosis[match(Diagnosis_temp_QCd, rownames(Diagnosis)),]

# assign color to group
Diagnosis_pca_color_QCd<-labels2colors(as.character(Diagnosis_pca_QCd))

# pca plot - color by disease - case/control
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Disease Status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomleft', unique(Diagnosis_pca_QCd), fill=unique(Diagnosis_pca_color_QCd))

#pca plot - color by brain region

# order of samples in expression data
brain_region_temp_QCd<-colnames(t(brain_data_QCd))
brain_region_temp_QCd

# match order
brain_region_pca_QCd<-brain_region_lookup[match(brain_region_temp_QCd, rownames(brain_region_lookup)),]

# assign color to group
brain_region_pca_color_QCd<-labels2colors(as.character(brain_region_pca_QCd))

# pca plot - color by disease - case/control
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by Brain Region after QC",col="black", pch=21,bg=brain_region_pca_color_QCd)
legend('bottomleft', unique(brain_region_pca_QCd), fill=unique(brain_region_pca_color_QCd))

#plot to pdf - before/after QC

setwd(pca_dir)

pdf("PCA_plot_prior_and_after_QC.pdf")
# pca plot by diagnosis
plot(pca$rotation[,1:2], main=" PCA plot coloured by disease status ",col="black", pch=21,bg=Diagnosis_pca_color)
legend('bottomright', unique(Diagnosis_pca), fill=unique(Diagnosis_pca_color))
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by disease status after QC",col="black", pch=21,bg=Diagnosis_pca_color_QCd)
legend('bottomright', unique(Diagnosis_pca_QCd), fill=unique(Diagnosis_pca_color_QCd))
#pca plot by brain region
plot(pca$rotation[,1:2], main=" PCA plot coloured by brain region ",col="black", pch=21,bg=brain_region_pca_color)
legend('bottomright', unique(brain_region_pca), fill=unique(brain_region_pca_color))
plot(pca_QCd$rotation[,1:2], main=" PCA plot coloured by brain region after QC",col="black", pch=21,bg=brain_region_pca_color_QCd)
legend('bottomright', unique(brain_region_pca_QCd), fill=unique(brain_region_pca_color_QCd))
dev.off()

setwd(work_dir)

##### CONVERT PROBE ID TO ENTREZ ID #####

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using illuminaHumanv3.db
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
  # transform dataset # - removed this step
  expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}

# using illuminaHumanv3.db_mapping_entrez_id_unique instead of Illumina_HT_12_v3_probe_list_entrez_id_unique from this point onwards

brain_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(brain_data_QCd, illuminaHumanv3.db_mapping_entrez_id_unique)
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

head(phenotype_data)
names(phenotype_data)

phenotype_to_attach<-phenotype_data[c(3, 5, 7, 44)]

head(phenotype_to_attach)

# rename columns

colnames(phenotype_to_attach)<-c("Tissue", "Age", "BRAAK", "APOE")

# standardise gender columns

table(phenotype_to_attach$Gender, exclude=NULL)

# add additional columns as unknown - AGE, GENDER, BRAAK

phenotype_to_attach$MMSE<-"Unknown"
phenotype_to_attach$APOE<-"Unknown"

# phenotype_to_attach - change Tissue name
table(phenotype_to_attach$Tissue)

phenotype_to_attach[phenotype_to_attach$Tissue=="hippocampus - CA1",1]<-"Hippocampus_CA1"
phenotype_to_attach[phenotype_to_attach$Tissue=="hippocampus - CA3",1]<-"Hippocampus_CA3"

# add diagnosis

phenotype_to_attach<-merge(phenotype_to_attach, Diagnosis, by="row.names")
rownames(phenotype_to_attach)<- phenotype_to_attach$Row.names
phenotype_to_attach$Row.names<-NULL

# add gender

phenotype_to_attach<-merge(phenotype_to_attach, gender_comparison, by="row.names")
rownames(phenotype_to_attach)<-phenotype_to_attach$Row.names
phenotype_to_attach$Row.names<-NULL
head(phenotype_to_attach)

# standardise diagnosis

phenotype_to_attach[phenotype_to_attach$Diagnosis=="CONTROL",6]<-"Control"

# convetr NA to Uknown
phenotype_to_attach[is.na(phenotype_to_attach)]<-"Unknown"

#check BRAAK available for all AD samples - convert NA to "Unknown"

table(phenotype_to_attach$BRAAK, exclude=NULL)
phenotype_to_attach$Diagnosis_sub_cat<-phenotype_to_attach$Diagnosis

# for AD_Category column, change Diagnosis AD and BRAAK = 4 to Moderate AD

head(phenotype_to_attach)
phenotype_to_attach[phenotype_to_attach$Diagnosis=="AD" & phenotype_to_attach$BRAAK=="4", 9]<-"Moderate_AD"

# for AD_Category column, change Diagnosis AD and BRAAK = 5/6 to Severe AD

phenotype_to_attach[phenotype_to_attach$Diagnosis=="AD" & phenotype_to_attach$BRAAK=="5" | phenotype_to_attach$BRAAK=="6" , 9]<-"Severe_AD"

table(phenotype_to_attach$Diagnosis_sub_cat, exclude=NULL)

# order pheno information

phenotype_to_attach<-phenotype_to_attach[,order(names(phenotype_to_attach), decreasing = T)]

# should be 9 colnames - Diagnosis", "Tissue", "Age", "Clinical_Gender", "Predicted_Gender", "MMSE", "BRAAK", "APOE", "Diagnosis_sub_cat"
colnames(phenotype_to_attach)

# add phenotype to expression

brain_data_QCd_entrez_id_unique_pheno<-merge(phenotype_to_attach, brain_data_QCd_entrez_id_unique, by="row.names")

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

rownames(brain_data_QCd_entrez_id_unique_pheno)<-brain_data_QCd_entrez_id_unique_pheno$Row.names

brain_data_QCd_entrez_id_unique_pheno$Row.names<-NULL

head(brain_data_QCd_entrez_id_unique_pheno)[1:5]

dim(brain_data_QCd_entrez_id_unique_pheno)
dim(brain_data_QCd_entrez_id_unique)

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(brain_data_QCd_entrez_id_unique_pheno, file="E-GEOD-29378_pre-processed_data2.txt", sep="\t")

##### SAVE PHENOTYPE DATAFRAME #####

setwd(clean_data_dir)

write.table(phenotype_data, file="E-GEOD-29378_phenotype_data2.txt", sep="\t")

##### WRITE LIST OF CONTROL GENES EXPRESSED #####

setwd("/media/hamel/1TB/Projects/Brain_expression/2.Expression_across_brain_regions_in_control_datasets/1.Data")

write(unique(sort(c(Hippocampus_CA1_control_expressed_probes_list_F, 
                    Hippocampus_CA1_control_expressed_probes_list_M))), file="E-GEOD-29378_Hippocampus_CA1.txt")

write(unique(sort(c(Hippocampus_CA2_control_expressed_probes_list_F, 
                    Hippocampus_CA2_control_expressed_probes_list_M))), file="E-GEOD-29378_Hippocampus_CA2.txt")

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="E-GEOD-29378_data_processing2.Rdata")
