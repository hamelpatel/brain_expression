
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

work_dir="/media/hamel/1TB/Projects/Brain_expression/1.Data/2.Re-process_with_new_pipeline/AD/AMP_MAyoeGWAS/Pre-processing"

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

##### PROBE ID DETECTION #####

# probe filtering already done to a level of 75% detection

##### PCA ####

# TEMPORAL CORTEX

# calculate pca

Temporal_cortex_pca<-prcomp(Temporal_cortex_exprs_dataframe)

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
Temporal_cortex_Diagnosis_temp<-colnames(Temporal_cortex_exprs_dataframe)

# match order
Temporal_cortex_Diagnosis_pca<-Temporal_cortex_Diagnosis[match(Temporal_cortex_Diagnosis_temp, rownames(Temporal_cortex_Diagnosis)),]
head(Temporal_cortex_Diagnosis_pca)

# assign color to group
Temporal_cortex_Diagnosis_pca_color<-labels2colors(as.character(Temporal_cortex_Diagnosis_pca))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Temporal_cortex_Diagnosis_pca_color)
legend('bottomleft', unique(Temporal_cortex_Diagnosis_pca), fill=unique(Temporal_cortex_Diagnosis_pca_color))


# CEREBEULLUM

# calculate pca

Cerebellum_pca<-prcomp(Cerebellum_exprs_dataframe)

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
Cerebellum_Diagnosis_temp<-colnames(Cerebellum_exprs_dataframe)

# match order
Cerebellum_Diagnosis_pca<-Cerebellum_Diagnosis[match(Cerebellum_Diagnosis_temp, rownames(Cerebellum_Diagnosis)),]
head(Cerebellum_Diagnosis_pca)

# assign color to group
Cerebellum_Diagnosis_pca_color<-labels2colors(as.character(Cerebellum_Diagnosis_pca))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Cerebellum PCA plot coloured by Disease Status",col="black", pch=21,bg=Cerebellum_Diagnosis_pca_color)
legend('bottomleft', unique(Cerebellum_Diagnosis_pca), fill=unique(Cerebellum_Diagnosis_pca_color))

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
rownames(Temporal_cortex_test)==colnames(Temporal_cortex_exprs_dataframe)

pdf("Temporal_Cortex_batch_effect_PCA_plot.pdf")

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Temporal_cortex_Diagnosis_pca_color)
legend('bottomleft', unique(Temporal_cortex_Diagnosis_pca), fill=unique(Temporal_cortex_Diagnosis_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,5]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_test[,5])), fill=unique(Temporal_cortex_test_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,6]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_test[,6])), fill=unique(Temporal_cortex_test_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,7]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_test[,7])), fill=unique(Temporal_cortex_test_pca_color))


# assign color to group
Temporal_cortex_test_pca_color<-labels2colors(as.character(Temporal_cortex_test[,8]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_test[,8])), fill=unique(Temporal_cortex_test_pca_color))

dev.off()


#
#
# Cerebellum
#
#


# match order
Cerebellum_test<-Cerebellum_pheno[match(Cerebellum_Diagnosis_temp, rownames(Cerebellum_pheno)),]
head(Cerebellum_test)
rownames(Cerebellum_test)==colnames(Cerebellum_exprs_dataframe)

pdf("Cerebellum_batch_effect_PCA_plot.pdf")

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Disease Status",col="black", pch=21,bg=Cerebellum_Diagnosis_pca_color)
legend('bottomleft', unique(Cerebellum_Diagnosis_pca), fill=unique(Cerebellum_Diagnosis_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,5]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate0 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_test[,5])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,6]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate1 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_test[,6])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,7]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_test[,7])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,8]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_test[,8])), fill=unique(Cerebellum_test_pca_color))


# assign color to group
Cerebellum_test_pca_color<-labels2colors(as.character(Cerebellum_test[,9]))

# pca plot - color by disease - case/control
plot(Cerebellum_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Cerebellum_test_pca_color)
legend('bottomleft', as.character(unique(Cerebellum_test[,9])), fill=unique(Cerebellum_test_pca_color))

dev.off()

##### SVA ####

head(Temporal_cortex_exprs)[1:5]

head(Temporal_cortex_pheno)

# convert dataframe to matrix
Temporal_cortex_exprs_matrix<-as.matrix(Temporal_cortex_exprs_dataframe)

#full model matrix for Diagnosis
mod = model.matrix(~Diagnosis, data=Temporal_cortex_pheno)
mod0 = model.matrix(~1, data=Temporal_cortex_pheno)

mod
mod0

# check number of SV in data
n.sv<-num.sv(Temporal_cortex_exprs_matrix, mod, method="leek")
n.sv

# apply sva - leave n.sv as NULL
svobj = sva(Temporal_cortex_exprs_matrix, mod, mod0, method="two-step")

# adjust for sva in dataframe
X = cbind(mod, svobj$sv)
Hat = solve(t(X) %*% X) %*% t(X)
beta = (Hat %*% t(Temporal_cortex_exprs_matrix))
P = ncol(mod)
Temporal_cortex_clean_data<-Temporal_cortex_exprs_matrix - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])

#check SV in clean data - shoould be zero
num.sv(Temporal_cortex_clean_data, mod, method="leek")

###
#
# Cerebellum
#
###

head(Cerebellum_exprs)[1:5]

head(Cerebellum_pheno)

# convert dataframe to matrix
Cerebellum_exprs_matrix<-as.matrix(Cerebellum_exprs_dataframe)

#full model matrix for Diagnosis
mod = model.matrix(~Diagnosis, data=Cerebellum_pheno)
mod0 = model.matrix(~1, data=Cerebellum_pheno)

mod
mod0

# check number of SV in data
n.sv<-num.sv(Cerebellum_exprs_matrix, mod, method="leek")
n.sv

# apply sva - leave n.sv as NULL
svobj = sva(Cerebellum_exprs_matrix, mod, mod0, method="two-step")

# adjust for sva in dataframe
X = cbind(mod, svobj$sv)
Hat = solve(t(X) %*% X) %*% t(X)
beta = (Hat %*% t(Cerebellum_exprs_matrix))
P = ncol(mod)
Cerebellum_clean_data<-Cerebellum_exprs_matrix - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])

#check SV in clean data - shoould be zero
num.sv(Cerebellum_clean_data, mod, method="leek")

##### PLOTS AFTER SVA #####

# plot density plot on clean data
plotDensity(Temporal_cortex_clean_data,logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )

# plot PCA plot on clean data

Temporal_cortex_clean_data_dataframe<-as.data.frame(Temporal_cortex_clean_data)
head(Temporal_cortex_clean_data_dataframe)[1:5]

Temporal_cortex_SVA_pca<-prcomp(Temporal_cortex_clean_data_dataframe)

# plot variance
plot(Temporal_cortex_SVA_pca, type="l")

#color

head(Temporal_cortex_pheno)

Temporal_cortex_SVA_Diagnosis<-Temporal_cortex_pheno[1]
head(Temporal_cortex_SVA_Diagnosis)
table(Temporal_cortex_SVA_Diagnosis)

# order of samples in expression data
Temporal_cortex_SVA_Diagnosis_temp<-colnames(Temporal_cortex_clean_data_dataframe)

# match order
Temporal_cortex_SVA_test<-Temporal_cortex_pheno[match(Temporal_cortex_SVA_Diagnosis_temp, rownames(Temporal_cortex_pheno)),]
head(Temporal_cortex_SVA_test)
rownames(Temporal_cortex_SVA_test)==colnames(Temporal_cortex_clean_data_dataframe)

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,1]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,1])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,5]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,5])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,6]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,6])), fill=unique(Temporal_cortex_SVA_test_pca_color))

# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,7]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,7])), fill=unique(Temporal_cortex_SVA_test_pca_color))


# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,8]))

# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,8])), fill=unique(Temporal_cortex_SVA_test_pca_color))


#plot to pdf
setwd(pca_dir)
pdf("Temporal_cortex_plots_after_SVA_adjustment.pdf")
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,1]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by Diagnosis samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,1])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate2 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,5])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,6]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate3 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,6])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,7]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate4 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,7])), fill=unique(Temporal_cortex_SVA_test_pca_color))
# assign color to group
Temporal_cortex_SVA_test_pca_color<-labels2colors(as.character(Temporal_cortex_SVA_test[,8]))
# pca plot - color by disease - case/control
plot(Temporal_cortex_SVA_pca$rotation[,1:2], main="Temporal cortex PCA plot coloured by plate5 samples",col="black", pch=21,bg=Temporal_cortex_SVA_test_pca_color)
legend('bottomleft', as.character(unique(Temporal_cortex_SVA_test[,8])), fill=unique(Temporal_cortex_SVA_test_pca_color))
dev.off()

#plot density plot
setwd(work_dir)
pdf("Temporal_Cortex_Density_plot_after_SVA.pdf")
plotDensity(Temporal_cortex_clean_data,logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )
dev.off()

###
#
# Cerebellum
#
###

# plot density plot on clean data
plotDensity(Cerebellum_clean_data,logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )

# plot PCA plot on clean data

Cerebellum_clean_data_dataframe<-as.data.frame(Cerebellum_clean_data)
head(Cerebellum_clean_data_dataframe)[1:5]

Cerebellum_SVA_pca<-prcomp(Cerebellum_clean_data_dataframe)

# plot variance
plot(Cerebellum_SVA_pca, type="l")

#color

head(Cerebellum_pheno)

Cerebellum_SVA_Diagnosis<-Cerebellum_pheno[1]
head(Cerebellum_SVA_Diagnosis)
table(Cerebellum_SVA_Diagnosis)

# order of samples in expression data
Cerebellum_SVA_Diagnosis_temp<-colnames(Cerebellum_clean_data_dataframe)

# match order
Cerebellum_SVA_test<-Cerebellum_pheno[match(Cerebellum_SVA_Diagnosis_temp, rownames(Cerebellum_pheno)),]
head(Cerebellum_SVA_test)
rownames(Cerebellum_SVA_test)==colnames(Cerebellum_clean_data_dataframe)

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
dev.off()

#plot density plot
setwd(work_dir)
pdf("Cerebellum_Density_plot_after_SVA.pdf")
plotDensity(Cerebellum_clean_data,logMode=F, addLegend=F, main="Density plot of Temporal Cortex Samples after SVA" )
dev.off()

##### NETWORK ANALYSIS #####

#separate case and control ID's

Temporal_cortex_case_IDs<-rownames(subset(Temporal_cortex_Diagnosis, Diagnosis=="AD"))

Temporal_cortex_control_IDs<-rownames(subset(Temporal_cortex_Diagnosis, Diagnosis=="Control"))

length(Temporal_cortex_case_IDs)
length(Temporal_cortex_control_IDs)
table(Temporal_cortex_Diagnosis$Diagnosis)


# order Diagnosis by order of expression table - from pca plot

head(Temporal_cortex_clean_data_dataframe)[1:5]

# match order to 
rownames(Temporal_cortex_Diagnosis)==colnames(Temporal_cortex_clean_data_dataframe)

# sample plot function - taken from steve expression pipeline

sampleNetwork_plot <- function(datExprs, diagnosis, colBy=c("chip","group") ) {
  gp_col <- colBy
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
  ## ClusterCoef
  #par(mar=c(5,5,4,2))
  #plot(Z.C,main="ClusterCoef", ylab="Z.C",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4,type="n")
  #text(Z.C,labels=samle_names,cex=0.8,col=colorvec)
  #abline(h=-2)
  #abline(h=-3)
  ## Connectivity vs ClusterCoef
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

names_of_outliers<-function(datExprs, threshold){
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

# run plot on full data - by Diease

sampleNetwork_plot(Temporal_cortex_clean_data_dataframe, Temporal_cortex_Diagnosis, colBy="group")

# case only

Temporal_cortex_clean_data_dataframe_cases<-Temporal_cortex_clean_data_dataframe[,colnames(Temporal_cortex_clean_data_dataframe)%in%Temporal_cortex_case_IDs]

sampleNetwork_plot(Temporal_cortex_clean_data_dataframe_cases, subset(Temporal_cortex_Diagnosis, Diagnosis=="AD"), colBy="group")

# control only

Temporal_cortex_clean_data_dataframe_controls<-Temporal_cortex_clean_data_dataframe[,colnames(Temporal_cortex_clean_data_dataframe)%in%Temporal_cortex_control_IDs]

sampleNetwork_plot(Temporal_cortex_clean_data_dataframe_controls, subset(Temporal_cortex_Diagnosis, Diagnosis=="Control"), colBy="group")


# create function to run network analysis on each expression dataset, plot and remove bad samples - used Z.K value of -3 as ISA ~1

run_sample_network_plot<-function(dataset, Diagnosis){
  #sample network plot
  sampleNetwork_plot(dataset, Diagnosis, colBy="group")
  #identify sample below Z.K -2  
  dataset_removal_1<-names_of_outliers(dataset, -3)
  #remove samples with ZK below -2
  dataset_QC<-dataset[,!(colnames(dataset)%in%dataset_removal_1)]
  #sample network plot
  sampleNetwork_plot(dataset_QC, Diagnosis, colBy="group")
  #create empty count list to record samples removed
  count<-dataset_removal_1
  # reiterate above till no samples fall below -2
  while (length(dataset_removal_1)>0) {
    # remove bad samples - 1st iteration removes none
    dataset_QC<-dataset_QC[,!(colnames(dataset_QC)%in%dataset_removal_1)]
    #identify sample below Z.K -2  
    dataset_removal_1<-names_of_outliers(dataset_QC, -3)
    #record samples removed
    count<-c(count, dataset_removal_1)
  }
  #final network plot
  sampleNetwork_plot(dataset_QC, Diagnosis, colBy="group")
  # print to screen number of samples removed
  cat("\n")
  print(c("Total number of samples removed...", length(count)))
  # return clean expression set
  return(dataset_QC)
}

# run sample plot on individual brain region by case and exclude outliers

Temporal_cortex_clean_data_dataframe_cases_QC<-run_sample_network_plot(Temporal_cortex_clean_data_dataframe_cases, subset(Temporal_cortex_Diagnosis, Diagnosis=="AD"))
dim(Temporal_cortex_clean_data_dataframe_cases)
dim(Temporal_cortex_clean_data_dataframe_cases_QC)

Temporal_cortex_clean_data_dataframe_controls_QC<-run_sample_network_plot(Temporal_cortex_clean_data_dataframe_controls, subset(Temporal_cortex_Diagnosis, Diagnosis=="Control"))
dim(Temporal_cortex_clean_data_dataframe_controls)
dim(Temporal_cortex_clean_data_dataframe_controls_QC)

##
#
# Cerebellum
#
##


#separate case and control ID's

Cerebellum_case_IDs<-rownames(subset(Cerebellum_Diagnosis, Diagnosis=="AD"))

Cerebellum_control_IDs<-rownames(subset(Cerebellum_Diagnosis, Diagnosis=="Control"))

length(Cerebellum_case_IDs)
length(Cerebellum_control_IDs)
table(Cerebellum_Diagnosis$Diagnosis)

# order Diagnosis by order of expression table - from pca plot

head(Cerebellum_clean_data_dataframe)[1:5]

# match order to 
rownames(Cerebellum_Diagnosis)==colnames(Cerebellum_clean_data_dataframe)

# run plot on full data - by Diease

sampleNetwork_plot(Cerebellum_clean_data_dataframe, Cerebellum_Diagnosis, colBy="group")

# case only

Cerebellum_clean_data_dataframe_cases<-Cerebellum_clean_data_dataframe[,colnames(Cerebellum_clean_data_dataframe)%in%Cerebellum_case_IDs]

sampleNetwork_plot(Cerebellum_clean_data_dataframe_cases, subset(Cerebellum_Diagnosis, Diagnosis=="AD"), colBy="group")

# control only

Cerebellum_clean_data_dataframe_controls<-Cerebellum_clean_data_dataframe[,colnames(Cerebellum_clean_data_dataframe)%in%Cerebellum_control_IDs]

sampleNetwork_plot(Cerebellum_clean_data_dataframe_controls, subset(Cerebellum_Diagnosis, Diagnosis=="Control"), colBy="group")

# run sample plot on individual brain region by case and exclude outliers

Cerebellum_clean_data_dataframe_cases_QC<-run_sample_network_plot(Cerebellum_clean_data_dataframe_cases, subset(Cerebellum_Diagnosis, Diagnosis=="AD"))
dim(Cerebellum_clean_data_dataframe_cases)
dim(Cerebellum_clean_data_dataframe_cases_QC)

Cerebellum_clean_data_dataframe_controls_QC<-run_sample_network_plot(Cerebellum_clean_data_dataframe_controls, subset(Cerebellum_Diagnosis, Diagnosis=="Control"))
dim(Cerebellum_clean_data_dataframe_controls)
dim(Cerebellum_clean_data_dataframe_controls_QC)

# plot to PDF

setwd(sample_network_dir)

pdf("Netowrk_Anlaysis_on_Temporal_Cortex_Cases.pdf")
Temporal_cortex_clean_data_dataframe_cases_QC<-run_sample_network_plot(Temporal_cortex_clean_data_dataframe_cases, subset(Temporal_cortex_Diagnosis, Diagnosis=="AD"))
dev.off()

pdf("Netowrk_Anlaysis_on_Temporal_Cortex_Controls.pdf")
Temporal_cortex_clean_data_dataframe_controls_QC<-run_sample_network_plot(Temporal_cortex_clean_data_dataframe_controls, subset(Temporal_cortex_Diagnosis, Diagnosis=="Control"))
dev.off()

pdf("Netowrk_Anlaysis_on_Cerebellum_Cases.pdf")
Cerebellum_clean_data_dataframe_cases_QC<-run_sample_network_plot(Cerebellum_clean_data_dataframe_cases, subset(Cerebellum_Diagnosis, Diagnosis=="AD"))
dev.off()

pdf("Netowrk_Anlaysis_on_Cerebellum_Control.pdf")
Cerebellum_clean_data_dataframe_controls_QC<-run_sample_network_plot(Cerebellum_clean_data_dataframe_controls, subset(Cerebellum_Diagnosis, Diagnosis=="Control"))
dev.off()

##### CREATE QC'd DATASET #####

Temporal_cortex_good_sample_ID<-c(colnames(Temporal_cortex_clean_data_dataframe_cases_QC),
                                  colnames(Temporal_cortex_clean_data_dataframe_controls_QC))

Cerebellum_good_sample_ID<-c(colnames(Cerebellum_clean_data_dataframe_cases_QC),
                             colnames(Cerebellum_clean_data_dataframe_controls_QC))

Temporal_cortex_QCd<-Temporal_cortex_clean_data_dataframe[,colnames(Temporal_cortex_clean_data_dataframe)%in%Temporal_cortex_good_sample_ID]
dim(Temporal_cortex_clean_data_dataframe)
dim(Temporal_cortex_QCd)

Cerebellum_QCd<-Cerebellum_clean_data_dataframe[,colnames(Cerebellum_clean_data_dataframe)%in%Cerebellum_good_sample_ID]
dim(Cerebellum_clean_data_dataframe)
dim(Cerebellum_QCd)

##### CONVERT PROBE ID TO ENTREZ ID #####

#check probes same in both datasets
any(rownames(Temporal_cortex_QCd)==rownames(Cerebellum_QCd))=="FLASE"

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
  expression_dataset_t<-as.data.frame(t(expression_dataset))
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

# arrange in order of

Cerebellum_pheno_to_attach<-Cerebellum_pheno_to_attach[c(5,6,2,7,1,8,3,4)]
head(Cerebellum_pheno_to_attach)

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

# arrange in order of

Temporal_cortex_pheno_to_attach<-Temporal_cortex_pheno_to_attach[c(5,6,2,7,1,8,3,4)]
head(Temporal_cortex_pheno_to_attach)

#subset phenotype to samples in Temporal_cortex_QCd_entrez_id_unique
dim(Temporal_cortex_pheno_to_attach)
Temporal_cortex_pheno_to_attach<-subset(Temporal_cortex_pheno_to_attach, rownames(Temporal_cortex_pheno_to_attach) %in% rownames(Temporal_cortex_QCd_entrez_id_unique))
dim(Temporal_cortex_pheno_to_attach)
dim(Temporal_cortex_QCd_entrez_id_unique)

#subset phenotype to samples in Cerebellum_QCd_entrez_id_unique
dim(Cerebellum_pheno_to_attach)
Cerebellum_pheno_to_attach<-subset(Cerebellum_pheno_to_attach, rownames(Cerebellum_pheno_to_attach) %in% rownames(Cerebellum_QCd_entrez_id_unique))
dim(Cerebellum_pheno_to_attach)
dim(Cerebellum_QCd_entrez_id_unique)

# merge Temporal_cortex pheno info with expr

Temporal_cortex_QCd_entrez_id_unique_pheno<-merge(Temporal_cortex_pheno_to_attach, Temporal_cortex_QCd_entrez_id_unique, by="row.names", all.x=T)
head(Temporal_cortex_QCd_entrez_id_unique_pheno)[1:5]
rownames(Temporal_cortex_QCd_entrez_id_unique_pheno)<-Temporal_cortex_QCd_entrez_id_unique_pheno$Row.names
Temporal_cortex_QCd_entrez_id_unique_pheno$Row.names<-NULL
head(Temporal_cortex_QCd_entrez_id_unique_pheno)[1:5]
dim(Temporal_cortex_QCd_entrez_id_unique)
dim(Temporal_cortex_pheno_to_attach)
dim(Temporal_cortex_QCd_entrez_id_unique_pheno)

# merge Cerebellum pheno info with expr

Cerebellum_QCd_entrez_id_unique_pheno<-merge(Cerebellum_pheno_to_attach, Cerebellum_QCd_entrez_id_unique, by="row.names", all.x=T)
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

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file="AMP_MayoeGWAS_data_processing2.Rdata")
