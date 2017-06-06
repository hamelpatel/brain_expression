########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                                   PATHWAY ANALYSIS - SPIA                                            #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 05/04/2017

##### DESCRIPTION OF ANALYSIS ####
## SPIA is a topology pathway analysis technique.
## requires logFC - so will convert p values to mimick logFC values
##
##
## NOTE: 
## 
## 
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARIES #####

library(SPIA)

##### SET DIRECTORIES ####

data_dir1="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/4.AD_DE_Meta_Analysis/AW_results/AW_results_with_logFC/"

data_dir2="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/7.Combine_AD_and_non_AD_Meta_results/DE_lists"

RNA_SEQ_data_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/AD/RNA_Seq"

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/8.Pathway_Analysis/8.4.SPIA/"

##### READ DEG LIST #####

setwd(data_dir1)

Analysis_1_AD_Cerebellum_full_DEG<-read.table("Cerebellum_DEG_sig_AW_results_with_logFC.txt", as.is=T, header=T, sep="\t")
Analysis_1_AD_Frontal_Lobe_full_DEG<-read.table("Frontal_lobe_sig_AW_results_with_logFC.txt", as.is=T, header=T, sep="\t")
Analysis_1_AD_Parietal_Lobe_full_DEG<-read.table("Parietal_lobe_sig_AW_results_with_logFC.txt", as.is=T, header=T, sep="\t")
Analysis_1_AD_Temporal_Lobe_full_DEG<-read.table("Temporal_lobe_sig_AW_results_with_logFC.txt", as.is=T, header=T, sep="\t")

head(Analysis_1_AD_Cerebellum_full_DEG)
head(Analysis_1_AD_Frontal_Lobe_full_DEG)
head(Analysis_1_AD_Parietal_Lobe_full_DEG)
head(Analysis_1_AD_Temporal_Lobe_full_DEG)

setwd(data_dir2)
Analysis_2_AD_specific_Cerebellum<-read.table("2.AD_Cerebellum_unique.txt", as.is=T, header=T, sep="\t")
Analysis_2_AD_specific_Frontal_Lobe<-read.table("2.AD_Frontal_Lobe_unique.txt", as.is=T, header=T, sep="\t")
Analysis_2_AD_specific_Parietal_Lobe<-read.table("2.AD_Parietal_Lobe_unique.txt", as.is=T, header=T, sep="\t")

Analysis_3_AD_Cerebellum_common<-read.table("4.Cerebellum_common_DEG.txt", as.is=T, header=T, sep="\t")
Analysis_3_AD_Frontal_Lobe_common<-read.table("4.Frontal_Lobe_common_DEG.txt", as.is=T, header=T, sep="\t")
Analysis_3_AD_Parietal_Lobe_common<-read.table("4.Parietal_Lobe_common_DEG.txt", as.is=T, header=T, sep="\t")

# full DEG list for background genes

full_DEG_list<-rownames(read.table("/media/hamel/Workspace/Dropbox/Projects/Brain_expression/3.AD_single_dataset_DE_analysis/Differential_Expression/Differential_expression_on_full_expression_data_AD_vs_CONTROL.txt", header=T, as.is=T, sep="\t"))

head(full_DEG_list)
length(full_DEG_list)

# RNA SEQ DATA

setwd(RNA_SEQ_data_dir)

AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG<-read.table("Frontal_Lobe/DE_results/Full_Differential_Expression_Results.txt", head=T, as.is=T)

AD_RNA_Seq_AD_Temporal_Lobe_Full_DEG<-read.table("Temporal_Lobe/DE_results/Full_Differential_Expression_Results.txt", head=T, as.is=T)

head(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG)
head(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG)

##### STANDARDISE DATA #####

#subset analysis 2 + 3 from analysis 1
Analysis_2_AD_specific_Cerebellum<-Analysis_1_AD_Cerebellum_full_DEG[Analysis_1_AD_Cerebellum_full_DEG$Entrez_Gene_ID %in% Analysis_2_AD_specific_Cerebellum$Entrez_Gene_ID,]
Analysis_2_AD_specific_Frontal_Lobe<-Analysis_1_AD_Frontal_Lobe_full_DEG[Analysis_1_AD_Frontal_Lobe_full_DEG$Entrez_Gene_ID %in% Analysis_2_AD_specific_Frontal_Lobe$Entrez_Gene_ID,]
Analysis_2_AD_specific_Parietal_Lobe<-Analysis_1_AD_Parietal_Lobe_full_DEG[Analysis_1_AD_Parietal_Lobe_full_DEG$Entrez_Gene_ID %in% Analysis_2_AD_specific_Parietal_Lobe$Entrez_Gene_ID,]

dim(Analysis_1_AD_Cerebellum_full_DEG)
dim(Analysis_1_AD_Frontal_Lobe_full_DEG)
dim(Analysis_1_AD_Parietal_Lobe_full_DEG)

dim(Analysis_2_AD_specific_Cerebellum)
dim(Analysis_2_AD_specific_Frontal_Lobe)
dim(Analysis_2_AD_specific_Parietal_Lobe)

Analysis_3_AD_Cerebellum_common<-Analysis_1_AD_Cerebellum_full_DEG[Analysis_1_AD_Cerebellum_full_DEG$Entrez_Gene_ID %in% Analysis_3_AD_Cerebellum_common$Entrez_Gene_ID,]
Analysis_3_AD_Frontal_Lobe_common<-Analysis_1_AD_Frontal_Lobe_full_DEG[Analysis_1_AD_Frontal_Lobe_full_DEG$Entrez_Gene_ID %in% Analysis_3_AD_Frontal_Lobe_common$Entrez_Gene_ID,]
Analysis_3_AD_Parietal_Lobe_common<-Analysis_1_AD_Parietal_Lobe_full_DEG[Analysis_1_AD_Parietal_Lobe_full_DEG$Entrez_Gene_ID %in% Analysis_3_AD_Parietal_Lobe_common$Entrez_Gene_ID,]

dim(Analysis_1_AD_Cerebellum_full_DEG)
dim(Analysis_1_AD_Frontal_Lobe_full_DEG)
dim(Analysis_1_AD_Parietal_Lobe_full_DEG)

dim(Analysis_3_AD_Cerebellum_common)
dim(Analysis_3_AD_Frontal_Lobe_common)
dim(Analysis_3_AD_Parietal_Lobe_common)

#create function to manipulate data
prep_data_for_heatmap<-function(dataset){
  new_object<-dataset$meta_summary
names(new_object)<-as.vector(dataset$Entrez_Gene_ID)
  #return object
  return(new_object)
}

# apply to data analysis 1

Analysis_1_AD_Cerebellum_full_DEG_SPIA_format<-prep_data_for_heatmap(Analysis_1_AD_Cerebellum_full_DEG)
Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_format<-prep_data_for_heatmap(Analysis_1_AD_Frontal_Lobe_full_DEG)
Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_format<-prep_data_for_heatmap(Analysis_1_AD_Parietal_Lobe_full_DEG)
Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_format<-prep_data_for_heatmap(Analysis_1_AD_Temporal_Lobe_full_DEG)

# apply to analysis 2
Analysis_2_AD_specific_Cerebellum_SPIA_format<-prep_data_for_heatmap(Analysis_2_AD_specific_Cerebellum)
Analysis_2_AD_specific_Frontal_Lobe_SPIA_format<-prep_data_for_heatmap(Analysis_2_AD_specific_Frontal_Lobe)
Analysis_2_AD_specific_Parietal_Lobe_SPIA_format<-prep_data_for_heatmap(Analysis_2_AD_specific_Parietal_Lobe)

# apply to analysis 3
Analysis_3_AD_Cerebellum_common_SPIA_format<-prep_data_for_heatmap(Analysis_3_AD_Cerebellum_common)
Analysis_3_AD_Frontal_Lobe_common_SPIA_format<-prep_data_for_heatmap(Analysis_3_AD_Frontal_Lobe_common)
Analysis_3_AD_Parietal_Lobe_common_SPIA_format<-prep_data_for_heatmap(Analysis_3_AD_Parietal_Lobe_common)

# manipulate RNA-Seq data

AD_RNA_Seq_AD_Frontal_Lobe_sig_DEG<-as.vector((subset(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG, adj.P.Val<=0.05))[,1])
names(AD_RNA_Seq_AD_Frontal_Lobe_sig_DEG)<-rownames(subset(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG, adj.P.Val<=0.05))


AD_RNA_Seq_AD_Temporal_Lobe_sig_DEG<-as.vector((subset(AD_RNA_Seq_AD_Temporal_Lobe_Full_DEG, adj.P.Val<=0.05))[,1])
names(AD_RNA_Seq_AD_Temporal_Lobe_sig_DEG)<-rownames(subset(AD_RNA_Seq_AD_Temporal_Lobe_Full_DEG, adj.P.Val<=0.05))

##### RUN SPIA #####

## ANALYSIS 1

Analysis_1_AD_Cerebellum_full_DEG_SPIA_results<-spia(de=Analysis_1_AD_Cerebellum_full_DEG_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_results<-spia(de=Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_results<-spia(de=Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_results<-spia(de=Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)

Analysis_1_AD_RNA_Seq_Frontal_Lobe_full_DEG_SPIA_results<-spia(de=AD_RNA_Seq_AD_Frontal_Lobe_sig_DEG, all=rownames(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG), organism="hsa", plots=F)
Analysis_1_AD_RNA_Seq_Temporal_Lobe_full_DEG_SPIA_results<-spia(de=AD_RNA_Seq_AD_Temporal_Lobe_sig_DEG, all=rownames(AD_RNA_Seq_AD_Temporal_Lobe_Full_DEG), organism="hsa", plots=F)

head(Analysis_1_AD_Cerebellum_full_DEG_SPIA_results)
head(Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_results)
head(Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_results)
head(Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_results)

head(Analysis_1_AD_RNA_Seq_Frontal_Lobe_full_DEG_SPIA_results)
head(Analysis_1_AD_RNA_Seq_Temporal_Lobe_full_DEG_SPIA_results)

# subset to significant

Analysis_1_AD_Cerebellum_full_DEG_SPIA_sig<-subset(Analysis_1_AD_Cerebellum_full_DEG_SPIA_results, pGFWER<=0.05)
Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_sig<-subset(Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_results, pGFWER<=0.05)
Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_sig<-subset(Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_results, pGFWER<=0.05)
Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_sig<-subset(Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_results, pGFWER<=0.05)

Analysis_1_AD_RNA_Seq_Frontal_Lobe_full_DEG_SPIA_sig<-subset(Analysis_1_AD_RNA_Seq_Frontal_Lobe_full_DEG_SPIA_results, pGFWER<=0.05)
Analysis_1_AD_RNA_Seq_Temporal_Lobe_full_DEG_SPIA_sig<-subset(Analysis_1_AD_RNA_Seq_Temporal_Lobe_full_DEG_SPIA_results, pGFWER<=0.05)

Analysis_1_AD_Cerebellum_full_DEG_SPIA_sig[c(1,11)]
Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_sig[c(1,11)]
Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_sig[c(1,11)]
Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_sig[c(1,11)]

Analysis_1_AD_RNA_Seq_Frontal_Lobe_full_DEG_SPIA_sig[1]
Analysis_1_AD_RNA_Seq_Temporal_Lobe_full_DEG_SPIA_sig[1]

# overlap across all brain regions

Analysis_1_microarray_overlap<-rbind(Analysis_1_AD_Cerebellum_full_DEG_SPIA_sig,
                          Analysis_1_AD_Frontal_Lobe_full_DEG_SPIA_sig,
                          Analysis_1_AD_Parietal_Lobe_full_DEG_SPIA_sig,
                          Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_sig)

Analysis_1_microarray_overlap_table<-as.data.frame(table(Analysis_1_microarray_overlap$Name, Analysis_1_microarray_overlap$Status))
Analysis_1_microarray_overlap_table<-Analysis_1_microarray_overlap_table[order(-Analysis_1_microarray_overlap_table$Freq),]

head(Analysis_1_microarray_overlap_table)

as.data.frame(table(subset(Analysis_1_microarray_overlap_table, Freq!=0)[1]))

table(Analysis_1_microarray_overlap_table$Var1)

# number oinhibited + activaed
table(Analysis_1_microarray_overlap$Status)


subset(Analysis_1_microarray_overlap, Analysis_1_microarray_overlap$Name=="Alzheimer's disease")
subset(Analysis_1_microarray_overlap, Analysis_1_microarray_overlap$Name=="Mineral absorption")

table(Analysis_1_microarray_overlap[c(1, 11)])

## ANALYSIS 2


Analysis_2_AD_specific_Cerebellum_SPIA_results<-spia(de=Analysis_2_AD_specific_Cerebellum_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_2_AD_specific_Frontal_Lobe_SPIA_results<-spia(de=Analysis_2_AD_specific_Frontal_Lobe_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_2_AD_specific_Parietal_Lobe_SPIA_results<-spia(de=Analysis_2_AD_specific_Parietal_Lobe_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)

head(Analysis_2_AD_specific_Cerebellum_SPIA_results)
head(Analysis_2_AD_specific_Frontal_Lobe_SPIA_results)
head(Analysis_2_AD_specific_Parietal_Lobe_SPIA_results)

# subset to significant

Analysis_2_AD_specific_Cerebellum_SPIA_sig<-subset(Analysis_2_AD_specific_Cerebellum_SPIA_results, pGFWER<=0.05)
Analysis_2_AD_specific_Frontal_Lobe_SPIA_sig<-subset(Analysis_2_AD_specific_Frontal_Lobe_SPIA_results, pGFWER<=0.05)
Analysis_2_AD_specific_Parietal_Lobe_SPIA_sig<-subset(Analysis_2_AD_specific_Parietal_Lobe_SPIA_results, pGFWER<=0.05)

Analysis_2_AD_specific_Cerebellum_SPIA_sig[c(1,11)]
Analysis_2_AD_specific_Frontal_Lobe_SPIA_sig[c(1,11)]
Analysis_2_AD_specific_Parietal_Lobe_SPIA_sig[c(1,11)]

# overlap across all brain regions

Analysis_2_microarray_overlap<-rbind(Analysis_2_AD_specific_Cerebellum_SPIA_sig,
                                     Analysis_2_AD_specific_Frontal_Lobe_SPIA_sig,
                                     Analysis_2_AD_specific_Parietal_Lobe_SPIA_sig,
                                     Analysis_1_AD_Temporal_Lobe_full_DEG_SPIA_sig)

Analysis_2_microarray_overlap_table<-as.data.frame(table(Analysis_2_microarray_overlap$Name, Analysis_2_microarray_overlap$Status))
Analysis_2_microarray_overlap_table<-Analysis_2_microarray_overlap_table[order(-Analysis_2_microarray_overlap_table$Freq),]

head(Analysis_2_microarray_overlap_table)

#pathways lost in anlsysis 2

unique(Analysis_2_microarray_overlap_table[1])
unique(Analysis_1_microarray_overlap_table[1])

subset(Analysis_1_microarray_overlap, !(Name %in% Analysis_2_microarray_overlap$Name))[1]

## ANALYSIS 3

Analysis_3_AD_Cerebellum_common_SPIA_results<-spia(de=Analysis_3_AD_Cerebellum_common_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_3_AD_Frontal_Lobe_SPIA_results<-spia(de=Analysis_3_AD_Frontal_Lobe_common_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)
Analysis_3_AD_Parietal_Lobe_common_SPIA_results<-spia(de=Analysis_3_AD_Parietal_Lobe_common_SPIA_format, all=full_DEG_list, organism="hsa", plots=F)

head(Analysis_3_AD_Cerebellum_common_SPIA_results)
head(Analysis_3_AD_Frontal_Lobe_SPIA_results)
head(Analysis_3_AD_Parietal_Lobe_common_SPIA_results)

# subset to significant

Analysis_3_AD_Cerebellum_common_SPIA_sig<-subset(Analysis_3_AD_Cerebellum_common_SPIA_results, pGFWER<=0.05)
Analysis_3_AD_Frontal_Lobe_SPIA_sig<-subset(Analysis_3_AD_Frontal_Lobe_SPIA_results, pGFWER<=0.05)
Analysis_3_AD_Parietal_Lobe_common_SPIA_sig<-subset(Analysis_3_AD_Parietal_Lobe_common_SPIA_results, pGFWER<=0.05)

Analysis_3_AD_Cerebellum_common_SPIA_sig[c(1,11)]
Analysis_3_AD_Frontal_Lobe_SPIA_sig[c(1,11)]
Analysis_3_AD_Parietal_Lobe_common_SPIA_sig[c(1,11)]

##### SAVE IMAGE #####

setwd(work_dir)

#write out list of background genes
write.table(as.data.frame(full_DEG_list), file="full_background_list.txt", row.names=F, quote=F, col.names=F)
write.table(as.data.frame(rownames(AD_RNA_Seq_AD_Frontal_Lobe_Full_DEG)), file="RNA_seq_Frontal_lobe_full_background_list.txt", row.names=F, quote=F, col.names=F)
write.table(as.data.frame(rownames(AD_RNA_Seq_AD_Temporal_Lobe_Full_DEG)), file="RNA_seq_Temporal_lobe_full_background_list.txt", row.names=F, quote=F, col.names=F)


save.image("PATHWAY_ANALYSIS_USING_SPIA.Rdata")
