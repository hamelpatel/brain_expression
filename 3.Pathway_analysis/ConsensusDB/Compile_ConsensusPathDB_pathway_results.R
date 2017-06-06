########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                             CONSENSUSPATHDB PATHWAYS RESULTS COMPILATION                             #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 26/04/2017

##### DESCRIPTION OF ANALYSIS ####
## each brain region in gene symbol form (taken from enrichR mapping) queried in ConsensusPathDB.
## Ontology + pathway + protein complex results saved
## 
## 
## NOTE: 
## 
## 
#####

##### LIBRARIES #####

library(gplots)
library(d3heatmap)
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### DIRECTORIES #####

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/8.Pathway_Analysis/8.5.ConsensusPathDB/8.5.1.PATHWAYS/"

dir.create(paste(work_dir, "Plots", sep="/"))
plots_dir=paste(work_dir, "Plots", sep="/")

dir.create(paste(work_dir, "Result_Tables", sep="/"))
Result_Tables_dir=paste(work_dir, "Result_Tables", sep="/")

##### READ DATA  - ONLY PATHWAY ######

setwd(work_dir)

# files in directory
file_names<-list.files()

# keep only .txt files
file_names<-file_names[grep("pathways_results.tab", file_names)]

for (i in 1:length(file_names)) 
  assign(file_names[i], 
         read.table(file_names[i], head=T, sep="\t", quote="", fill=T))

ls()

#check ones with warnings
Analysis_1_AD_RNA_SEQ_Frontal_Lobe_all_DEG_protein_complex_results.tab
Analysis_1_AD_RNA_SEQ_Temporal_Lobe_all_DEG_protein_complex_results.tab
Analysis_3_AD_Frontal_Lobe_common_DEG_protein_complex_results.tab

#general check
head(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results.tab)[1:6]
head(Analysis_1_AD_Cerebellum_all_DEG_pathways_results.tab)[1:5]
head(Analysis_1_AD_Cerebellum_all_DEG_protein_complex_results.tab)[1:5]

##### SUBSET TO SIG ONLY #####

# subset to significant data to where q value<=0.05 - (FDR adjusted p value)

for (x in 1:length(file_names)) {
  # subset to FDR sig only and assign a new object name
  assign(paste(gsub(".tab", "_sig_only", file_names[x])), subset(get(file_names[x]), q.value<=0.05))
}

dim(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results.tab)
dim(Analysis_1_AD_Cerebellum_all_DEG_pathways_results.tab)
dim(Analysis_1_AD_Cerebellum_all_DEG_protein_complex_results.tab)

dim(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only)
dim(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only)
dim(Analysis_1_AD_Cerebellum_all_DEG_protein_complex_results_sig_only)

#dim each sig file
sig_file_names<-gsub(".tab", "_sig_only", file_names)

for (x in 1:length(sig_file_names)) {
  print(sig_file_names[x])
  print(dim(get(sig_file_names[x])))
}

# create seperate lists for ontology, pathway and proteins

ontotlogy_files<-sig_file_names[grep("ontology_results_sig_only", sig_file_names)]
pathway_files<-sig_file_names[grep("pathways_results_sig_only", sig_file_names)]
protein_files<-sig_file_names[grep("protein_complex_results_sig_only", sig_file_names)]

ontotlogy_files
pathway_files
protein_files

##### ANALYSIS 1.1: MICROARRAY + RNA_SEQ OVERLAP #####

#Temporal Lobe - microarray
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[3]

#Temporal Lobe - RNA-Seq
Analysis_1_AD_RNA_SEQ_Temporal_Lobe_all_DEG_pathways_results_sig_only[3]

# any microarray pathway in RNA-seq data
subset(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only, Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only$pathway %in% Analysis_1_AD_RNA_SEQ_Temporal_Lobe_all_DEG_pathways_results_sig_only$pathway) 

#Frontal Lobe - microarray
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[3]

#Frontal Lobe - RNA-Seq
Analysis_1_AD_RNA_SEQ_Frontal_Lobe_all_DEG_pathways_results_sig_only[3]

# any microarray pathway in RNA-seq data
subset(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only, Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only$pathway %in% Analysis_1_AD_RNA_SEQ_Frontal_Lobe_all_DEG_pathways_results_sig_only$pathway) 

##### ANALYSIS 1.2: MOST SIGNIFCANT PATHWAYS PER BRAIN REGION  - FULL DEG LIST #####

# Cerebellum

head(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only)[3]
dim(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only)

# Temporal_Lobe
head(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only)[3]
dim(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only)

# Parietal_Lobe
head(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only)[3]
dim(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only)

# Frontal_Lobe
head(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only)[3]
dim(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only)

##### ANALYSIS 1.2: EXPORT FULL TABLE OF PATHWAYS ACROSS BRAIN REGIONS WITH Q-value q<=0.01 #####

#create function to subset dataset to certain q Vlue and print results with specific columns

subset_dataset<-function(dataset, q_value){
  #subset by q value
  subset_of_data<-subset(dataset, q.value<=q_value)
  # create column to count number of genes matching
  subset_of_data$gene_count<-0
  for (x in 1:nrow(subset_of_data)){
    subset_of_data$gene_count[x]<-length(unlist(strsplit(as.character(subset_of_data$members_input_overlap[x]), ";")))
  }
  print(subset_of_data[c(3,2,10,4)])
  return(subset_of_data[c(3,2,10,4)])
}

Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only, 0.01)
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only, 0.01)
Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only, 0.01)
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only, 0.01)

setwd(Result_Tables_dir)

write.table(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table, "Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")
write.table(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table, "Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")
write.table(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table, "Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")
write.table(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table, "Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")

setwd(work_dir)

##### ANALYSIS 1.3: COMMON PATHWAYS ACROSS BRAIN REGION  - FULL DEG LIST #####

# merge results
Analysis_1_all_DEG_across_brain_regions<-as.data.frame(table(rbind(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[3],
                                                                   Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                   Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                   Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[3])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions<-Analysis_1_all_DEG_across_brain_regions[order(-Analysis_1_all_DEG_across_brain_regions$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions)<-1:nrow(Analysis_1_all_DEG_across_brain_regions)

head(Analysis_1_all_DEG_across_brain_regions)

subset(Analysis_1_all_DEG_across_brain_regions, Freq==4)


## repeat with only 3 brain regions - removing cerebllum

# merge results
Analysis_1_all_DEG_across_brain_regions_no_cerebellum<-as.data.frame(table(rbind(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                                 Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                                 Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[3])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions_no_cerebellum<-Analysis_1_all_DEG_across_brain_regions_no_cerebellum[order(-Analysis_1_all_DEG_across_brain_regions_no_cerebellum$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)

head(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)

subset(Analysis_1_all_DEG_across_brain_regions_no_cerebellum, Freq==3)


## repeat using q value of 0.01

Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_001<-subset(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only, q.value<=0.01)

Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_001<-subset(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only, q.value<=0.01)
Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_001<-subset(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only, q.value<=0.01)
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_001<-subset(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only, q.value<=0.01)

# merge results
Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001<-as.data.frame(table(rbind(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_001[3],
                                                                                 Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_001[3],
                                                                                 Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_001[3])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001<-Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001[order(-Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001)

head(Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001)

subset(Analysis_1_all_DEG_across_brain_regions_no_cerebellum_001, Freq==3)

##### ANALYSIS 1.3: - HEATMAP - q=0.01 #####

# using subset of results where q<=0.01

# create frunction to rearrange data for
rearrange_pathway_data_for_heatmap<-function(data, column_name){
  data<-data[c(2,3)]
  rownames(data)<-data$pathway
  data$pathway<-NULL
  colnames(data)<-column_name
  return(data)
}

Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_001, "Cerebellum")
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_001, "Temporal_Lobe")
Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_001, "parietal_Lobe")
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_001, "Frontal_Lobe")

# merge results

MyMerge <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= TRUE, all.y=T)
  rownames(df)<-df$Row.names
  df$Row.names<-NULL
  return(df)
}

Analysis_1_all_DEG_across_brain_regions_heatmap<-as.matrix(Reduce(MyMerge, list(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_heatmap,
                                                                                Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_heatmap,
                                                                                Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_heatmap,
                                                                                Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_heatmap)))


head(Analysis_1_all_DEG_across_brain_regions_heatmap)
dim(Analysis_1_all_DEG_across_brain_regions_heatmap)

# plot
setwd(plots_dir)

heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_1_all_DEG_across_all_brain_regions_pathway_heatmap.pdf", height=25, width=15)
heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()

##### ANALYSIS 1.3: - WORD CLOUD  - q0.05 #####

# create function for word cloud plot
plot_create_wordcloud<-function(data_as_vector, words_to_remove){
  text_corpus<-Corpus(VectorSource(data_as_vector))
  
  #convert everything to lowercase
  text_corpus <- tm_map(text_corpus, content_transformer(tolower))
  
  # stem words - i.e learning > learn, walked > walk
  #text_corpus <- tm_map(text_corpus, stemDocument)
  
  # remove unwanted grammer
  text_corpus <- tm_map(text_corpus, removePunctuation)
  text_corpus <- tm_map(text_corpus, removeWords, stopwords('english'))
  
  # remove unwanted owrds
  text_corpus <- tm_map(text_corpus, removeWords, words_to_remove)
  
  pal2 <- brewer.pal(8,"Dark2")
  #plot
  wordcloud(text_corpus, min.freq=1, max.words = 150, random.order = FALSE, colors=pal2)
}

# apply function

plot_create_wordcloud(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))

# plot

setwd(plots_dir)

#create function for all brain wordcloud comparison
# plot_wordcloud_brain<-function(dataset1, dataset1_name, dataset2, dataset2_name, dataset3, dataset3_name, dataset4, dataset4_name, words_to_remove){
#   # concatentate data
#   all_temp<-c(paste(dataset1, collapse =" "),
#               paste(dataset2, collapse =" "),
#               paste(dataset3, collapse =" "),
#               paste(dataset4, collapse =" "))
#   #convert to corpus
#   all<-Corpus(VectorSource(all_temp))
#   #convert everything to lowercase
#   all <- tm_map(all, content_transformer(tolower))
#   # stem words - i.e learning > learn, walked > walk
#   #all <- tm_map(all, stemDocument)
#   # remove unwanted grammer
#   all <- tm_map(all, removePunctuation)
#   all <- tm_map(all, removeWords, stopwords('english'))
#   # remove unwanted words
#   all <- tm_map(all, removeWords, words_to_remove)
#   #convert to matrix
#   all_matrix<-as.matrix(TermDocumentMatrix(all))
#   # rename columns
#   colnames(all_matrix)<-c(dataset1_name, dataset2_name, dataset3_name, dataset4_name)
#   comparison.cloud(all_matrix, max.words=1000,random.order=FALSE,colors=c("#1F497D","#C0504D", "#458B00", "#8B008B"),rot.per=.15, title.size = 2)
# }
# 
# plot_wordcloud_brain(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Frontal Lobe",
#                      Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Parietal Lobe",
#                      Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[,3],
#                      "Cerebellum",
#                      Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Temporal Lobe",
#                      c("sapiens", "homo", "human", "pathway", "action"))

# setwd(plots_dir)
# 
# pdf("Analysis_1_wordcloud_all_brain_regions_comparison.pdf", height=10, width =10)
# plot_wordcloud_brain(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Frontal Lobe",
#                      Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Parietal Lobe",
#                      Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[,3],
#                      "Cerebellum",
#                      Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3],
#                      "Temporal Lobe",
#                      c("sapiens", "homo", "human", "pathway", "action"))
# dev.off()


pdf("Analysis_1_comparison_of_pathways_across_brain_regions.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "pathway", "action", "system", "human"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "pathway", "action", "system", "human"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "pathway", "action", "system", "human"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "pathway", "action", "system", "human"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

##### ANALYSIS 1.3: - WORD CLOUD - superclass - q0.01 ######

# pathways mapped manually to parent pathways

setwd(Result_Tables_dir)

list.files()

Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")
Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep=",")

head(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass)

plot_create_wordcloud(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))

#word cloud changes font size of pathways in one brain region if other brain regions has sig more mention of the word.
#want to keep each brain word clud independent as cerebellum overpowers remaining brain regions.
# plot_wordcloud_brain(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2],
#                      "Frontal Lobe",
#                      Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2],
#                      "Parietal Lobe",
#                      Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass[,2],
#                      "Cerebellum",
#                      Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2],
#                      "Temporal Lobe",
#                      c("sapiens", "homo", "human", "pathway", "action", "system"))

layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))


setwd(plots_dir)

pdf("Analysis_1_comparison_of_superclass_pathways_across_brain_regions.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

##### ANALYSIS 1.4: KNOWN AD ASSOCIATED PATHWAYS #####

#grepping th efollowing pathways from results to see if they exist:

# immune
# tranlsation/transcription
# calcium
# RA
# synapse
# metabloism
# neurotransmitters
# wnt

Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_1_AD_Cerebellum_all_DEG_pathways_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_1_AD_Parietal_Lobe_all_DEG_pathways_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_1_AD_Frontal_Lobe_all_DEG_pathways_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]

##### ANLYSIS 1.5: BRAIN REGION UNIQUE PATHWAYS #####



##### ANALYSIS 2.1: MOST SIGNIFCANT PATHWAYS PER BRAIN REGION - AD SPECIFIC GENE LIST #####

# Cerebellum

head(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only)[3]
dim(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only)

# Parietal_Lobe
head(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only)[3]
dim(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only)

# Frontal_Lobe
head(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only)[3]
dim(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only)

##### ANALYSIS 2.2: COMMON PATHWAYS ACROSS BRAIN REGION -  AD SPECIFIC GENE LIST #####

# merge results
Analysis_2_AD_specific_DEG_across_brain_regions<-as.data.frame(table(rbind(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only[3],
                                                                           Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                           Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[3],
                                                                           Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[3])))

#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions<-Analysis_2_AD_specific_DEG_across_brain_regions[order(-Analysis_2_AD_specific_DEG_across_brain_regions$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions)

head(Analysis_2_AD_specific_DEG_across_brain_regions)

subset(Analysis_2_AD_specific_DEG_across_brain_regions, Freq==4)

## repeat with only 3 brain regions - removing cerebellum

# merge results
Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum<-as.data.frame(table(rbind(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[3],
                                                                                         Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[3],
                                                                                         Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[3])))

#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum<-Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum[order(-Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum)

head(Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_no_cerebellum, Freq==3)

##### ANALYSIS 2.2: EXPORT TABLE OF PATHWAYS ACROSS BRAIN REGIONS WITH Q-value q<=0.01 #####

#apply function to subset dataset to certain q Value and print results with specific columns
Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only, 0.01)
Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only, 0.01)
Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only, 0.01)

setwd(Result_Tables_dir)

write.table(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table, "Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")
write.table(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table, "Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")
write.table(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table, "Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")

setwd(work_dir)

##### ANALYSIS 2.2: - HEATMAP - q=0.01 ####

# subset datasets to 0.01 sig

Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_001<-subset(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only, q.value<=0.01)
Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_001<-subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only, q.value<=0.01)
Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_001<-subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only, q.value<=0.01)



Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_001, "Cerebellum")
Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_heatmap
Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_001, "parietal_Lobe")
Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_heatmap<-rearrange_pathway_data_for_heatmap(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_001, "Frontal_Lobe")

# merge results

Analysis_2_AD_specific_DEG_across_brain_regions_heatmap<-as.matrix(Reduce(MyMerge, list(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_heatmap,
                                                                                        Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_heatmap,
                                                                                        Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_heatmap,
                                                                                        Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_heatmap)))


head(Analysis_2_AD_specific_DEG_across_brain_regions_heatmap)
dim(Analysis_2_AD_specific_DEG_across_brain_regions_heatmap)

# plot
setwd(plots_dir)

heatmap.2(Analysis_2_AD_specific_DEG_across_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_2_AD_specific_DEG_across_all_brain_regions_pathway_heatmap.pdf", height=25, width=15)
heatmap.2(Analysis_2_AD_specific_DEG_across_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()

##### ANALYSIS 2.2: - WORD CLOUD  - q0.05 #####

# create function for word cloud plot

# apply function

plot_create_wordcloud(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))

# plot

setwd(plots_dir)

pdf("Analysis_2_wordcloud_Cerebellum.pdf")
plot_create_wordcloud(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
dev.off()

pdf("Analysis_2_wordcloud_Parietal_Lobe.pdf")
plot_create_wordcloud(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
dev.off()

pdf("Analysis_2_wordcloud_Frontal_Lobe.pdf")
plot_create_wordcloud(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
dev.off()

# brain plot
plot_wordcloud_brain(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[,3],
                     "Frontal Lobe",
                     Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[,3],
                     "Parietal Lobe",
                     Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only[,3],
                     "Cerebellum",
                     Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3],
                     "Temporal Lobe",
                     c("sapiens", "homo", "human", "pathway", "action"))

setwd(plots_dir)

pdf("Analysis_2_wordcloud_all_brain_regions_comparison.pdf", height=10, width =10)
plot_wordcloud_brain(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only[,3],
                     "Frontal Lobe",
                     Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only[,3],
                     "Parietal Lobe",
                     Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only[,3],
                     "Cerebellum",
                     Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only[,3],
                     "Temporal Lobe",
                     c("sapiens", "homo", "human", "pathway", "action", "action"))
dev.off()

##### ANALYSIS 2.2: - WORD CLOUD - superclass - q0.01 ####

# pathways mapped manually to parent pathways

setwd(Result_Tables_dir)

list.files()

Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")
Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")
Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")

head(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass)
head(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass)

plot_create_wordcloud(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system"))


layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))


setwd(plots_dir)

pdf("Analysis_2_comparison_of_superclass_pathways_across_brain_regions.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(Analysis_1_AD_Temporal_Lobe_all_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

##### ANALYSIS 2.3: KNOWN AD ASSOCIATED PATHWAYS #####

#grepping th efollowing pathways from results to see if they exist:

# immune
# tranlsation/transcription
# calcium
# RA
# synapse
# metabloism
# neurotransmitters
# wnt

Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_2_AD_Cerebellum_unique_DEG_pathways_results_sig_only_table[,1],ignore.case=TRUE),]
Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_2_AD_Parietal_Lobe_unique_DEG_pathways_results_sig_only_table[,1],ignore.case=TRUE),]
Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table[grep("translation|transcription|immune|calcium|RA|rheumatoid|chemical|synapse|neuro|wnt|meta", Analysis_2_AD_Frontal_Lobe_unique_DEG_pathways_results_sig_only_table[,1],ignore.case=TRUE),]

##### ANALYSIS 3.1: COMMPON PATHWAYS - BRAIN REGION - COMMON DEG LIST #####

# Cerebellum

head(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only)[3]
dim(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only)

# Parietal_Lobe
head(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only)[3]
dim(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only)

# Frontal_Lobe
head(Analysis_3_AD_Frontal_Lobe_common_DEG_pathways_results_sig_only)[3]
dim(Analysis_3_AD_Frontal_Lobe_common_DEG_pathways_results_sig_only)

##### ANALYSIS 3.1: EXPORT TABLE OF PATHWAYS ACROSS BRAIN REGIONS WITH Q-value q<=0.01 #####

#apply function to subset dataset to certain q Value and print results with specific columns
Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table<-subset_dataset(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only, 0.01)

setwd(Result_Tables_dir)

write.table(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table, "Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table.txt", quote=F, row.names = F, sep="\t")

setwd(work_dir)

##### ANALYSIS 3.1. - HEATMAP q0.05 #####

# subset datasets to 0.01 sig - leaves no pathways

#format

Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp<-Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only[c(2,3)]
rownames(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp)<-Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp$pathway
Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp$pathway<-NULL
colnames(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp)<-"Cerebellum"

Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp

Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp<-Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only[c(2,3)]
rownames(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp)<-Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp$pathway
Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp$pathway<-NULL
colnames(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp)<-"Parietal_Lobe"

Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp

# merge results

Analysis_3_common_DEG_brain_regions_heatmap<-as.matrix(Reduce(MyMerge, list(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only_temp,
                                                                            Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_temp)))

setwd(plots_dir)

heatmap.2(Analysis_3_common_DEG_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="none",
          Rowv = FALSE, 
          Colv=F,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_3_common_DEG_across_all_brain_regions_pathway_heatmap.pdf", height=25, width=15)
heatmap.2(Analysis_3_common_DEG_brain_regions_heatmap,
          margins = c(12,50),
          dendrogram="none",
          Rowv = FALSE, 
          Colv=F,
          cexCol=1.3,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()

##### ANALYSIS 3.1: - WORD CLOUD  - q0.05 #####

# create function for word cloud plot

# apply function

plot_create_wordcloud(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
plot_create_wordcloud(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))


# plot

setwd(plots_dir)

pdf("Analysis_3_wordcloud_Cerebellum.pdf")
plot_create_wordcloud(Analysis_3_AD_Cerebellum_common_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
dev.off()

pdf("Analysis_3_wordcloud_Parietal_Lobe.pdf")
plot_create_wordcloud(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only[,3], c("sapiens", "homo", "human", "pathway", "action"))
dev.off()

##### ANALYSIS 3.2: WORD CLOUD - superclass - q 0.01 #####

#pathways mapped manually to parent pathways

setwd(Result_Tables_dir)

list.files()

Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table_with_superclass<-read.table("Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table_with_superclass.csv", head=T, as.is=T, sep="\t")

head(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table_with_superclass)

plot_create_wordcloud(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))

setwd(plots_dir)

pdf("Analysis_3_Parietal_Lobe_common_superclass_pathways.pdf")
plot_create_wordcloud(Analysis_3_AD_Parietal_Lobe_common_DEG_pathways_results_sig_only_table_with_superclass[,2], c("sapiens", "homo", "pathway", "action", "system", "human", "organismal"))
dev.off()

##### SAVE IMAGE #####

setwd(work_dir)

save.image("CONSENSUSPATHDB_pathway_RESULTS_COMPILATION.Rdata")
