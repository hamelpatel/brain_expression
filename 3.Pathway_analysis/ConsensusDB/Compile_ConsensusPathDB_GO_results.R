########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                         CONSENSUSPATHDB GENE ONTOLOGY RESULTS COMPILATION                            #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 16/05/2017

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

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/8.Pathway_Analysis/8.5.ConsensusPathDB/8.5.2.GO_onotologies"

dir.create(paste(work_dir, "Plots", sep="/"))
plots_dir=paste(work_dir, "Plots", sep="/")

dir.create(paste(work_dir, "Result_Tables", sep="/"))
Result_Tables_dir=paste(work_dir, "Result_Tables", sep="/")

##### READ DATA  - ONLY PATHWAY ######

setwd(work_dir)

# files in directory
file_names<-list.files()

# keep only .txt files
file_names<-file_names[grep("ontology_results.tab", file_names)]

for (i in 1:length(file_names)) 
  assign(file_names[i], 
         read.table(file_names[i], head=T, sep="\t", quote="", fill=T))

ls()

#general check
head(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results.tab)[1:6]

##### SUBSET TO SIG ONLY #####

# subset to significant data to where q value<=0.05 - (FDR adjusted p value)

for (x in 1:length(file_names)) {
  # subset to FDR sig only and assign a new object name
  assign(paste(gsub(".tab", "_sig_only", file_names[x])), subset(get(file_names[x]), q.value<=0.05))
}

dim(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results.tab)

dim(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only)

#dim each sig file
sig_file_names<-gsub(".tab", "_sig_only", file_names)

for (x in 1:length(sig_file_names)) {
  print(sig_file_names[x])
  print(dim(get(sig_file_names[x])))
}

##### ANALYSIS 1.1: MICROARRAY + RNA_SEQ OVERLAP #####

#Temporal Lobe - microarray
Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only[6]

#Temporal Lobe - RNA-Seq
Analysis_1_AD_RNA_SEQ_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only[6]

# any microarray pathway in RNA-seq data
subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only$pathway %in% Analysis_1_AD_RNA_SEQ_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only$pathway) 

#Frontal Lobe - microarray
Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only[6]

#Frontal Lobe - RNA-Seq
Analysis_1_AD_RNA_SEQ_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only[6]

# any microarray pathway in RNA-seq data
subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only$pathway %in% Analysis_1_AD_RNA_SEQ_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only$pathway) 

##### ANALYSIS 1.2: MOST SIGNIFCANT gene_ontology PER BRAIN REGION  - FULL DEG LIST #####

# Cerebellum

head(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Temporal_Lobe
head(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Parietal_Lobe
head(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Frontal_Lobe
head(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

##### ANALYSIS 1.3: COMMON gene_ontology ACROSS BRAIN REGION q=0.01 - FULL DEG LIST #####

# merge results (CELLULAR) 

Analysis_1_all_DEG_across_brain_regions_c<-as.data.frame(table(rbind(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions_c<-Analysis_1_all_DEG_across_brain_regions_c[order(-Analysis_1_all_DEG_across_brain_regions_c$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_c)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_c)

head(Analysis_1_all_DEG_across_brain_regions_c)

subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==4)

# merge results (MOLECULAR)
Analysis_1_all_DEG_across_brain_regions_m<-as.data.frame(table(rbind(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6])))


#arrange by count
Analysis_1_all_DEG_across_brain_regions_m<-Analysis_1_all_DEG_across_brain_regions_m[order(-Analysis_1_all_DEG_across_brain_regions_m$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_m)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_m)

head(Analysis_1_all_DEG_across_brain_regions_m)

subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==4)


# merge results (BIOLOGICAL)
Analysis_1_all_DEG_across_brain_regions_b<-as.data.frame(table(rbind(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                     subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6])))


#arrange by count
Analysis_1_all_DEG_across_brain_regions_b<-Analysis_1_all_DEG_across_brain_regions_b[order(-Analysis_1_all_DEG_across_brain_regions_b$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_b)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_b)

head(Analysis_1_all_DEG_across_brain_regions_b)

subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==4)

##### ANALYSIS 1.3: - HEATMAP - q=0.01 #####

# using subset of results where q<=0.01

# create frunction to rearrange data for
rearrange_pathway_data_for_heatmap<-function(data, column_name, component){
  #subset data to component
  data<-subset(data, term_category==component & q.value<=0.01)
  # move term name to rownames
  rownames(data)<-data$term_name
  #keep only q values
  data<-data[2]
  colnames(data)<-column_name
  return(data)
}

Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_c<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, "Cerebellum", "c")
Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, "Temporal_Lobe", "c")
Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, "parietal_Lobe", "c")
Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, "Frontal_Lobe", "c")

# merge results

MyMerge <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= TRUE, all.y=T)
  rownames(df)<-df$Row.names
  df$Row.names<-NULL
  return(df)
}

Analysis_1_all_DEG_across_brain_regions_heatmap_c<-as.matrix(Reduce(MyMerge, list(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_c,
                                                                                  Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c,
                                                                                  Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c,
                                                                                  Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_c)))


head(Analysis_1_all_DEG_across_brain_regions_heatmap_c)
dim(Analysis_1_all_DEG_across_brain_regions_heatmap_c)

# plot
setwd(plots_dir)

heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_c,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_1_all_DEG_across_all_brain_regions_GO_CELLULAR_heatmap.pdf", height=25, width=15)
heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_c,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()

##
## MOLECULAR 
##

Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_m<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, "Cerebellum", "m")
Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, "Temporal_Lobe", "m")
Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, "parietal_Lobe", "m")
Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, "Frontal_Lobe", "m")

# merge results

Analysis_1_all_DEG_across_brain_regions_heatmap_m<-as.matrix(Reduce(MyMerge, list(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_m,
                                                                                  Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m,
                                                                                  Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m,
                                                                                  Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_m)))


head(Analysis_1_all_DEG_across_brain_regions_heatmap_m)
dim(Analysis_1_all_DEG_across_brain_regions_heatmap_m)

# plot
setwd(plots_dir)

heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_m,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_1_all_DEG_across_all_brain_regions_GO_MOLECULAR_heatmap.pdf", height=25, width=15)
heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_m,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()


##
## BIOLOGICAL
##

Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_b<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, "Cerebellum", "b")
Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, "Temporal_Lobe", "b")
Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, "parietal_Lobe", "b")
Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b<-rearrange_pathway_data_for_heatmap(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, "Frontal_Lobe", "b")

# merge results

Analysis_1_all_DEG_across_brain_regions_heatmap_b<-as.matrix(Reduce(MyMerge, list(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only_heatmap_b,
                                                                                  Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b,
                                                                                  Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b,
                                                                                  Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only_heatmap_b)))


head(Analysis_1_all_DEG_across_brain_regions_heatmap_b)
dim(Analysis_1_all_DEG_across_brain_regions_heatmap_b)

# plot
setwd(plots_dir)

heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_b,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")

pdf("Analysis_1_all_DEG_across_all_brain_regions_GO_BIOLOGICAL_heatmap.pdf", height=50, width=15)
heatmap.2(Analysis_1_all_DEG_across_brain_regions_heatmap_b,
          margins = c(12,50),
          dendrogram="column",
          Rowv = FALSE, 
          Colv=T,
          cexCol=1.6,
          cexRow=1,
          keysize = 1,
          trace="none")
dev.off()

##### ANALYSIS 1.3: - WORD CLOUD  - q0.01 #####

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

## CELLULAR

# apply function

plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))

# plot

setwd(plots_dir)

pdf("Analysis_1_comparison_of_CELLULAR_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

## MOLECULAR

# apply function

plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))

# plot

setwd(plots_dir)

pdf("Analysis_1_comparison_of_MOLECULAR_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

## BIOLOGICAL

# apply function

plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))

# plot

setwd(plots_dir)

pdf("Analysis_1_comparison_of_BIOLOGICAL_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

##### ANLYSIS 1.4: BRAIN REGION UNIQUE gene_ontology #####

subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)

# cerebellum - 

(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Cerebellum_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6]

subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)

# Temporal_Lobe - 

(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6]

dim((subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6])

# Parietal_Lobe - 

(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6]

dim((subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6])

subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)
subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)

# Frontal_Lobe - 

(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6]
(subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6]

dim((subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_c, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_m, Freq==1)[,1]))[6])
dim((subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_name %in%  subset(Analysis_1_all_DEG_across_brain_regions_b, Freq==1)[,1]))[6])

##### ANALYSIS 1.5: COMMON GO IN BRAIN WITHOUT CEREBELLUM #####

# merge results (CELLULAR) 

Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c<-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c[order(-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c)

head(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c)

subset(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_c, Freq==3)

# merge results (MOLECULAR)
Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6])))


#arrange by count
Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m<-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m[order(-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m)

head(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m)

subset(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_m, Freq==3)


# merge results (BIOLOGICAL)
Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Parietal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                                    subset(Analysis_1_AD_Frontal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6])))


#arrange by count
Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b<-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b[order(-Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b)

head(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b)

subset(Analysis_1_all_DEG_across_brain_regions_NO_CEREBBELLUM_b, Freq==3)

##### ANALYSIS 2.1: MOST SIGNIFCANT gene_ontology PER BRAIN REGION - AD SPECIFIC GENE LIST #####

# Cerebellum

head(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Parietal_Lobe
head(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Frontal_Lobe
head(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

##### ANALYSIS 2.2: COMMON gene_ontology ACROSS BRAIN REGION q=0.01 -  AD SPECIFIC GENE LIST #####

# merge results (CELLULAR) 

Analysis_2_AD_specific_DEG_across_brain_regions_c<-as.data.frame(table(rbind(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                             subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6])))

#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_c<-Analysis_2_AD_specific_DEG_across_brain_regions_c[order(-Analysis_2_AD_specific_DEG_across_brain_regions_c$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_c)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_c)

head(Analysis_2_AD_specific_DEG_across_brain_regions_c)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_c, Freq==4)

# merge results (MOLECULAR)
Analysis_2_AD_specific_DEG_across_brain_regions_m<-as.data.frame(table(rbind(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                             subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6])))


#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_m<-Analysis_2_AD_specific_DEG_across_brain_regions_m[order(-Analysis_2_AD_specific_DEG_across_brain_regions_m$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_m)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_m)

head(Analysis_2_AD_specific_DEG_across_brain_regions_m)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_m, Freq==4)


# merge results (BIOLOGICAL)
Analysis_2_AD_specific_DEG_across_brain_regions_b<-as.data.frame(table(rbind(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                             subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                             subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6])))


#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_b<-Analysis_2_AD_specific_DEG_across_brain_regions_b[order(-Analysis_2_AD_specific_DEG_across_brain_regions_b$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_b)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_b)

head(Analysis_2_AD_specific_DEG_across_brain_regions_b)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_b, Freq==4)

##### ANALYSIS 2.2: - WORD CLOUD  - q0.01 #####

## CELLULAR

# apply function

plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))

# plot

setwd(plots_dir)

pdf("Analysis_2_comparison_of_CELLULAR_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[,6], c("part", "complex"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

## MOLECULAR

# apply function

plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))

# plot

setwd(plots_dir)

pdf("Analysis_2_comparison_of_MOLECULAR_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[,6], c(""))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

## BIOLOGICAL

# apply function

plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))

# plot

setwd(plots_dir)

pdf("Analysis_2_comparison_of_BIOLOGICAL_gene_ontology_across_brain_regions_wordcloud.pdf")
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mfrow=c(2,2))
plot_create_wordcloud(subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Parietal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Frontal Lobe")
plot_create_wordcloud(subset(Analysis_2_AD_Cerebellum_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Cerebellum")
plot_create_wordcloud(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[,6], c("process", "regulation"))
text(x=0.5, y=1, cex=2, "Temporal Lobe")
par(mfrow=c(1,1))
dev.off()

##### ANALYSIS 2.3: COMMON GO IN BRAIN WITHOUT CEREBELLUM - AD SPECIFIC LIST #####

# merge results (CELLULAR) 

Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="c" & q.value<=0.01)[6])))

#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c<-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c[order(-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c)

head(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_c, Freq==3)

# merge results (MOLECULAR)
Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="m" & q.value<=0.01)[6])))


#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m<-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m[order(-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m)

head(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_m, Freq==3)


# merge results (BIOLOGICAL)
Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b<-as.data.frame(table(rbind(subset(Analysis_1_AD_Temporal_Lobe_all_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Parietal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6],
                                                                                            subset(Analysis_2_AD_Frontal_Lobe_unique_DEG_gene_ontology_results_sig_only, term_category=="b" & q.value<=0.01)[6])))


#arrange by count
Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b<-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b[order(-Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b$Freq),]
rownames(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b)<-1:nrow(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b)

head(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b)

subset(Analysis_2_AD_specific_DEG_across_brain_regions_NO_CEREBBELLUM_b, Freq==3)

##### ANALYSIS 3.1: COMMPON gene_ontology - BRAIN REGION - COMMON DEG LIST #####

# Cerebellum

head(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_3_AD_Cerebellum_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Parietal_Lobe
head(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_3_AD_Parietal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

# Frontal_Lobe
head(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
head(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
head(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])


dim(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="c")[6])
dim(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="m")[6])
dim(subset(Analysis_3_AD_Frontal_Lobe_common_DEG_gene_ontology_results_sig_only, term_category=="b")[6])

##### SAVE IMAGE #####

setwd(work_dir)

save.image("CONSENSUSPATHDB_GO_RESULTS_COMPILATION.Rdata")
