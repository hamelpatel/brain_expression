##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                            CREATE CMAP INPUT FILES AND QUERY IN ENRICHR                                                #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################


# DESCRIPTION
# READ IN MERGED de RESULTS
# CONVERT ENTREZ ID TO GENE SYMBOL
# EXTRACT COMMON BRAIN REGIONS BY EXPRESSION CHIP
# IDENTIFY COMMON DE BRAIN GENES USING DIFFERENT THRESHOLDS
# QUERY IN CMAP
# 


##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/1TB/Projects/Brain_expression/3.DE_analysis"

work_dir="/media/hamel/1TB/Projects/Brain_expression/4.Create_cMAP_input/"

setwd(work_dir)

dir.create(paste(work_dir, "cMAP_input_files", sep="/"))

cmap_input_dir=paste(work_dir, "cMAP_input_files", sep="/")

setwd(cmap_input_dir)

dir.create(paste(cmap_input_dir, "up_regulated_genes", sep="/"))

dir.create(paste(cmap_input_dir, "down_regulated_genes", sep="/"))

up_genes_dir=paste(cmap_input_dir, "up_regulated_genes", sep="/")

down_genes_dir=paste(cmap_input_dir, "down_regulated_genes", sep="/")

##### LOAD LIBRARIES ####

library(org.Hs.eg.db)

##### READ IN DE RESULTS #####

setwd(data_dir)

load("Create_common_features_and_DE.Rdata")

##### CREATE ENTREZ GENE ID TO GENE SYMBOL LOOK UP TABLE#####

# can't convert results to gene symbol as there are 632 dupliacted gene symbols when converting all entrez gene IDs

# convert entrez gene id to gene symbol
entrez_gene_symbol <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(entrez_gene_symbol)
# Convert to a list
gene_symbol_lookup_table <- as.data.frame(entrez_gene_symbol[mapped_genes])

head(gene_symbol_lookup_table)

##### CREATE FUNCTION TO EXTRACT SIG DOWN-REGULATED DE GENES #####

create_cmap_down_gene_list<-function(DE_table, adj.p.val_threshold, logFC_threshold){
  # subset DE result using user threshold threshold
  filtered<-DE_table[DE_table[7]<=adj.p.val_threshold & DE_table[1]<logFC_threshold,]
  # convert list of entrez id to tgene symbols
  gene_symbol<-gene_symbol_lookup_table[gene_symbol_lookup_table$gene_id %in% rownames(filtered),]
  # print to user number of sig genes
  print(c("number of genes", length(unique(gene_symbol$symbol))))
  # return list of genes
  return(gene_symbol$symbol)
}

# apply function to all datasets - each brain region separate
E_GEOD_28146_top_genes_Hippocampus_CA1_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_28146_top_genes_Hippocampus_CA1, 0.05, 0)
E_GEOD_29378_top_genes_CA1_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_29378_top_genes_CA1, 0.05, 0)
E_GEOD_29378_top_genes_CA3_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_29378_top_genes_CA3, 0.05, 0)
E_GEOD_36980_top_genes_Frontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_36980_top_genes_Frontal_Cortex, 0.05, 0)
E_GEOD_36980_top_genes_Hippocampus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_36980_top_genes_Hippocampus, 0.05, 0)
E_GEOD_36980_top_genes_Temporal_Cortex_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_36980_top_genes_Temporal_Cortex, 0.05, 0)
E_GEOD_48350_top_genes_Entorhinal_Cortex_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_48350_top_genes_Entorhinal_Cortex, 0.05, 0)
E_GEOD_48350_top_genes_Hippocampus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_48350_top_genes_Hippocampus, 0.05, 0)
E_GEOD_48350_top_genes_Postcentral_Gyrus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_48350_top_genes_Postcentral_Gyrus, 0.05, 0)
E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus, 0.05, 0)
E_GEOD_1297_top_genes_Hippocampus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_1297_top_genes_Hippocampus, 0.05, 0)
E_GEOD_5281_top_genes_Entorhinal_Cortex_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Entorhinal_Cortex, 0.05, 0)
E_GEOD_5281_top_genes_Hippocampus_CA1_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Hippocampus_CA1, 0.05, 0)
E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus, 0.05, 0)
E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus, 0.05, 0)
E_GEOD_5281_top_genes_Posterior_Singulate_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Posterior_Singulate, 0.05, 0)
E_GEOD_5281_top_genes_Primary_Visual_Cortex_sig_genes_down<-create_cmap_down_gene_list(E_GEOD_5281_top_genes_Primary_Visual_Cortex, 0.05, 0)
inhouse_data_top_genes_Entorhinal_Cortex_sig_genes_down<-create_cmap_down_gene_list(inhouse_data_clean_top_genes_Entorhinal_Cortex, 0.05, 0)
inhouse_data_top_genes_Cerebellum_sig_genes_down<-create_cmap_down_gene_list(inhouse_data_clean_top_genes_Cerebellum, 0.05, 0)
inhouse_data_top_genes_Frontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(inhouse_data_clean_top_genes_Frontal_Cortex, 0.05, 0)
inhouse_data_top_genes_Temporal_Cortex_sig_genes_down<-create_cmap_down_gene_list(inhouse_data_clean_top_genes_Temporal_Cortex, 0.05, 0)
AMP_MAYOeGWAS_top_genes_Temporal_Cortex_sig_genes_down<-create_cmap_down_gene_list(AMP_MAYOeGWAS_top_genes_Temporal_Cortex, 0.05, 0)
AMP_MAYOeGWAS_top_genes_Cerebellum_sig_genes_down<-create_cmap_down_gene_list(AMP_MAYOeGWAS_top_genes_Cerebellum, 0.05, 0)
AMP_MSBB_U133A_top_genes_Frontal_Pole_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Frontal_Pole, 0.05, 0)
AMP_MSBB_U133A_top_genes_Precentral_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Precentral_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule, 0.05, 0)
AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Hippocampus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Hippocampus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Temporal_Pole_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133A_top_genes_Temporal_Pole, 0.05, 0)
AMP_MSBB_U133B_top_genes_Frontal_Pole_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Frontal_Pole, 0.05, 0)
AMP_MSBB_U133B_top_genes_Precentral_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Precentral_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule, 0.05, 0)
AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Hippocampus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Hippocampus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Temporal_Pole_sig_genes_down<-create_cmap_down_gene_list(AMP_MSBB_U133B_top_genes_Temporal_Pole, 0.05, 0)

# write out to cmap input

setwd(down_genes_dir)

write.table(as.data.frame(unique(E_GEOD_28146_top_genes_Hippocampus_CA1_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_28146_top_genes_Hippocampus_CA1_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_29378_top_genes_CA1_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_29378_top_genes_CA1_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_29378_top_genes_CA3_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_29378_top_genes_CA3_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Frontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Frontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Hippocampus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Hippocampus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Temporal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Temporal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Entorhinal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Hippocampus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Hippocampus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Postcentral_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Postcentral_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_1297_top_genes_Hippocampus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_1297_top_genes_Hippocampus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Entorhinal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Hippocampus_CA1_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Hippocampus_CA1_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Posterior_Singulate_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Posterior_Singulate_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Primary_Visual_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Primary_Visual_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Entorhinal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Cerebellum_sig_genes_down)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Cerebellum_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Frontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Frontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Temporal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Temporal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MAYOeGWAS_top_genes_Temporal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MAYOeGWAS_top_genes_Temporal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MAYOeGWAS_top_genes_Cerebellum_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MAYOeGWAS_top_genes_Cerebellum_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Frontal_Pole_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Frontal_Pole_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Precentral_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Precentral_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Hippocampus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Hippocampus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Temporal_Pole_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Temporal_Pole_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Frontal_Pole_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Frontal_Pole_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Precentral_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Precentral_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Hippocampus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Hippocampus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_p_0.05_logFC_0_down.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Temporal_Pole_sig_genes_down)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Temporal_Pole_p_0.05_logFC_0_down.txt")

##### CREATE FUNCTION TO EXTRACT SIG UPREGULATED DE GENES #####

create_cmap_up_gene_list<-function(DE_table, adj.p.val_threshold, logFC_threshold){
  # subset DE result using user threshold threshold
  filtered<-DE_table[DE_table[7]<=adj.p.val_threshold & DE_table[1]>logFC_threshold,]
  # convert list of entrez id to tgene symbols
  gene_symbol<-gene_symbol_lookup_table[gene_symbol_lookup_table$gene_id %in% rownames(filtered),]
  # print to user number of sig genes
  print(c("number of genes", length(unique(gene_symbol$symbol))))
  # return list of genes
  return(gene_symbol$symbol)
}

# apply function to all datasets - each brain region separate
E_GEOD_28146_top_genes_Hippocampus_CA1_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_28146_top_genes_Hippocampus_CA1, 0.05, 0)
E_GEOD_29378_top_genes_CA1_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_29378_top_genes_CA1, 0.05, 0)
E_GEOD_29378_top_genes_CA3_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_29378_top_genes_CA3, 0.05, 0)
E_GEOD_36980_top_genes_Frontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_36980_top_genes_Frontal_Cortex, 0.05, 0)
E_GEOD_36980_top_genes_Hippocampus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_36980_top_genes_Hippocampus, 0.05, 0)
E_GEOD_36980_top_genes_Temporal_Cortex_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_36980_top_genes_Temporal_Cortex, 0.05, 0)
E_GEOD_48350_top_genes_Entorhinal_Cortex_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_48350_top_genes_Entorhinal_Cortex, 0.05, 0)
E_GEOD_48350_top_genes_Hippocampus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_48350_top_genes_Hippocampus, 0.05, 0)
E_GEOD_48350_top_genes_Postcentral_Gyrus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_48350_top_genes_Postcentral_Gyrus, 0.05, 0)
E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus, 0.05, 0)
E_GEOD_1297_top_genes_Hippocampus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_1297_top_genes_Hippocampus, 0.05, 0)
E_GEOD_5281_top_genes_Entorhinal_Cortex_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Entorhinal_Cortex, 0.05, 0)
E_GEOD_5281_top_genes_Hippocampus_CA1_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Hippocampus_CA1, 0.05, 0)
E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus, 0.05, 0)
E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus, 0.05, 0)
E_GEOD_5281_top_genes_Posterior_Singulate_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Posterior_Singulate, 0.05, 0)
E_GEOD_5281_top_genes_Primary_Visual_Cortex_sig_genes_up<-create_cmap_up_gene_list(E_GEOD_5281_top_genes_Primary_Visual_Cortex, 0.05, 0)
inhouse_data_top_genes_Entorhinal_Cortex_sig_genes_up<-create_cmap_up_gene_list(inhouse_data_clean_top_genes_Entorhinal_Cortex, 0.05, 0)
inhouse_data_top_genes_Cerebellum_sig_genes_up<-create_cmap_up_gene_list(inhouse_data_clean_top_genes_Cerebellum, 0.05, 0)
inhouse_data_top_genes_Frontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(inhouse_data_clean_top_genes_Frontal_Cortex, 0.05, 0)
inhouse_data_top_genes_Temporal_Cortex_sig_genes_up<-create_cmap_up_gene_list(inhouse_data_clean_top_genes_Temporal_Cortex, 0.05, 0)
AMP_MAYOeGWAS_top_genes_Temporal_Cortex_sig_genes_up<-create_cmap_up_gene_list(AMP_MAYOeGWAS_top_genes_Temporal_Cortex, 0.05, 0)
AMP_MAYOeGWAS_top_genes_Cerebellum_sig_genes_up<-create_cmap_up_gene_list(AMP_MAYOeGWAS_top_genes_Cerebellum, 0.05, 0)
AMP_MSBB_U133A_top_genes_Frontal_Pole_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Frontal_Pole, 0.05, 0)
AMP_MSBB_U133A_top_genes_Precentral_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Precentral_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule, 0.05, 0)
AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Hippocampus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Hippocampus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133A_top_genes_Temporal_Pole_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133A_top_genes_Temporal_Pole, 0.05, 0)
AMP_MSBB_U133B_top_genes_Frontal_Pole_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Frontal_Pole, 0.05, 0)
AMP_MSBB_U133B_top_genes_Precentral_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Precentral_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule, 0.05, 0)
AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex, 0.05, 0)
AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Hippocampus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Hippocampus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus, 0.05, 0)
AMP_MSBB_U133B_top_genes_Temporal_Pole_sig_genes_up<-create_cmap_up_gene_list(AMP_MSBB_U133B_top_genes_Temporal_Pole, 0.05, 0)

# write out to cmap input

setwd(up_genes_dir)

write.table(as.data.frame(unique(E_GEOD_28146_top_genes_Hippocampus_CA1_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_28146_top_genes_Hippocampus_CA1_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_29378_top_genes_CA1_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_29378_top_genes_CA1_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_29378_top_genes_CA3_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_29378_top_genes_CA3_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Frontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Frontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Hippocampus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Hippocampus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_36980_top_genes_Temporal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_36980_top_genes_Temporal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Entorhinal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Hippocampus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Hippocampus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Postcentral_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Postcentral_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_48350_top_genes_Superior_Frontal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_1297_top_genes_Hippocampus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_1297_top_genes_Hippocampus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Entorhinal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Hippocampus_CA1_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Hippocampus_CA1_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Medial_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Superior_Frontal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Posterior_Singulate_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Posterior_Singulate_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(E_GEOD_5281_top_genes_Primary_Visual_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="E_GEOD_5281_top_genes_Primary_Visual_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Entorhinal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Entorhinal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Cerebellum_sig_genes_up)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Cerebellum_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Frontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Frontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(inhouse_data_top_genes_Temporal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="inhouse_data_top_genes_Temporal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MAYOeGWAS_top_genes_Temporal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MAYOeGWAS_top_genes_Temporal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MAYOeGWAS_top_genes_Cerebellum_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MAYOeGWAS_top_genes_Cerebellum_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Frontal_Pole_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Frontal_Pole_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Precentral_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Precentral_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Inferior_Frontal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Dorsolateral_Prefrontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Superior_Parietal_Lobule_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Prefrontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Parahippocampal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Hippocampus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Hippocampus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Inferior_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Middle_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Superior_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133A_top_genes_Temporal_Pole_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133A_top_genes_Temporal_Pole_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Frontal_Pole_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Frontal_Pole_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Precentral_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Precentral_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Inferior_Frontal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Dorsolateral_Prefrontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Superior_Parietal_Lobule_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Prefrontal_Cortex_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Parahippocampal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Hippocampus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Hippocampus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Inferior_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Middle_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Superior_Temporal_Gyrus_p_0.05_logFC_0_up.txt")
write.table(as.data.frame(unique(AMP_MSBB_U133B_top_genes_Temporal_Pole_sig_genes_up)), row.names=F, col.names=F, quote=F, file="AMP_MSBB_U133B_top_genes_Temporal_Pole_p_0.05_logFC_0_up.txt")

##### MANUAL CHECK IN enrichR cmap and DOWNLOAD FULL old CMAP DOWN table #####

##### CHECK DAPH DRIVERS in ALL DATASETS #####

# read cmap files - downloaded from enrichR

cmap_down<-read.table("/media/hamel/1TB/Projects/Brain_expression/5.cMAP_result_analysis/Old_CMAP_down", fill=T)
cmap_up<-read.table("/media/hamel/1TB/Projects/Brain_expression/5.cMAP_result_analysis/Old_CMAP_up", fill=T)

head(cmap_down)[1:5]
head(cmap_up)[1:5]

#subset to "4,5-dianilinophthalimide -624"
subset(cmap_down, V1=="4,5-dianilinophthalimide" & V2=="-624")
subset(cmap_up, V1=="4,5-dianilinophthalimide" & V2=="-624")

DAPH_down_regulated_genes<-as.data.frame(t(subset(cmap_up, V1=="4,5-dianilinophthalimide" & V2=="-624")))
DAPH_up_regulated_genes<-as.data.frame(t(subset(cmap_down, V1=="4,5-dianilinophthalimide" & V2=="-624")))

head(DAPH_down_regulated_genes)
head(DAPH_up_regulated_genes)

#remove 1st two rows
dim(DAPH_down_regulated_genes)
dim(DAPH_up_regulated_genes)

DAPH_down_regulated_genes<-as.data.frame(DAPH_down_regulated_genes[3:dim(DAPH_down_regulated_genes)[1],])
DAPH_up_regulated_genes<-as.data.frame(DAPH_up_regulated_genes[3:dim(DAPH_up_regulated_genes)[1],])

dim(DAPH_down_regulated_genes)
dim(DAPH_up_regulated_genes)

colnames(DAPH_down_regulated_genes)<-"DAPH_down"
colnames(DAPH_up_regulated_genes)<-"DAPH_up"

head(DAPH_down_regulated_genes)
head(DAPH_up_regulated_genes)

# convert list to entrez gene ID

head(gene_symbol_lookup_table)

# convert list of entrez id to tgene symbols
DAPH_down_regulated_genes_entrez<-gene_symbol_lookup_table[gene_symbol_lookup_table$symbol %in% DAPH_down_regulated_genes$DAPH_down ,]
head(DAPH_down_regulated_genes_entrez)
dim(DAPH_down_regulated_genes_entrez)

DAPH_up_regulated_genes_entrez<-gene_symbol_lookup_table[gene_symbol_lookup_table$symbol %in% DAPH_up_regulated_genes$DAPH_up ,]
head(DAPH_up_regulated_genes_entrez)
dim(DAPH_up_regulated_genes_entrez)

# extract logFC values
LogFC_only<-merged_diff_exprs[grep("logFC", colnames(merged_diff_exprs))]
head(LogFC_only)[1:10]

# keep only DAPH entrez
LogFC_DAPH_up<-LogFC_only[rownames(LogFC_only) %in% DAPH_up_regulated_genes_entrez$gene_id,]
head(LogFC_DAPH_up)

rowSums(LogFC_DAPH_up)
boxplot(LogFC_DAPH_up, las=2, par(mar=c(14,3,3,2)), cex.axis=0.5)
abline(h=0)

# keep only DAPH entrez
LogFC_DAPH_down<-LogFC_only[rownames(LogFC_only) %in% DAPH_down_regulated_genes_entrez$gene_id,]
head(LogFC_DAPH_down)

rowSums(LogFC_DAPH_down)
boxplot(LogFC_DAPH_down, las=2, par(mar=c(14,3,3,2)), cex.axis=0.5)
abline(h=0)

##### SAVE IMAGE #####

setwd(work_dir)

save.image("create_cmap_input.Rdata")
