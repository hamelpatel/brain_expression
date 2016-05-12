##### LIBRARIES #####

library(sva)
library(limma)
library("org.Hs.eg.db")

##### LOAD ORIGINAL INHOUSE DATA - EXTRACT ENTORHINAL CORTEX GOOD PROBES #####

head(eset_bg_log2_rsn_combat_adj)[,1:5]

colnames(pData(eset_bg_log2_rsn_combat_adj))

#extract exprs results
exprs<-as.data.frame(exprs(eset_bg_log2_rsn_combat_adj))

head(exprs)[1:5]
dim(exprs)

head(good_probes)
length(good_probes)

#extract good probes
exprs_good_probes<-subset(exprs, rownames(exprs)%in%Entorhinal_Cortex_probes)
dim(exprs_good_probes)
head(exprs_good_probes)[1:5]

ls()

head(Entorhinal_Cortex_samples)
length(Entorhinal_Cortex_samples)

#extract EC samples
exprs_good_probes_EC<-exprs_good_probes[,colnames(exprs_good_probes)%in%Entorhinal_Cortex_samples]
dim(exprs_good_probes_EC)
head(exprs_good_probes_EC)[1:5]

#diagnosis
diagnosis<-as.data.frame(pData(eset_bg_log2_rsn_combat_adj)[25])
head(diagnosis)
table(diagnosis)

#subset diagnosis to EC samples and AD vs CO
dim(diagnosis)
head(diagnosis)
diagnosis_CO_AD<-subset(diagnosis, PHENOTYPE=="CO" | PHENOTYPE=="AD")
dim(diagnosis_CO_AD)
table(diagnosis_CO_AD)
diagnosis_CO_AD<-droplevels(diagnosis_CO_AD)
table(diagnosis_CO_AD)

# subset exprs table to case control

exprs_good_probes_EC_AD_CO<- exprs_good_probes_EC[,colnames(exprs_good_probes_EC) %in% rownames(diagnosis_CO_AD)]
dim(exprs_good_probes_EC_AD_CO)

# merge exprs + diagnosis

head(exprs_good_probes_EC_AD_CO)[1:5]
head(diagnosis_CO_AD)

exprs_good_probes_EC_AD_CO_diag<-merge(diagnosis_CO_AD, t(exprs_good_probes_EC_AD_CO), by="row.names")
dim(exprs_good_probes_EC_AD_CO_diag)
head(exprs_good_probes_EC_AD_CO_diag)[1:5]
colnames(exprs_good_probes_EC_AD_CO_diag)[2]<-"Diagnosis"
rownames(exprs_good_probes_EC_AD_CO_diag)<-exprs_good_probes_EC_AD_CO_diag$Row.names
exprs_good_probes_EC_AD_CO_diag$Row.names<-NULL

head(exprs_good_probes_EC_AD_CO_diag)[1:5]
dim(exprs_good_probes_EC_AD_CO_diag)
table(exprs_good_probes_EC_AD_CO_diag$Diagnosis)

##### RUN LIMMA ON ORIGINAL DATA #####

# limma
library("limma")

# create function to run DE analysis

run_diff_exprs_analysis<-function(dataset){
  #sort dataset by Diagnosis column - want AD samples 1st
  dataset_sorted<-dataset[order(dataset$Diagnosis, decreasing=T),]
  #split dataset into expres and diagnosis
  dataset_exprs<-dataset_sorted[2:dim(dataset_sorted)[2]]
  dataset_pheno<-dataset_sorted[,1]
  #replace sample name to numbers
  rownames(dataset_exprs)<-c(1:dim(dataset)[1])
  # setup experimental desgin
  design <- model.matrix(~0 + dataset_pheno)
  # change colnames to AD and Control
  colnames(design)<-c("CO", "AD")
  # transpose dataset, convert to numeric 
  transposed_dataset_exprs<-t(dataset_exprs)
  #run diff expression
  dataset_exprs_fit <- lmFit(transposed_dataset_exprs, design, method="robust")
  dataset_exprs_contrast_matrix<- makeContrasts(AD-CO, levels=design)
  dataset_exprs_contrast_fit <-contrasts.fit(dataset_exprs_fit, dataset_exprs_contrast_matrix)
  dataset_exprs_contrast_ebayes <- eBayes(dataset_exprs_contrast_fit, robust=T)
  dataset_exprs_top_genes <- topTable(dataset_exprs_contrast_ebayes, number=(dim(transposed_dataset_exprs)[1]), coef=1, adjust.method="fdr", confint=TRUE) 
  return(dataset_exprs_top_genes)
}

# run DE analysis on data
exprs_good_probes_EC_AD_CO_diag_exprs_top_genes<-run_diff_exprs_analysis(exprs_good_probes_EC_AD_CO_diag)

# add nuID
exprs_good_probes_EC_AD_CO_diag_exprs_top_genes$nuID <- rownames(exprs_good_probes_EC_AD_CO_diag_exprs_top_genes)

head(exprs_good_probes_EC_AD_CO_diag_exprs_top_genes)

fData(eset_bg_log2_rsn_combat_adj[rownames(exprs_good_probes_EC_AD_CO_diag_exprs_top_genes)[1:20],])[,1:5]

# merge fData with DE results
limma_ec_advco_inhouse <- merge(exprs_good_probes_EC_AD_CO_diag_exprs_top_genes,
            fData(eset_bg_log2_rsn_combat_adj),
            by.x="nuID",
            by.y="nuID",
            all.x=TRUE,
            sort=FALSE)

head(limma_ec_advco_inhouse)
dim(limma_ec_advco_inhouse)

#save DE results

setwd("/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/inhouse_brain/original_inhouse_limma_redone")
save(limma_ec_advco_inhouse, file="limma_ec_advco_inhouse.RData")

head(limma_ec_advco_inhouse)
dim(limma_ec_advco_inhouse)

sig_up<-subset(limma_ec_advco_inhouse, adj.P.Val<=0.01 & logFC>=0)
dim(sig_up)
head(sig_up)

length(unique(sig_up$TargetID))
as.data.frame(sig_up$TargetID)

setwd("/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/inhouse_brain/original_inhouse_limma_redone")

write.table(as.data.frame(sig_up$TargetID), row.names=F, col.names=F, quote=F, file="original_inhouse_p_0.01_up.txt")

##### CHECK SV ON DATA #####

# create function to check for SV in each dataset

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

# create function to adjust for SV in each dataset and return clean dataset - only run on dataset if has > 0 SV

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

# check SV in data

head(exprs_good_probes_EC_AD_CO_diag)[1:5]

check_SV_in_data(exprs_good_probes_EC_AD_CO_diag)

##### RUN LIMMA WITHOUT SVA ON RE_PROCESSED INHOUSSE DATA ####

# data taken from create_common_features_and_differential...... .Rdata

head(inhouse_data_Entorhinal_Cortex)[1:10]
inhouse_data_Entorhinal_Cortex_with_diag<-subset(inhouse_data_Entorhinal_Cortex, select = -c(1,2,3,4,6,7,8))
head(inhouse_data_Entorhinal_Cortex_with_diag)[1:5]

# check SVA

check_SV_in_data(inhouse_data_Entorhinal_Cortex_with_diag)

# run limma without SVA

re_processed_inhouse_data_Entorhinal_Cortex_top_genes<-run_diff_exprs_analysis(inhouse_data_Entorhinal_Cortex_with_diag)

head(re_processed_inhouse_data_Entorhinal_Cortex_top_genes)

# aconvert entrez gene id to gene symbol
entrez_gene_symbol <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(entrez_gene_symbol)
# Convert to a list
gene_symbol_lookup_table <- as.data.frame(entrez_gene_symbol[mapped_genes])

#subset to p<0.01 and logFC>0

re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up<-subset(re_processed_inhouse_data_Entorhinal_Cortex_top_genes, adj.P.Val<=0.01 & logFC>=0)
dim(re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up)

# extract entrez id only 
re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up_entrez_ID<-rownames(re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up)

#convert entrez id to gene symbol
re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up_gene_symbol<-gene_symbol_lookup_table[gene_symbol_lookup_table$gene_id %in% re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up_entrez_ID,]

head(re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up_gene_symbol)

setwd("/media/hamel/1TB/Projects/Brain_expression/1.Data/AD/inhouse_brain/original_inhouse_limma_redone")

write.table(as.data.frame(unique(re_processed_inhouse_data_Entorhinal_Cortex_top_genes_sig_up_gene_symbol$symbol)),
            row.names=F, 
            col.names=F, 
            quote=F, 
            file="re-processed_inhouse_p_0.01_up.txt")




