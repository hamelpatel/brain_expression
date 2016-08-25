##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                         ANALYSE CMAP RESULTS                                                           #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################


# DESCRIPTION
# READ IN CMAP RESULTS
# DOES daph APPEAR AS TOP HIT IN ENTORIHNAL CORTEX AND REALTED REGIONS - 
# COMPILE LIST OF COMMON DRUGS PER BRAIN REGION
# 

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

cmap_down_dir="/media/hamel/1TB/Projects/Brain_expression/5.cMAP_result_analysis/EnrichR_output/cMAP_down/"

cmap_up_dir="/media/hamel/1TB/Projects/Brain_expression/5.cMAP_result_analysis/EnrichR_output/cMAP_up/"

work_dir="/media/hamel/1TB/Projects/Brain_expression/5.cMAP_result_analysis"

##### READ IN CMAP RESULTS #####

# convert file names using bash to replace hyphens to underscore - run on command line within directory
#   
# for i in `ls *-*`
# do
# NEW=`echo $i|tr '-' '_'`
# mv $i $NEW
# 
# done

setwd(cmap_down_dir)

list.files()

files_to_read<-list.files()

for (i in 1:length(files_to_read)) assign(files_to_read[i], read.csv(files_to_read[i], header=T, sep="\t"))

setwd(cmap_up_dir)

list.files()

files_to_read<-list.files()

for (i in 1:length(files_to_read)) assign(files_to_read[i], read.csv(files_to_read[i], header=T, sep="\t"))

ls()

head(inhouse_Entorhinal_Cortex_p0.05_up.txt[order(inhouse_Entorhinal_Cortex_p0.05_up.txt$P.value),])
head(inhouse_Entorhinal_Cortex_p0.05_up.txt)

head(inhouse_Entorhinal_Cortex_p0.05_down.txt[order(inhouse_Entorhinal_Cortex_p0.05_down.txt$P.value),])
head()

##### SORT DATAFRAMES BY P.values #####

AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt<-AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt[order(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt$P.value),]
AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt<-AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt[order(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt$P.value),]
AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt<-AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt[order(AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt$P.value),]
AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt<-AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt[order(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt$P.value),]
E_GEOD_1297_Hippocampus_p0.05_up.txt<-E_GEOD_1297_Hippocampus_p0.05_up.txt[order(E_GEOD_1297_Hippocampus_p0.05_up.txt$P.value),]
E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt<-E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt[order(E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt$P.value),]
E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt<-E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt[order(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt$P.value),]
E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt<-E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt[order(E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt$P.value),]
E_GEOD_36980_Frontal_Cortex_p0.05_up.txt<-E_GEOD_36980_Frontal_Cortex_p0.05_up.txt[order(E_GEOD_36980_Frontal_Cortex_p0.05_up.txt$P.value),]
E_GEOD_36980_Hippocampus_p0.05_up.txt<-E_GEOD_36980_Hippocampus_p0.05_up.txt[order(E_GEOD_36980_Hippocampus_p0.05_up.txt$P.value),]
E_GEOD_36980_Temporal_Cortex_p0.05_up.txt<-E_GEOD_36980_Temporal_Cortex_p0.05_up.txt[order(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt$P.value),]
E_GEOD_48350_Hippocampus_p0.05_up.txt<-E_GEOD_48350_Hippocampus_p0.05_up.txt[order(E_GEOD_48350_Hippocampus_p0.05_up.txt$P.value),]
E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt<-E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt[order(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt$P.value),]
E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt<-E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt[order(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt$P.value),]
E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt<-E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt[order(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt$P.value),]
E_GEOD_5281_Hippocampus_p0.05_up.txt<-E_GEOD_5281_Hippocampus_p0.05_up.txt[order(E_GEOD_5281_Hippocampus_p0.05_up.txt$P.value),]
E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt<-E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt[order(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt$P.value),]
E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt<-E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt[order(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt$P.value),]
E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt<-E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt[order(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt$P.value),]
E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt<-E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt[order(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt$P.value),]
inhouse_Cerebellum_p0.05_up.txt<-inhouse_Cerebellum_p0.05_up.txt[order(inhouse_Cerebellum_p0.05_up.txt$P.value),]
inhouse_Entorhinal_Cortex_p0.05_up.txt<-inhouse_Entorhinal_Cortex_p0.05_up.txt[order(inhouse_Entorhinal_Cortex_p0.05_up.txt$P.value),]
inhouse_Frontal_Cortex_p0.05_up.txt<-inhouse_Frontal_Cortex_p0.05_up.txt[order(inhouse_Frontal_Cortex_p0.05_up.txt$P.value),]
inhouse_Temporal_Cortex_p0.05_up.txt<-inhouse_Temporal_Cortex_p0.05_up.txt[order(inhouse_Temporal_Cortex_p0.05_up.txt$P.value),]

AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt<-AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt[order(AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt$P.value),]
AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt<-AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt[order(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt$P.value),]
AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt<-AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt[order(AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt$P.value),]
E_GEOD_1297_Hippocampus_p0.05_down.txt<-E_GEOD_1297_Hippocampus_p0.05_down.txt[order(E_GEOD_1297_Hippocampus_p0.05_down.txt$P.value),]
E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt<-E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt[order(E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt$P.value),]
E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt<-E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt[order(E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt$P.value),]
E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt<-E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt[order(E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt$P.value),]
E_GEOD_36980_Frontal_Cortex_p0.05_down.txt<-E_GEOD_36980_Frontal_Cortex_p0.05_down.txt[order(E_GEOD_36980_Frontal_Cortex_p0.05_down.txt$P.value),]
E_GEOD_36980_Hippocampus_p0.05_down.txt<-E_GEOD_36980_Hippocampus_p0.05_down.txt[order(E_GEOD_36980_Hippocampus_p0.05_down.txt$P.value),]
E_GEOD_36980_Temporal_Cortex_p0.05_down.txt<-E_GEOD_36980_Temporal_Cortex_p0.05_down.txt[order(E_GEOD_36980_Temporal_Cortex_p0.05_down.txt$P.value),]
E_GEOD_48350_Hippocampus_p0.05_down.txt<-E_GEOD_48350_Hippocampus_p0.05_down.txt[order(E_GEOD_48350_Hippocampus_p0.05_down.txt$P.value),]
E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt<-E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt[order(E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt$P.value),]
E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt<-E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt[order(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt$P.value),]
E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt<-E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt[order(E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt$P.value),]
E_GEOD_5281_Hippocampus_p0.05_down.txt<-E_GEOD_5281_Hippocampus_p0.05_down.txt[order(E_GEOD_5281_Hippocampus_p0.05_down.txt$P.value),]
E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt<-E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt[order(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt$P.value),]
E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt<-E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt[order(E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt$P.value),]
E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt<-E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt[order(E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt$P.value),]
E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt<-E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt[order(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt$P.value),]
inhouse_Cerebellum_p0.05_down.txt<-inhouse_Cerebellum_p0.05_down.txt[order(inhouse_Cerebellum_p0.05_down.txt$P.value),]
inhouse_Entorhinal_Cortex_p0.05_down.txt<-inhouse_Entorhinal_Cortex_p0.05_down.txt[order(inhouse_Entorhinal_Cortex_p0.05_down.txt$P.value),]
inhouse_Frontal_Cortex_p0.05_down.txt<-inhouse_Frontal_Cortex_p0.05_down.txt[order(inhouse_Frontal_Cortex_p0.05_down.txt$P.value),]
inhouse_Temporal_Cortex_p0.05_down.txt<-inhouse_Temporal_Cortex_p0.05_down.txt[order(inhouse_Temporal_Cortex_p0.05_down.txt$P.value),]

# rename rownames to numeric

rownames(AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt)<-1:dim(AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt)[1]
rownames(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt)<-1:dim(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt)[1]
rownames(AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt)<-1:dim(AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt)[1]
rownames(E_GEOD_1297_Hippocampus_p0.05_down.txt)<-1:dim(E_GEOD_1297_Hippocampus_p0.05_down.txt)[1]
rownames(E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt)<-1:dim(E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt)[1]
rownames(E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt)<-1:dim(E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt)[1]
rownames(E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt)<-1:dim(E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt)[1]
rownames(E_GEOD_36980_Frontal_Cortex_p0.05_down.txt)<-1:dim(E_GEOD_36980_Frontal_Cortex_p0.05_down.txt)[1]
rownames(E_GEOD_36980_Hippocampus_p0.05_down.txt)<-1:dim(E_GEOD_36980_Hippocampus_p0.05_down.txt)[1]
rownames(E_GEOD_36980_Temporal_Cortex_p0.05_down.txt)<-1:dim(E_GEOD_36980_Temporal_Cortex_p0.05_down.txt)[1]
rownames(E_GEOD_48350_Hippocampus_p0.05_down.txt)<-1:dim(E_GEOD_48350_Hippocampus_p0.05_down.txt)[1]
rownames(E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt)<-1:dim(E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt)[1]
rownames(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt)<-1:dim(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt)<-1:dim(E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Hippocampus_p0.05_down.txt)<-1:dim(E_GEOD_5281_Hippocampus_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt)<-1:dim(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt)<-1:dim(E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt)<-1:dim(E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt)[1]
rownames(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt)<-1:dim(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt)[1]
rownames(inhouse_Cerebellum_p0.05_down.txt)<-1:dim(inhouse_Cerebellum_p0.05_down.txt)[1]
rownames(inhouse_Entorhinal_Cortex_p0.05_down.txt)<-1:dim(inhouse_Entorhinal_Cortex_p0.05_down.txt)[1]
rownames(inhouse_Frontal_Cortex_p0.05_down.txt)<-1:dim(inhouse_Frontal_Cortex_p0.05_down.txt)[1]
rownames(inhouse_Temporal_Cortex_p0.05_down.txt)<-1:dim(inhouse_Temporal_Cortex_p0.05_down.txt)[1]

rownames(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt)<-1:dim(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt)[1]
rownames(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt)<-1:dim(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt)[1]
rownames(AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt)<-1:dim(AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt)[1]
rownames(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt)<-1:dim(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt)[1]
rownames(E_GEOD_1297_Hippocampus_p0.05_up.txt)<-1:dim(E_GEOD_1297_Hippocampus_p0.05_up.txt)[1]
rownames(E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt)<-1:dim(E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt)[1]
rownames(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt)<-1:dim(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt)[1]
rownames(E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt)<-1:dim(E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt)[1]
rownames(E_GEOD_36980_Frontal_Cortex_p0.05_up.txt)<-1:dim(E_GEOD_36980_Frontal_Cortex_p0.05_up.txt)[1]
rownames(E_GEOD_36980_Hippocampus_p0.05_up.txt)<-1:dim(E_GEOD_36980_Hippocampus_p0.05_up.txt)[1]
rownames(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt)<-1:dim(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt)[1]
rownames(E_GEOD_48350_Hippocampus_p0.05_up.txt)<-1:dim(E_GEOD_48350_Hippocampus_p0.05_up.txt)[1]
rownames(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt)<-1:dim(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt)[1]
rownames(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt)<-1:dim(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt)<-1:dim(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Hippocampus_p0.05_up.txt)<-1:dim(E_GEOD_5281_Hippocampus_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt)<-1:dim(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt)<-1:dim(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt)<-1:dim(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt)[1]
rownames(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt)<-1:dim(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt)[1]
rownames(inhouse_Cerebellum_p0.05_up.txt)<-1:dim(inhouse_Cerebellum_p0.05_up.txt)[1]
rownames(inhouse_Entorhinal_Cortex_p0.05_up.txt)<-1:dim(inhouse_Entorhinal_Cortex_p0.05_up.txt)[1]
rownames(inhouse_Frontal_Cortex_p0.05_up.txt)<-1:dim(inhouse_Frontal_Cortex_p0.05_up.txt)[1]
rownames(inhouse_Temporal_Cortex_p0.05_up.txt)<-1:dim(inhouse_Temporal_Cortex_p0.05_up.txt)[1]

head(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt)
head(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt)
head(AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt)
head(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt)
head(E_GEOD_1297_Hippocampus_p0.05_up.txt)
head(E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt)
head(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt)
head(E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt)
head(E_GEOD_36980_Frontal_Cortex_p0.05_up.txt)
head(E_GEOD_36980_Hippocampus_p0.05_up.txt)
head(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt)
head(E_GEOD_48350_Hippocampus_p0.05_up.txt)
head(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt)
head(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt)
head(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt)
head(E_GEOD_5281_Hippocampus_p0.05_up.txt)
head(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt)
head(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt)
head(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt)
head(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt)
head(inhouse_Cerebellum_p0.05_up.txt)
head(inhouse_Entorhinal_Cortex_p0.05_up.txt)
head(inhouse_Frontal_Cortex_p0.05_up.txt)
head(inhouse_Temporal_Cortex_p0.05_up.txt)

head(AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt)
head(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt)
head(AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt)
head(E_GEOD_1297_Hippocampus_p0.05_down.txt)
head(E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt)
head(E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt)
head(E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt)
head(E_GEOD_36980_Frontal_Cortex_p0.05_down.txt)
head(E_GEOD_36980_Hippocampus_p0.05_down.txt)
head(E_GEOD_36980_Temporal_Cortex_p0.05_down.txt)
head(E_GEOD_48350_Hippocampus_p0.05_down.txt)
head(E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt)
head(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt)
head(E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt)
head(E_GEOD_5281_Hippocampus_p0.05_down.txt)
head(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt)
head(E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt)
head(E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt)
head(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt)
head(inhouse_Cerebellum_p0.05_down.txt)
head(inhouse_Entorhinal_Cortex_p0.05_down.txt)
head(inhouse_Frontal_Cortex_p0.05_down.txt)
head(inhouse_Temporal_Cortex_p0.05_down.txt)

##### DAPH p Value in all datasets #####

AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt_DAPH<-AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt[grep("4,5-dianilinophthalimide",AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt$Term),]
AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt_DAPH<-AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt$Term),]
AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt_DAPH<-AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt[grep("4,5-dianilinophthalimide",AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt$Term),]
E_GEOD_1297_Hippocampus_p0.05_down.txt_DAPH<-E_GEOD_1297_Hippocampus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_1297_Hippocampus_p0.05_down.txt$Term),]
E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt_DAPH<-E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt$Term),]
E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt_DAPH<-E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt$Term),]
E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt_DAPH<-E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt$Term),]
E_GEOD_36980_Frontal_Cortex_p0.05_down.txt_DAPH<-E_GEOD_36980_Frontal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Frontal_Cortex_p0.05_down.txt$Term),]
E_GEOD_36980_Hippocampus_p0.05_down.txt_DAPH<-E_GEOD_36980_Hippocampus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Hippocampus_p0.05_down.txt$Term),]
E_GEOD_36980_Temporal_Cortex_p0.05_down.txt_DAPH<-E_GEOD_36980_Temporal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Temporal_Cortex_p0.05_down.txt$Term),]
E_GEOD_48350_Hippocampus_p0.05_down.txt_DAPH<-E_GEOD_48350_Hippocampus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Hippocampus_p0.05_down.txt$Term),]
E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt_DAPH<-E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt$Term),]
E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt_DAPH<-E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt$Term),]
E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt_DAPH<-E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt$Term),]
E_GEOD_5281_Hippocampus_p0.05_down.txt_DAPH<-E_GEOD_5281_Hippocampus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Hippocampus_p0.05_down.txt$Term),]
E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt_DAPH<-E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt$Term),]
E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt_DAPH<-E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt$Term),]
E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt_DAPH<-E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt$Term),]
E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt_DAPH<-E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt$Term),]
inhouse_Cerebellum_p0.05_down.txt_DAPH<-inhouse_Cerebellum_p0.05_down.txt[grep("4,5-dianilinophthalimide",inhouse_Cerebellum_p0.05_down.txt$Term),]
inhouse_Entorhinal_Cortex_p0.05_down.txt_DAPH<-inhouse_Entorhinal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",inhouse_Entorhinal_Cortex_p0.05_down.txt$Term),]
inhouse_Frontal_Cortex_p0.05_down.txt_DAPH<-inhouse_Frontal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",inhouse_Frontal_Cortex_p0.05_down.txt$Term),]
inhouse_Temporal_Cortex_p0.05_down.txt_DAPH<-inhouse_Temporal_Cortex_p0.05_down.txt[grep("4,5-dianilinophthalimide",inhouse_Temporal_Cortex_p0.05_down.txt$Term),]

AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt_DAPH<-AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt[grep("4,5-dianilinophthalimide",AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt$Term),]
AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt_DAPH<-AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt$Term),]
AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt_DAPH<-AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt[grep("4,5-dianilinophthalimide",AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt$Term),]
AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt_DAPH<-AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide",AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt$Term),]
E_GEOD_1297_Hippocampus_p0.05_up.txt_DAPH<-E_GEOD_1297_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_1297_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt_DAPH<-E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt$Term),]
E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt_DAPH<-E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt$Term),]
E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt_DAPH<-E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt$Term),]
E_GEOD_36980_Frontal_Cortex_p0.05_up.txt_DAPH<-E_GEOD_36980_Frontal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Frontal_Cortex_p0.05_up.txt$Term),]
E_GEOD_36980_Hippocampus_p0.05_up.txt_DAPH<-E_GEOD_36980_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_36980_Temporal_Cortex_p0.05_up.txt_DAPH<-E_GEOD_36980_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_36980_Temporal_Cortex_p0.05_up.txt$Term),]
E_GEOD_48350_Hippocampus_p0.05_up.txt_DAPH<-E_GEOD_48350_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt_DAPH<-E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt$Term),]
E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH<-E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt$Term),]
E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt_DAPH<-E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt$Term),]
E_GEOD_5281_Hippocampus_p0.05_up.txt_DAPH<-E_GEOD_5281_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt_DAPH<-E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt$Term),]
E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt_DAPH<-E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt$Term),]
E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt_DAPH<-E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt$Term),]
E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH<-E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide",E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt$Term),]
inhouse_Cerebellum_p0.05_up.txt_DAPH<-inhouse_Cerebellum_p0.05_up.txt[grep("4,5-dianilinophthalimide",inhouse_Cerebellum_p0.05_up.txt$Term),]
inhouse_Entorhinal_Cortex_p0.05_up.txt_DAPH<-inhouse_Entorhinal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",inhouse_Entorhinal_Cortex_p0.05_up.txt$Term),]
inhouse_Frontal_Cortex_p0.05_up.txt_DAPH<-inhouse_Frontal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",inhouse_Frontal_Cortex_p0.05_up.txt$Term),]
inhouse_Temporal_Cortex_p0.05_up.txt_DAPH<-inhouse_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide",inhouse_Temporal_Cortex_p0.05_up.txt$Term),]

# keep only dataframes with info

# DAPH up
head(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt_DAPH)
head(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt_DAPH)
head(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt_DAPH)
head(E_GEOD_29378_Hippocampus_CA3_p0.05_up.txt_DAPH)
head(E_GEOD_36980_Hippocampus_p0.05_up.txt_DAPH)
head(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt_DAPH)
head(E_GEOD_48350_Hippocampus_p0.05_up.txt_DAPH)
head(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt_DAPH)
head(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Hippocampus_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt_DAPH)
head(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH)
head(inhouse_Cerebellum_p0.05_up.txt_DAPH)
head(inhouse_Entorhinal_Cortex_p0.05_up.txt_DAPH)
head(inhouse_Frontal_Cortex_p0.05_up.txt_DAPH)
head(inhouse_Temporal_Cortex_p0.05_up.txt_DAPH)

# DAPH down - found
head(AMP_MAYOeGWAS_Cerebellum_p0.05_down.txt_DAPH)
head(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_down.txt_DAPH)
head(E_GEOD_29378_Hippocampus_CA1_p0.05_down.txt_DAPH)
head(E_GEOD_29378_Hippocampus_CA3_p0.05_down.txt_DAPH)
head(E_GEOD_36980_Frontal_Cortex_p0.05_down.txt_DAPH)
head(E_GEOD_36980_Hippocampus_p0.05_down.txt_DAPH)
head(E_GEOD_36980_Temporal_Cortex_p0.05_down.txt_DAPH)
head(E_GEOD_48350_Postcentral_Gyrus_p0.05_down.txt_DAPH)
head(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Entorhinal_Cortex_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Hippocampus_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Posterior_Cingulate_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Primary_Visual_Cortex_p0.05_down.txt_DAPH)
head(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_down.txt_DAPH)
head(inhouse_Cerebellum_p0.05_down.txt_DAPH)
head(inhouse_Entorhinal_Cortex_p0.05_down.txt_DAPH)
head(inhouse_Frontal_Cortex_p0.05_down.txt_DAPH)
head(inhouse_Temporal_Cortex_p0.05_down.txt_DAPH)

#DAPH up - not found
head(AMP_MSBB_U133A_Frontal_Pole_p0.05_up.txt_DAPH)
head(AMP_MSBB_U133A_Inferior_Temporal_Gyrus_p0.05_up.txt_DAPH)
head(E_GEOD_1297_Hippocampus_p0.05_up.txt_DAPH)
head(E_GEOD_28146_Hippocampus_CA1_p0.05_up.txt_DAPH)
head(E_GEOD_36980_Frontal_Cortex_p0.05_up.txt_DAPH)

#DAPH down - not found
head(AMP_MSBB_U133A_Frontal_Pole_p0.05_down.txt_DAPH)
head(E_GEOD_1297_Hippocampus_p0.05_down.txt_DAPH)
head(E_GEOD_28146_Hippocampus_CA1_p0.05_down.txt_DAPH)
head(E_GEOD_48350_Hippocampus_p0.05_down.txt_DAPH)

##### DAPH PLOT ####

#extract DAPH 624 from only dataframes with it
AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt_DAPH_624<-AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt$Term),]
AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt_DAPH_624<-AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt$Term),]
E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt_DAPH_624<-E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt$Term),]
E_GEOD_36980_Hippocampus_p0.05_up.txt_DAPH_624<-E_GEOD_36980_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_36980_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_36980_Temporal_Cortex_p0.05_up.txt_DAPH_624<-E_GEOD_36980_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_36980_Temporal_Cortex_p0.05_up.txt$Term),]
E_GEOD_48350_Hippocampus_p0.05_up.txt_DAPH_624<-E_GEOD_48350_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_48350_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt_DAPH_624<-E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt$Term),]
E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624<-E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt$Term),]
E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt$Term),]
E_GEOD_5281_Hippocampus_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Hippocampus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Hippocampus_p0.05_up.txt$Term),]
E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt$Term),]
E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt$Term),]
E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt$Term),]
E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624<-E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt$Term),]
inhouse_Entorhinal_Cortex_p0.05_up.txt_DAPH_624<-inhouse_Entorhinal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",inhouse_Entorhinal_Cortex_p0.05_up.txt$Term),]
inhouse_Frontal_Cortex_p0.05_up.txt_DAPH_624<-inhouse_Frontal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",inhouse_Frontal_Cortex_p0.05_up.txt$Term),]
inhouse_Temporal_Cortex_p0.05_up.txt_DAPH_624<-inhouse_Temporal_Cortex_p0.05_up.txt[grep("4,5-dianilinophthalimide-624",inhouse_Temporal_Cortex_p0.05_up.txt$Term),]

# change rownames
rownames(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt_DAPH_624)<-"AMP_MAYOeGWAS_Cerebellum"
rownames(AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt_DAPH_624)<-"AMP_MAYOeGWAS_Temporal_Cortex"
rownames(E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt_DAPH_624)<-"E_GEOD_29378_Hippocampus_CA1"
rownames(E_GEOD_36980_Hippocampus_p0.05_up.txt_DAPH_624)<-"E_GEOD_36980_Hippocampus"
rownames(E_GEOD_36980_Temporal_Cortex_p0.05_up.txt_DAPH_624)<-"E_GEOD_36980_Temporal_Cortex"
rownames(E_GEOD_48350_Hippocampus_p0.05_up.txt_DAPH_624)<-"E_GEOD_48350_Hippocampus"
rownames(E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt_DAPH_624)<-"E_GEOD_48350_Postcentral_Gyrus"
rownames(E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624)<-"E_GEOD_48350_Superior_Frontal_Gyrus"
rownames(E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Entorhinal_Cortex"
rownames(E_GEOD_5281_Hippocampus_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Hippocampus"
rownames(E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Medial_Temporal_Gyrus"
rownames(E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Posterior_Cingulate"
rownames(E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Primary_Visual_Cortex"
rownames(E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624)<-"E_GEOD_5281_Superior_Frontal_Gyrus"
rownames(inhouse_Entorhinal_Cortex_p0.05_up.txt_DAPH_624)<-"inhouse_Entorhinal_Cortex"
rownames(inhouse_Frontal_Cortex_p0.05_up.txt_DAPH_624)<-"inhouse_Frontal_Cortex"
rownames(inhouse_Temporal_Cortex_p0.05_up.txt_DAPH_624)<-"inhouse_Temporal_Cortex"

# merge
up_regulated_genes_DAPH<-rbind(AMP_MAYOeGWAS_Cerebellum_p0.05_up.txt_DAPH_624,
                               AMP_MAYOeGWAS_Temporal_Cortex_p0.05_up.txt_DAPH_624,
                               E_GEOD_29378_Hippocampus_CA1_p0.05_up.txt_DAPH_624,
                               E_GEOD_36980_Hippocampus_p0.05_up.txt_DAPH_624,
                               E_GEOD_36980_Temporal_Cortex_p0.05_up.txt_DAPH_624,
                               E_GEOD_48350_Hippocampus_p0.05_up.txt_DAPH_624,
                               E_GEOD_48350_Postcentral_Gyrus_p0.05_up.txt_DAPH_624,
                               E_GEOD_48350_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Entorhinal_Cortex_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Hippocampus_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Medial_Temporal_Gyrus_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Posterior_Cingulate_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Primary_Visual_Cortex_p0.05_up.txt_DAPH_624,
                               E_GEOD_5281_Superior_Frontal_Gyrus_p0.05_up.txt_DAPH_624,
                               inhouse_Entorhinal_Cortex_p0.05_up.txt_DAPH_624,
                               inhouse_Frontal_Cortex_p0.05_up.txt_DAPH_624,
                               inhouse_Temporal_Cortex_p0.05_up.txt_DAPH_624)

#order by P.val
up_regulated_genes_DAPH<-up_regulated_genes_DAPH[order(up_regulated_genes_DAPH$P.value),]
head(up_regulated_genes_DAPH)

#plot
par(mar=c(16,5,5,1))
barplot(t(up_regulated_genes_DAPH[3]), las=2, cex.names=0.75, main="DAPH P value in upregulated DE genes", ylab="P Value")
abline(h=0.05)
mtext("Brain Region", side=1, line=13.5)

setwd(work_dir)

pdf("Up_regulated_DE_genes_DAPH_624_barplot.pdf")
par(mar=c(16,5,5,1))
barplot(t(up_regulated_genes_DAPH[3]), las=2, cex.names=0.75, main="DAPH P value in upregulated DE genes", ylab="P Value")
abline(h=0.05)
mtext("Brain Region", side=1, line=13.5)
dev.off()

##### SAVE IMAGE #####

setwd(work_dir)
save.image("Analyse_cMAP_results.Rdata")
