
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                  DOWNLOAD E-GEOD-3790                                                                  #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# QC PIPELINE VERSION: 1.0
# DATE: 30/01/2017
# ARRAY EXPRESS NUMBER: E-GEOD-3790
# DISORDER: Huntingdon's Disease
# MICROARRAY PLATFORM: Affymetrix
# EXPRESSION CHIP: HG-U133A + HG-U133B
# NUMBER OF SAMPLES: 
# TISSUE: Cerebellum, Frontal Cortex, Caudate nucleus
#
# NOTES - 
# ONLY READING RAW DATA AND CREATEING R OBJECT
# 
# 

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/NON-AD-NEUROLOGICAL/Huntingdons_Disease/E-GEOD-3790"

setwd(work_dir)

# create directory for raw data
dir.create(paste(work_dir,"Raw_Data", sep="/"))
raw_dir=paste(work_dir,"Raw_Data", sep="/")

##### LOAD LIBRARIES ####

library(ArrayExpress)
library(affy)

##### DOWNLOAD RAW DATA #####

setwd(raw_dir)

#data_raw=getAE("E-GEOD-3790", type = "raw")

#data_raw=getAE("E-GEOD-3790", type = "processed")

data_raw=getAE("E-GEOD-3790", type = "full")

##### CREATE R EXPRESSION OBJECT FROM RAW DATA #####

# METHOD 1 - convert MAGE-TAB files into expresssion set  - USING RAW DATA

expression_data = ae2bioc(mageFiles = data_raw)

expression_data

# # METHOD 2 - convert MAGE-TAB files into expresssion set - USING RAW DATA
# 
# expression_data<-ReadAffy()
# 
# METHOD 3 - convert MAGE-TAB files into expresssion set - USING PROCESSED DATA
# 
# cnames=getcolproc(data_raw)
# 
# cnames
# 
# expression_data=procset(data_raw, cnames[2])
# 
# expression_data

##### SPLIT DATA INTO CHIPS ######

table((pData(expression_data[[1]]))$Factor.Value..DiseaseState.)
table((pData(expression_data[[1]]))$Characteristics..OrganismPart.)

table((pData(expression_data[[2]]))$Factor.Value..DiseaseState.)
table((pData(expression_data[[2]]))$Characteristics..OrganismPart.)

#create objects
Affy_U133A<-expression_data[[1]]

Affy_U133B<-expression_data[[2]]

##### BACKGROUND CORRECT DATA #####

Affy_U133A_bc<-mas5(Affy_U133A)

Affy_U133B_bc<-mas5(Affy_U133B)

##### CREATE SUB DIRECTORIES FOR EACH CHIP AN TISSUE AND SAVE #####

# create directory for each chip and write expression table
setwd(work_dir)

dir.create(paste(work_dir,"Affy_U133A", sep="/"))
Affy_U133A_dir=paste(work_dir,"Affy_U133A", sep="/")

dir.create(paste(work_dir,"Affy_U133B", sep="/"))
Affy_U133B_dir=paste(work_dir,"Affy_U133B", sep="/")

setwd(Affy_U133A_dir)
save(Affy_U133A_bc, file="E-GEOD-3790_Affy_U133A_bc.Rdata")

setwd(Affy_U133B_dir)
save(Affy_U133B_bc, file="E-GEOD-3790_Affy_U133B_bc.Rdata")

# create directory for cerebellum + frontal lobe in each chip

setwd(Affy_U133A_dir)
dir.create(paste(Affy_U133A_dir,"Cerebellum", sep="/"))
dir.create(paste(Affy_U133A_dir,"Frontal_Lobe", sep="/"))

setwd(Affy_U133B_dir)
dir.create(paste(Affy_U133B_dir,"Cerebellum", sep="/"))
dir.create(paste(Affy_U133B_dir,"Frontal_Lobe", sep="/"))

