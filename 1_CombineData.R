rm(list=ls())
#--------------------------------------
# load libraries
## reading/writing data

library(readr) 

# tidyverse
library(dplyr)
#--------------------------------------

# define directory
setwd("D:/BCBL/Thalamus/Output")
# setwd("G:/Amepa/ABCDA/phenotypes/data/working_data/imaging/fs_output/")
if (dir.exists("figures")==FALSE) { dir.create("figures") }
if (dir.exists("summary_tables")==FALSE) { dir.create("summary_tables") }
#--------------------------------------

# read data
aseg <- read_csv("allsubjects_aseg_stats_vol_mm3.csv") #BCBL's aseg segmentation (output by Liu, Rodriguez, and Carrion)
thal <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv") #BCBL's thalamic nuclei volumes (output by Liu, Rodriguez, and Carrion)
ABCDv4_thalamus_data2merge_new <- readRDS("C:/Users/Alec/Downloads/ABCDv4_thalamus_data2merge_new.Rds") #output by ABCD dataset (version 4)
# ABCDv4_thalamus_data2merge_new <- readRDS("G:/Amepa/ABCDA/phenotypes/data/working_data/phenotypes/all/ABCDv4_thalamus_data2merge.Rds")
#--------------------------------------

#1. Editing column names
# rename column "subjectid" in "subject_session"
colnames(ABCDv4_thalamus_data2merge_new)[colnames(ABCDv4_thalamus_data2merge_new) == 'subjectid'] <- 'subject_session'

# Delete "_" to subject and timepoint columns
ABCDv4_thalamus_data2merge_new$subject_session <- gsub("_", "", ABCDv4_thalamus_data2merge_new$subject_session)
ABCDv4_thalamus_data2merge_new$timepoint <- gsub("_","", ABCDv4_thalamus_data2merge_new$timepoint)

# Add an underscore at the beginning of the items in timepoint column
# ABCDv4_thalamus_data2merge_new$timepoint <- paste0("_", ABCDv4_thalamus_data2merge_new$timepoint)

# merge names of column timepoint and subject_session as in "aseg" and "thal"
ABCDv4_thalamus_data2merge_new$subject_session <- paste(ABCDv4_thalamus_data2merge_new$subject_session, ABCDv4_thalamus_data2merge_new$timepoint,sep="_")

# fix upper case differences
ABCDv4_thalamus_data2merge_new$subject_session <- gsub("year", "Year", ABCDv4_thalamus_data2merge_new$subject_session) %>%
  gsub("followup", "FollowUp",.) %>%
  gsub("yarm1", "YArm1",.) %>%
  gsub("baselineyear", "baselineYear",.) %>%
  gsub("arm1", "Arm1", .)
#--------------------------------------

# 2. merge thal,aseg, and ABCD dataset
## merge information from aseg segmentation with thalamic nuclei segmentation volumes (BCBL's outputs)
data_1 <- merge(aseg, thal, by = "subject_session")
# generate new columns for subject and timepoint
data_1$subject <- strsplit(data_1$subject_session,"_") %>% sapply("[[",1)
data_1$timepoint <- strsplit(data_1$subject_session,"_") %>% sapply("[[",2) %>%   tolower()

# combine BCBL's outputs (data_1) with ABCD dataset into "data"
data <- merge(data_1, ABCDv4_thalamus_data2merge_new, by = c("subject_session","timepoint"),all.x=TRUE)

rm(data_1)

# 3. editing data
## replace "-" "(" ")" or " " with "."
colnames(data) <- make.names(colnames(data), unique = TRUE)

# create the new variables MD (mediodorsal and pulvinal, by combining their own subnuclei defined by Iglesias et al. 2018)
data <- data %>%
  mutate(Left.totalMD=rowSums(select(., contains("Left.MD"))),
         Right.totalMD=rowSums(select(., contains("Right.MD"))),
         #Left.totalPul=rowSums(select(., contains("Left.Pu"))),
         #Right.totalPul=rowSums(select(., contains("Right.Pu")))
         Left.totalPul = rowSums(select(., contains(c("Left.PuI", "Left.PuA", "Left.PuL", "Left.PuM")))),
         Right.totalPul = rowSums(select(., contains(c("Right.PuI", "Right.PuA", "Right.PuL", "Right.PuM"))))
  )
  
  
# define colnames for nuclei
thal_cols<-colnames(thal)[grep("subject|session|timepoint",colnames(thal),invert=TRUE)] %>% make.names()

# 4. save data
#save combined data 
saveRDS(data, file = "data_AMEPA.rds")













