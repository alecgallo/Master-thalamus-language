rm(list=ls())
#--------------------------------------
# load libraries

## basic, should not be necessary
library(base)
library(stats)

# reading/writing data
library(readr) 
library(writexl)

# tidyverse
library(dplyr)
library(tidyr)
library(tidyverse)
# stats, kurtosis and skewness... currently not used?
library(moments)

# plotting
library(ggplot2)
library(cowplot)
library(ggpubr)
#--------------------------------------

# define directories and create if needed
# setwd("G:/Amepa/ABCDA/phenotypes/data/working_data/imaging/fs_output/")
setwd("D:/BCBL/Thalamus/Output/")
if (dir.exists("figures")==FALSE) { dir.create("figures") }
if (dir.exists("summary_tables")==FALSE) { dir.create("summary_tables") }
#--------------------------------------
# read combined data:
data<-readRDS("data_AMEPA_Normalized.rds") 
data_raw<-readRDS("data_AMEPA.rds")
thal<-readRDS("thal_AMEPA.rds")

thal_cols <- colnames(thal)[-1]
# get rid of "-" in colnames with "."
thal_cols<- make.names(thal_cols)
thal_cols_of_interest<-thal_cols[grep("MD|Pul|LGN|MGN|AV",thal_cols)] 

# define dependent variables
dv_cols<-colnames(data)[grep("nihtbx",colnames(data))]

#------------------

#5- QC: we would apply the thresholds for defined in descriptive stats here. --> generate a clean dataset that we can use for modelling. 
#--> we could also discard all the variables that we are not interested in (e.g. the other aseg volumes?) so that we focus only on the thalamus.

#apply the threshold to define outliers (outlier value = mean ± 4sd)

outliers_df <- data.frame()

# it loops over the thalamic volumes
for (vol in c(thal_cols,dv_cols, paste0("adj_",thal_cols))) { # ACC: I would only use thal_cols of interest and dvs here
  
  # loops over the timepoints
  for (timepoint in c("baselineyear1arm1", "2yearfollowupyarm1")) {
    
    # subset the data based on the volume and the tiempoint (baseline and 2yearsfollowup)
    subset_data <- data[data$timepoint == timepoint, vol]
    
    if (any(is.na(subset_data))) { # check if vol_data contains NA
      subset_data <- na.omit(subset_data) # remove missing values
    }
    
    # outlier stats
    mean_val <- mean(subset_data)
    sd_val <- sd(subset_data)
    cutoff_upper <- mean_val + 4 * sd_val
    cutoff_lower <- mean_val - 4 * sd_val
    outliers <- subset(subset_data, subset_data > cutoff_upper | subset_data < cutoff_lower) #outlier values 
    
    if (length(outliers) == 0) { # if outlier vector is empty, skip to the next vol
      next
    }
    
    # outlier-subjects
    #AG: I REALIZED THIS IS WRONG (ERROR!) because it's giving the row number of the subset_data rather than data
    outlier_corresp <- which(subset_data > cutoff_upper | subset_data < cutoff_lower) # row corresponding to the outlier values in the df "data"
    outlier_participants <- data$subject_session[data$timepoint == timepoint][outlier_corresp] # participants related to the outlier values
    
    # outlier related to which participant
    outlier2participants <- data.frame(vol, timepoint, cutoff_lower, cutoff_upper,outliers, outlier_corresp, outlier_participants)
    
    # ACC if outliers_df does not exist, just take outlier2participants, otherwise combine
    if(!exists("outliers_df")){
      outliers_df <- outlier2participants
    } else {
      outliers_df <- rbind(outliers_df, outlier2participants)
    }
    
    rm(cutoff_lower,cutoff_upper,mean_val,sd_val)
    rm(outlier2participants, outlier_participants, outlier_corresp, outliers)
    
  }
  
}

outliers_df <- outliers_df %>% mutate(id=outlier_participants %>% strsplit(.,"_") %>% sapply("[[",1)) %>%
  mutate(var_of_interest= if_else(vol %in% thal_cols_of_interest,1,0)) %>%
  dplyr::select(-outlier_corresp)

# check how many outliers per participant
table(outliers_df$outlier_participants) %>% table()
table(outliers_df$id) %>% table()

write_xlsx(outliers_df, 'summary_tables/outliers_table.xlsx')
##############################################################
outliers_df2<-subset(outliers_df,var_of_interest==1)
# examine whether a participant occurs frequently in the outliers_df
participant_counts <- table(outliers_df2$outlier_participants)

most_common_participant <- names(participant_counts)[which.max(participant_counts)] # "NDARINVFMHWAFKJ_2YearFollowUpYArm1"
most_common_participant_count <- participant_counts[most_common_participant] #occurs 6 times


## Define the clean dataset (thalamic nuclei, thalamic volumes, dependent variables, coviariate)
# clear the data off from the outlier values

# subset considering the thalamic nuclei
data_clean <- select(data, subject, # ACC: include subject ID
                     subject_session, age, timepoint, sex,                            #subjects, age, timepoints, and sex to birth
                     nihtbx_picvocab_uncorrected, nihtbx_reading_uncorrected,               #dependent variables: vocabulary and reading values
                     thal_cols # ACC: these ones will be used in the regression analyses
                     # c(paste0("adj_",thal_cols))
                     )                                           #thalamic volumes and nuclei volumes (thalamic volume adjusted)
                     #c(thal_cols, paste0("ICV_adj_))



# delete outliers (delete each participant having even a single outlier)

participants_outlier <- unique(outliers_df$outlier_participants) # 480
participants_outlier <- unique(outliers_df$id) # 426
#exclude participants_outlier from the data_clean
data_cleaned <- data_clean[!data_clean$subject %in% participants_outlier, ] # 10656 - 408 = 10248 

# ACC: make sure all subjects have both datapoints
subjects2timepoints<-names(which(table(data_cleaned$subject)==2))
data_cleaned <- data_clean[data_clean$subject %in% subjects2timepoints, ] # 
table(table(data_cleaned$subject)) # ACC: 4853 subjects with 2 datapoints each

#save the cleaned dataset
saveRDS(data_cleaned, file = "data_AMEPA_cleaned.rds")


########################################################
#------------------------------------------------------
## computing AI
#------------------------------------------------------
# compute AIs for each roi *could move down...*
computeAI<-function(lh,rh){
  ai=(lh-rh)/((lh+rh)/2)
}

rois4ai<-gsub("Left.|Right.","",thal_cols) %>% unique()
# compute AI for each ROI
d_ai<-sapply(rois4ai,function(r){
  lh=d[,paste0(r,"Left.")]
  rh=d[,paste0(r,"Right.")]
  ai=computeAI(lh,rh)
  return(ai)
})
colnames(d_ai)<-paste0(colnames(d_ai),".ai")
# combine with data
d<-cbind(d,d_ai)
rm(d_ai)

#------------------------------------------------------
## z-scores
#------------------------------------------------------
base_dir <- "D:/BCBL/Thalamus/Output/summary_tables/"
# scale (z-transform)
data<-data %>% mutate_if(is.numeric,scale)
f=paste(base_dir,"scaled4analysis.csv",sep="")
f2=gsub(".csv",".RDS",f)
write.csv(data,f,row.names=FALSE)


