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
# source general functions
# source("G:/Amepa/ABCDA/phenotypes/scripts/phenotypes/imaging/AlecGallo/clean/general_functions.R")
source("D:/BCBL/Thalamus/Output/general_functions.R") 
#--------------------------------------
# define directories and create if needed
# setwd("G:/Amepa/ABCDA/phenotypes/data/working_data/imaging/fs_output/")
setwd("D:/BCBL/Thalamus/Output/")
if (dir.exists("figures")==FALSE) { dir.create("figures") }
if (dir.exists("summary_tables")==FALSE) { dir.create("summary_tables") }

#--------------------------------------
# read combined data:
data<-readRDS("data_AMEPA.rds")
thal<-read.csv("thal.csv")

# thal_cols <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv",n_max = 0)[-1] %>% colnames()
thal_cols <- colnames(thal)[-1]

# get rid of "-" in colnames with "."
thal_cols<- make.names(thal_cols)

# define dependent variables
dv_cols<-colnames(data)[grep("nihtbx",colnames(data))]

# covariates
covs<- c("timepoint","sex","age",
        # "high.educ.bl","household.income.bl", 
        # "mri_info_device.serial.number" , 
        # "rel_family_id", 
        "abcd_site")

#--------------------------------------
# define all variables for which we want "raw" desriptive statistics
vars2check<-c(thal_cols,dv_cols,covs)
# define which ones are numeric and which are factor --> they will be treated differently when computing summary stats
n_vars<-vars2check[grep("numeric", data %>% dplyr::select(all_of(vars2check)) %>% sapply(.,class) )]
c_vars<-vars2check[grep("character|factor", data %>% dplyr::select(all_of(vars2check)) %>% sapply(.,class) )]
#--------------------------------------
#1. Descriptive Statistics of raw thalamic nuclei
descriptives <- data %>%
  group_by(timepoint) %>%
  skimr::skim_without_charts(all_of(vars2check))
colnames(descriptives)<-gsub("numeric.","",colnames(descriptives))

# summary_vars function to get all the descriptives we want, adds some additional columns such as skewness and kurtosis
descriptives_all<-summary_vars(d=data,n_vars=n_vars,c_vars=c_vars) %>% mutate(subset="all")
descriptives_t1<-summary_vars(d=data %>% filter(timepoint=="baselineyear1arm1"),n_vars=n_vars,c_vars=c_vars) %>%
  mutate(subset="baselineyear1arm1")
descriptives_t2<-summary_vars(d=data %>% filter(timepoint=="2yearfollowupyarm1"),n_vars=n_vars,c_vars=c_vars) %>% 
  mutate(subset="2yearfollowupyarm1")

descriptives2<-rbind(descriptives_all,descriptives_t1,descriptives_t2)
rm(descriptives_all,descriptives_t1,descriptives_t2)
# save
write.csv(descriptives,file="summary_tables/thalamicNuclei_descriptives.csv")
write.csv(descriptives2,file="summary_tables/thalamicNuclei_descriptives2.csv",row.names=FALSE)
# openxlsx::write.xlsx(descriptives,file="summary_tables/thalamicNuclei_descriptives.xlsx")


#########################################
##PLOTS
# Histograms of the raw thalamic nuclei and dvs
if (dir.exists("figures/histograms")==FALSE) { dir.create("figures/histograms") }
if (dir.exists("figures/histograms/rawnuclei")==FALSE) { dir.create("figures/histograms/rawnuclei") }

# add additional columns needed for plot 
descriptive_df <- descriptives %>% #filter(skim_variable %in% n_vars) %>% 
  # group_by(variable) %>%
  mutate(mean2sd_lower = mean - 2 * sd,
         mean2sd_upper = mean + 2 * sd,
         mean4sd_lower = mean - 4 * sd,
         mean4sd_upper = mean + 4 * sd,
         measure = skim_variable)

# in order to access the content of descriptive_df, you can now use "volume" for filtering.
for (vol in c(thal_cols)) {
  
  plot_name<-paste("dist",vol,"histogram.png",sep="_")
  
  if (!file.exists(  paste0("figures/histograms/rawnuclei/", plot_name) ) ){ 
  tmp<-subset(descriptive_df,measure==vol)
  p <- ggplot(data) + 
    geom_histogram(aes_string(vol,color="timepoint",fill="timepoint")) +
    geom_vline(data=tmp, aes(xintercept = mean,color=timepoint, fill=timepoint)) + # note that we don't need "string_aes" here, we can use "aes" directly because the column names are the same across volumes, not variable dependent.
    geom_vline(data=tmp, aes(xintercept = mean2sd_lower,color=timepoint),linetype="dashed") +
    geom_vline(data=tmp, aes(xintercept = mean2sd_upper,color=timepoint),linetype="dashed") +
    geom_vline(data=tmp, aes(xintercept = mean4sd_lower,color=timepoint),linetype="dotted") +
    geom_vline(data=tmp, aes(xintercept = mean4sd_upper,color=timepoint),linetype="dotted") +
    labs(title = paste0(vol, "_Distribution")) +
    theme(legend.position = "bottom",legend.direction = "horizontal") +
    theme(plot.title=element_text(family="Calibri",hjust = 0.5,size=16,face="bold")) +
    theme(axis.title.x=element_text(family = "Calibri",size=13)) +
    theme(axis.title.y=element_text(family = "Calibri",size=13)) +
    theme(legend.text=element_text(family="Calabri",size=11)) +
    #skewedness and kurtosis
    geom_text(aes(x = Inf, y = 50, label = paste0("Skewness= ", round(skewness(data[[vol]], na.rm = TRUE),2))),
              hjust = 1, size = 3, color = "black") + # x = Inf means put the lable at the very end
    geom_text(aes(x = Inf, y = 150, label = paste0("Kurtosis= ", round(kurtosis(data[[vol]], na.rm = TRUE),2))),
              hjust = 1, size = 3, color = "black") 
  # print(p)
  # you can also save the plots as independent objects

  # assign(plot_name,p)
  
  # save as figure
  ggsave(p, file = paste0("figures/histograms/rawnuclei/", plot_name, "_histogram.png"),
         width=5,height=4) # to edit/adjust
  # remove intermediate objects
  rm(p)
  }
  rm(plot_name)
}
#-------------------------------------------------

