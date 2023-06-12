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
# Remember to load the "general.functions.R" before running the code

#--------------------------------------
# define directories and create if needed
# setwd("G:/Amepa/ABCDA/phenotypes/data/working_data/imaging/fs_output/")
setwd("D:/BCBL/Thalamus/Output/")
if (dir.exists("figures")==FALSE) { dir.create("figures") }
if (dir.exists("summary_tables")==FALSE) { dir.create("summary_tables") }
#--------------------------------------
# read combined data:
data<-readRDS("data_AMEPA_Normalized.rds")
aseg <- read_csv("allsubjects_aseg_stats_vol_mm3.csv") #BCBL's aseg segmentation (output by Liu, Rodriguez, and Carrion)
thal<-read.csv("thal.csv")
# thal <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv") #BCBL's thalamic nuclei volumes (output by Liu, Rodriguez, and Carrion)

thal_cols <- colnames(thal)[-1]

# # get rid of "-" in colnames with "."
# thal_cols<- make.names(thal_cols)
#------------------

# Create descriptive statistics (maybe we can delete one of them: The skimr is finally working on my computer!)
# define aseg cols
aseg_cols<-colnames(aseg)[grep("subject|session|timepoint",colnames(aseg),invert=TRUE)] %>% make.names()

cols2check<-c(thal_cols,aseg_cols,paste0("adj_",thal_cols), paste0("ICV_adj_",thal_cols))
cols2check %in% colnames(data)

descriptive_df <- data.frame()
for (vol in cols2check){
  descriptive_stats <- data %>%
    group_by(timepoint) %>%
    # using sym to access the content of variable vol
    summarise(mean = mean(!!sym(vol),na.rm=TRUE), sd = sd(!!sym(vol),na.rm=TRUE)) %>%
    mutate(mean2sd_lower = mean - 2 * sd,
           mean2sd_upper = mean + 2 * sd,
           mean4sd_lower = mean - 4 * sd,
           mean4sd_upper = mean + 4 * sd,
           measure = vol) # add a new column with the volume name
  descriptive_stats <- select(descriptive_stats, timepoint, measure, everything()) # reorder the columns with the new column (nucleus) in the second position
  # print(descriptive_stats)
  descriptive_df <- rbind(descriptive_df, descriptive_stats)
  rm(descriptive_stats)
}

# add column to define the atlas
descriptive_df <- descriptive_df %>%
  mutate(atlas= if_else(measure %in% c(thal_cols, paste0("adj_",thal_cols), paste0("ICV_adj_",thal_cols)),
                        "thalamicNuclei", "aseg" )
  )
# add a column to flag nuclei of interest
## which are...
nuclei<-c("LGN","MGN","Pu","MD","AV")
grep("LGN|MGN|Pu|MD|AV",descriptive_df$measure)
paste(nuclei,collapse = "|")
w<-grep(paste(nuclei,collapse = "|"),descriptive_df$measure) # create an index of all the rows for nuclei of interest
descriptive_df$targetNucleus<-0
descriptive_df$targetNucleus[w]<-1
print(descriptive_df)
print(subset(descriptive_df,targetNucleus==1))


## create new columns such as hemi, nucleus, adjustement...

descriptive_normalized_complete<-subset(descriptive_df,atlas=="thalamicNuclei"&targetNucleus==1) # raw, adjusted_thalamus, adjusted_ICV

# add a column including normalization method (raw, ThalVolume_norm, ICV_norm)

descriptive_normalized_complete

normalization_method<-c(rep("raw",44), rep("ThalVolume_norm ",44), rep("ICV_norm",44)) #rep() = repeat the strings x36 times
descriptive_normalized_complete<-cbind(normalization_method,descriptive_normalized_complete)

descriptive_normalized_complete <- descriptive_normalized_complete %>%
  select("timepoint", "measure", "normalization_method", "mean", "sd", "mean2sd_lower", "mean2sd_upper", "mean4sd_lower", "mean4sd_upper", "atlas", "targetNucleus")


# write_xlsx(descriptive_normalized_complete, 'D:/BCBL/Thalamus/Output/summary_tables/descriptive_statistics_normalization.xlsx')
write_xlsx(descriptive_normalized_complete, 'summary_tables/descriptive_statistics_normalization.xlsx')

#--------------------------------------
# define all variables for which we want adjusted descriptive statistics

dv_cols<-colnames(data)[grep("nihtbx",colnames(data))]

# covariates
covs<- c("timepoint","sex","age",
         # "high.educ.bl","household.income.bl", 
         # "mri_info_device.serial.number" , 
         # "rel_family_id", 
         "abcd_site")

vars2check<-c(paste0("adj_", thal_cols), paste0("ICV_adj_", thal_cols), dv_cols, covs)
# define which ones are numeric and which are factor --> they will be treated differently when computing summary stats
n_vars<-vars2check[grep("numeric", data %>% dplyr::select(all_of(vars2check)) %>% sapply(.,class) )]
c_vars<-vars2check[grep("character|factor", data %>% dplyr::select(all_of(vars2check)) %>% sapply(.,class) )]
#--------------------------------------

# Descriptive Statistics of adjusted thalamic nuclei
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
write.csv(descriptives,file="summary_tables/AdjustedthalamicNuclei_descriptives.csv")
write.csv(descriptives2,file="summary_tables/AdjustedthalamicNuclei_descriptives2.csv",row.names=FALSE)
# openxlsx::write.xlsx(descriptives,file="summary_tables/AdjustedthalamicNuclei_descriptives.xlsx")



## Histograms
if (dir.exists("figures/histograms/adjustedNuclei")==FALSE) { dir.create("figures/histograms/adjustedNuclei") }

# in order to access the content of descriptive_df, you can now use "volume" for filtering.
for (vol in c(paste0("adj_",thal_cols), paste0("ICV_adj_"))) { 
  plot_name<-paste("dist",vol,"histogram.png",sep="_")
  if (!file.exists(paste0("figures/histograms/adjustedNuclei/", plot_name))){ 
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
      geom_text(aes(x = Inf, y = 100, label = paste0("Kurtosis= ", round(kurtosis(data[[vol]], na.rm = TRUE),2))),
                hjust = 1, size = 3, color = "black") #I wished I could put the labels at the top, but the "Inf - 100" does not work.
    
    # you can also save the plots as independent objects
    # save as figure
    ggsave(p, file = paste0("figures/histograms/adjustedNuclei/", plot_name),
           width=4,height=4) # to edit/adjust
    rm(p)
  }

    # remove intermediate objects
  rm(plot_name)
}
#-------------------------------------------------
