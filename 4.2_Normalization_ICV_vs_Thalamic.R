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
library(cowplot)
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

# thal_cols <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv",n_max = 0)[-1] %>% colnames()
thal_cols <- colnames(thal)[-1]
#------------------
################################
## Scatterplots: raw vs adjusted
if (dir.exists("figures/adjustment")==FALSE) { dir.create("figures/adjustment") }


for (vol in thal_cols) {
  if (!file.exists(paste0("figures/adjustment/","comparison_",vol,".png"))){
  
  hemi=strsplit(vol,"\\.") %>% sapply("[[",1)
  global=paste0(hemi,".Whole_thalamus") #paste0() = concatenate
  icv="EstimatedTotalIntraCranialVol"
  
  p1<-data %>%
    ggplot(aes_string(vol,global)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    stat_cor()
  
  # after normalization 
  p2<-data%>%
    ggplot(aes_string(paste0("adj_",vol), global)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm") +
    stat_cor()
  
  p3<-data %>%
    ggplot(aes_string(vol,icv)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    stat_cor()
  
  # after normalization 
  p4<-data%>%
    ggplot(aes_string(paste0("ICV_adj_",vol), icv)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm") +
    stat_cor()
  
  
  library(cowplot)
  combined_plot <- plot_grid(plotlist = list(p1, p2, p3, p4), ncol = 2,
                             labels = c("A", "B", "C", "D"), label_size = 12) # add common title... with nuclei name
  # hints: https://wilkelab.org/cowplot/articles/plot_grid.html
  
  # save it!
  # save
  ggsave(combined_plot,file=paste0("figures/adjustment/","comparison_",vol,".png"),
         width=6,height=5)
  
  # remove intermediate objects
  rm(p1,p2,p3,p4,combined_plot)

} 
}

################################################################################


# check which normalization shows stronger correlations to the thalamic nuclei
correlation_table <- data.frame()

for (vol in thal_cols) {
  hemi <- strsplit(vol, "\\.") %>% sapply("[[", 1)
  global <- paste0(hemi, ".Whole_thalamus")
  icv <- "EstimatedTotalIntraCranialVol"
  
  
  if (vol == global) {
    next  # Skip the current iteration and move to the next vol
  }
  
  # calculate correlation between vol and thalamus 
  cor_vol_thal <- cor.test(data[[vol]], data[[global]])
  cor_adj_vol_thal <- cor.test(data[[paste0("adj_", vol)]], data[[global]])
  
  # calculate correlation between vol and icv
  cor_vol_icv <- cor.test(data[[vol]], data[[icv]])
  cor_adj_vol_icv <- cor.test(data[[paste0("ICV_adj_", vol)]], data[[icv]])

  vol_results <- data.frame(Nuclei = vol,
                            Iv = c("Thal (before normalization)",
                                   "Thal (after normalization)",
                                   "ICV (before normaliation)",
                                   "ICV (after normalization)"),
                            r = c(cor_vol_thal$estimate,
                                            cor_adj_vol_thal$estimate,
                                            cor_vol_icv$estimate,
                                            cor_adj_vol_icv$estimate))
  
  correlation_table <- rbind(correlation_table, vol_results)
}

write_xlsx(correlation_table, 'summary_tables/ICV_vs_Thalamic_normalization_correlations.xlsx')

mean(correlation_table$r[correlation_table$Iv == "Thal (before normalization)"])
mean(correlation_table$r[correlation_table$Iv == "ICV (before normaliation)"])

#-------------------------------------------------
