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
data<-readRDS("data_AMEPA.rds")
thal<-read.csv("thal.csv")

# define columns for thalamic nuclei
# thal_cols <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv",n_max = 0)[-1] %>% colnames()
thal_cols <- colnames(thal)[-1]
# get rid of "-" in colnames with "."
thal_cols<- make.names(thal_cols)
# define aseg cols
aseg_cols <- read_csv("allsubjects_aseg_stats_vol_mm3.csv")[-1] %>% colnames()
aseg_cols <- make.names(aseg_cols)
# define aseg cols from abcd run, which contain "aseg" as keyword
aseg_abcd_cols<-colnames(data)[grep("aseg",colnames(data))]

################################
# QC: comparing aseg volumes from two runs (for the same subject_session)
################################

# Create a dataframe with columns for the ABCD aseg... and edit them so that the names match the BCBL aseg columns
df <- data.frame(
  ABCD_original = aseg_abcd_cols,
  ABCD_matched=aseg_abcd_cols
)

### match names across datasets

# replace lh with Left in correct position
df$ABCD_matched[grep("lh$",df$ABCD_matched)] <- gsub("smri_vol_subcort.aseg_","ABCD.Left.",df$ABCD_matched[grep("lh$",df$ABCD_matched)])  %>% 
  gsub(".lh$","",.)

df$ABCD_matched[grep("rh$",df$ABCD_matched)] <- gsub("smri_vol_subcort.aseg_","ABCD.Right.",df$ABCD_matched[grep("rh$",df$ABCD_matched)]) %>%
  gsub(".rh$","",.)

df$ABCD_matched[grep("^s",df$ABCD_matched)] <- gsub("smri_vol_subcort.aseg_","ABCD.",df$ABCD_matched[grep("^s",df$ABCD_matched)]) %>%
  gsub("^s","",.)

df$ABCD_matched[grep("^A",df$ABCD_matched)] <- gsub("ABCD.", "", df$ABCD_matched[grep("^A",df$ABCD_matched)]) %>% gsub("$",".ABCD",.)

# make every letter after a "." in upper case
df$ABCD_matched <- gsub("(?<=\\.)([a-z])", "\\U\\1", df$ABCD_matched, perl = TRUE)

# change "X3/4/5" for "3/4/5"
df$ABCD_matched[11]<-"X3rd.Ventricle.ABCD"
df$ABCD_matched[12]<-"X4th.Ventricle.ABCD"
df$ABCD_matched[1]<-"EstimatedTotalIntraCranialVol.ABCD"

# other manual changes
df$ABCD_matched[16] <- "CSF.ABCD"
df$ABCD_matched[7] <- "Left.Thalamus.ABCD"
df$ABCD_matched[18] <- "Left.VentralDC.ABCD"
df$ABCD_matched[24] <- "Right.Thalamus.ABCD"
df$ABCD_matched[31] <- "Right.VentralDC.ABCD"
df$ABCD_matched[32] <- "WM.Hypointensities.ABCD"
df$ABCD_matched[42] <- "SubCortGrayVol.ABCD"
df$ABCD_matched[41] <- "SupraTentorialVol.ABCD"
df$ABCD_matched[2] <- "LhCerebralWhiteMatterVol.ABCD"
df$ABCD_matched[19] <- "RhCerebralWhiteMatterVol.ABCD"

# create new column with actual matching name
df$ABCD_matched<-gsub("(?<=\\.)([a-z])", "\\U\\1", df$ABCD_matched, perl = TRUE) %>%
  gsub("(?<=^)([a-z])", "\\U\\1",., perl = TRUE) %>%
  gsub("Cc.","CC_",.) %>%
  gsub(".Anterior","_Anterior",.) %>%
  gsub(".Posterior","_Posterior",.)
# double check that all ABCD_matched names have the suffix ".ABCD"
df$ABCD_matched[grep("ABCD",df$ABCD_matched,invert=TRUE)] # ok, all have it!

df$match<-gsub(".ABCD","",df$ABCD_matched)

# BCBL's aseg
df2 <- data.frame(
  BCBL_original = aseg_cols,
  match = gsub("(?<=\\.)([a-z])", "\\U\\1", aseg_cols, perl = TRUE) %>%
    gsub("(?<=^)([a-z])", "\\U\\1",., perl = TRUE)) %>%
  mutate(BCBL_matched=paste0(match,".BCBL") )

# combine both column names to get key of matching cols
df_aseg_cols<-merge(df,df2,by="match",all=TRUE)

write_xlsx(df_aseg_cols, 'summary_tables/aseg_correspondence.xlsx')

# remove intermediate objects
rm(df,df2)

################################
# change aseg colnames from data so that they are consistent across ABCD and BCBL runs
## ABCD
w1<-match(df_aseg_cols$ABCD_original,colnames(data))
w1<-w1[!is.na(w1)] # select all the non-missing values
colnames(data)[w1]<-df_aseg_cols$ABCD_matched[!is.na(df_aseg_cols$ABCD_original)]

## BCBL
table(df_aseg_cols$BCBL_original %in% colnames(data))
# there are three that are not present, double check those
df_aseg_cols$BCBL_original[!df_aseg_cols$BCBL_original %in% colnames(data)]  # ok, just NA's
# get indices to replace names
w2<-match(df_aseg_cols$BCBL_original,colnames(data))
w2<-w2[!is.na(w2)]
colnames(data)[w2]<-df_aseg_cols$BCBL_matched[!is.na(df_aseg_cols$BCBL_original)]

# remove intermediate objects w1 and w2
rm(w1,w2)

################################
# ASEG volume comparison across BCBL and ABCD runs
################################

# aseg_cols<-colnames(aseg)[grep("subject|session|timepoint",colnames(aseg),invert=TRUE)] %>% make.names()
# vol="Left.Thalamus"

subset(df_aseg_cols,!is.na(ABCD_matched)&!is.na(BCBL_matched))
aseg_cols<-subset(df_aseg_cols,!is.na(ABCD_matched)&!is.na(BCBL_matched))$match

#-------------------------------------------------------
# Get correlations for each volume computed by BCBL or ABCD
df_aseg_cols$r<-NA # creat new col to save correlation values
for(vol in aseg_cols) {
  plot_name<-paste0("plot_comparison_ABCDvsBCBL_",vol)
  if (!file.exists(paste0("correlation_plots/", vol, ".png"))){ 
  
    v1=paste0(vol,".BCBL")
    v2=paste0(vol,".ABCD")
    r<-cor(data[,v1],data[v2]) # this will give you the r, which we can save
    df_aseg_cols$r[df_aseg_cols$match==vol]<-r
    
    # using plot from ggpubr: https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
    # assign to access it later

    if (!file.exists(paste0("correlation_plots/", plot_name, ".png"))){
    p<-ggscatter(data=data,
                 x=v1,y=v2,
                 color="timepoint",
                 # label.x = -2, label.y = 2,
                 palette = "jco",
                 rug=TRUE,
                 add = "reg.line",
                 # label="id",
                 # label.select = d %>% filter(grupo=="clinico") %>% pull(id),
                 conf.int = FALSE, alpha=0.5) +
      stat_cor() + # include the correlation in the plot
      theme(legend.position="bottom")
    
    # assign(plot_name,p) 
    # save plot
    ggsave(p, file = paste0("correlation_plots/", plot_name, ".png"),
           width=4,height=4)
    rm(p)    
  }
  # remove intermediate
  rm(v1,v2,r)
  }
  rm(plot_name)
}

# you can check the correlations in the df_aseg_cols table...
subset(df_aseg_cols,!is.na(r)) %>% 
  write.csv(., file="summary_tables/aseg_comparison_ABCDvsBCBL_correlations.csv",
            row.names=FALSE)

#-------------------------------------------------

