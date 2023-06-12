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
aseg <- read_csv("allsubjects_aseg_stats_vol_mm3.csv") #BCBL's aseg segmentation (output by Liu, Rodriguez, and Carrion)
thal<-read.csv("thal.csv")
# thal <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv") #BCBL's thalamic nuclei volumes (output by Liu, Rodriguez, and Carrion)


# thal_cols <- read_csv("allsubjects_ThalamicNuclei.T1.volumes.csv",n_max = 0)[-1] %>% colnames()
thal_cols <- colnames(thal)[-1]
# get rid of "-" in colnames with "."
thal_cols<- make.names(thal_cols)
#------------------

#THALAMIC NUCLEI NORMALIZATION
##1.1 [Hemi] thalamus normalization: regression / normalization of the thalamic nuclei (and add the new adjusted nuclei to the "data") using whole thalamic volume

for (vol in thal_cols){
  # define variables
  hemi=strsplit(vol,"\\.") %>% sapply("[[",1)
  global=paste0(hemi,".Whole_thalamus") #paste0() = concatenate
  # names(data) <- make.names(names(data), unique=TRUE)  #I had to add this line because the global_mean line was complaining, that there were duplicated columns without names # now redundant?
  global_mean= data %>% pull(global) %>% mean()
  # define model
  model=paste(vol,"~", global)
  print(model)
  m0<-do.call("lm",list(as.formula(model),data=data))
  # summary(m0)
  b=coefficients(m0)[global] #global refers to the Left/Right Whole Thalamus, depending on the vol selected in the loop
  # generate new column using the formula
  # volAdj=volInd-b(globalInd-globalMean)
  volAdj=data[,vol]-b*(data[,global]-global_mean) #volAdj are the normalized values. 10656 values of each nucleus normalized across all participants
  volAdj<-as.data.frame(volAdj) # create a variable to then being able to add each of the normalized volumes in the df
  colnames(volAdj)<- paste0("adj_",vol) # call the new column, creating a new string ("adj" `+ name of the nucleus`)
  # add new adjusted volume data into dataframe
  data<-cbind(data,volAdj)
  # data$adj_Left_AV

  # remove intermediate objects
  rm(hemi,global,global_mean,model,m0,b,volAdj)
}
rm(vol)

# before normalization (Left_LGN/Left_Whole_thalamus)
data %>%
  ggplot(aes(Left.LGN,Left.Whole_thalamus)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor()


# after normalization
data%>%
  ggplot(aes(adj_Left.LGN,Left.Whole_thalamus)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor()


## 1.2 ICV normalization: regression / normalization of the thalamic nuclei (and add the new adjusted nuclei to the "data") using IntraCranial Volume
global="EstimatedTotalIntraCranialVol" #AG: I changed "EstimatedTotalIntraCranialVol" to "EstimatedTotalIntraCranialVol.BCBL" because in the previous script we used the EstimatedTotalIntraCranialVol which was contained in the aseg segmentation
global_mean = data %>% pull(global) %>% mean()

for (vol in thal_cols){
  # define model
  model=paste(vol,"~", global)
  print(model)
  m0<-do.call("lm",list(as.formula(model),data=data))
  # summary(m0)
  b=coefficients(m0)[global]
  # generate new column using the formula
  # volAdj=volInd-b(globalInd-globalMean)
  # data<-cbind(data,ICV)
  ICV_volAdj=data[,vol]-b*(data[,global]-global_mean)
  ICV_volAdj<-as.data.frame(ICV_volAdj)
  colnames(ICV_volAdj)<- paste0("ICV_adj_",vol)
  # add new adjusted volume data into dataframe
  data<-cbind(data,ICV_volAdj)

  # remove intermediate objects
  rm(model,m0,b,ICV_volAdj)
}
rm(vol)

colnames(data) <- make.names(colnames(data), unique = TRUE)

saveRDS(data, file = "data_AMEPA_Normalized.rds")

