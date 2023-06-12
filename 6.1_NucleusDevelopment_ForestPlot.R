rm(list=ls())
#--------------------------------------
# load libraries

# reading/writing data
library(readr) 
library(writexl)

# tidyverse
library(dplyr)
library(tidyr)
library(tidyverse)

# plotting
library(ggplot2)
library(cowplot)
library(ggpubr)

# regression
library(lme4) # mixed models
library(lmerTest) # pvalues for factors in mixed models, estimating dfs...
library(car) # to get vif
library(MuMIn) # for r.squaredGLMM

# visualization of models
library(sjPlot)
library(ggeffects)
library(glmmTMB)
#---------------------
# define directories and create if needed
# setwd("G:/Amepa/ABCDA/phenotypes/data/working_data/imaging/fs_output/")
setwd("D:/BCBL/Thalamus/Output/")
if (dir.exists("models")==FALSE) { dir.create("models") }
out_dir=paste0(getwd(),"/models/")
if (dir.exists("models/tables")==FALSE) { dir.create("models/tables") }
#--------------------------------------
# read clean data:
data<-readRDS("dataz_AMEPA_cleaned.rds")

#read estimates tables
##baseline model output
estimates0_table <- read_csv("models/tables/estimates0_table.csv")

##full model output
estimates_table <- read_csv("models/tables/estimates_table.csv")

#--------------------------------------
# forest plot for intercept
#--------------------------------------
#get the 95% CI (upper and lower CI) of estimates at the baseline
estimates0_table<- estimates0_table %>%
  mutate(lCI95 = Estimate - qt(0.975, df) * Std..Error, #97.5th percentile of the t-distribution with df degrees of freedom
         uCI95 = Estimate + qt(0.975, df) * Std..Error) #it gives the t-value for a two-tailed test with a 95% confidence level

# estimate0_df<-estimates0_df %>% mutate(iv=gsub("Left.|Right.","",.))  # either whole_thalamus, (intercept)

estimates0_table <- separate(estimates0_table, dv, into = c("Hemisphere", "Nuclei of Interest"), sep = "\\.")
estimates0_table <- separate(estimates0_table, iv, into = c("hemisphere", "region"), sep = "\\.")
estimates0_table$region[is.na(estimates0_table$region)]<-"Intercept"


p0<-estimates0_table %>% 
  ggplot(data=., aes(y=`Nuclei of Interest`,x=Estimate, color=Hemisphere)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh(aes(xmin = lCI95, xmax = uCI95), height = 0.25) +
  facet_grid(.~region,scales="free_x") +
  # scale_alpha_manual(values=c(0.2,1)) +
  theme(legend.position="bottom") +
  panel_border() +
  NULL

ggsave(p0, file = "PlotModels/baselineModel.png", width = 10, height = 9)

#get the 95% CI (upper and lower CI) of estimates at the full model
estimates_table <- estimates_table %>%
  mutate(lCI95 = Estimate - qt(0.975, df) * Std..Error, #97.5th percentile of the t-distribution with df degrees of freedom
         uCI95 = Estimate + qt(0.975, df) * Std..Error) 

estimates_table <- separate(estimates_table, dv, into = c("Hemisphere", "Nuclei of Interest"), sep = "\\.")
estimates_table <- separate(estimates_table, iv, into = c("hemi", "region"), sep = "\\.")
estimates_table$region[is.na(estimates_table$region)]<-"Intercept"


estimates_table <- estimates_table  %>% mutate(hemi=gsub("Left|Right","Whole.thalamus",hemi))

estimates_table$Pfdr<-p.adjust(estimates_table$Pr...t..,method = "fdr")
estimates_table$Pbonf<-p.adjust(estimates_table$Pr...t..,method = "bonferroni")

# to visualize only significant
estimates_table <- estimates_table %>% mutate(sig=if_else(Pbonf<0.05,"*","n.s."))

write.csv(estimates_table, file=paste0(out_dir,"tables/estimates_table_bonferroni.csv"), row.names = FALSE)


# reorder the nuclei: whole_thal, first-order and higher-order relay thal_nuclei
nucleus_order <- c("PuM","PuL","PuI","PuA","totalPul",
                   "MDm","MDl","totalMD",
                   "AV",
                   "MGN","LGN",
                   "Whole_thalamus")

# Reorder the "nucleus" variable 
estimates_table$nucleus <- factor(estimates_table$`Nuclei of Interest`, levels = nucleus_order)

#rename column useful to the plot labels
# Rename levels for the "hemi" variable
estimates_table$hemi <- factor(estimates_table$hemi,
                               levels = c("(Intercept)", "sexM", "age", "age:sexM", "Whole.thalamus"),
                               labels = c("Intercept", "Sex M", "Age", "Age*Sex M", "Total Thalamus"))

# Define the order of facets
facet_order <- c("Intercept", "Age", "Sex M", "Age*Sex M", "Total Thalamus")

# order the "hemi"
estimates_table$hemi <- factor(estimates_table$hemi, levels = facet_order)

# edit the names of the y-label text
estimates_table$`Nuclei of Interest` <- factor(estimates_table$`Nuclei of Interest`,
                               levels = c("PuM","PuL","PuI","PuA","totalPul","MDm","MDl","totalMD","AV","MGN","LGN","Whole_thalamus"),
                               labels = c("PuM","PuL","PuI","PuA","Total Pulvinar","MDm","MDl","Total MD","AV","MGN","LGN","Total Thalamus"))


p1 <- estimates_table %>% 
  ggplot(aes(y = `Nuclei of Interest`, x = Estimate, color = Hemisphere, alpha = sig, shape = sig)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh(aes(xmin = lCI95, xmax = uCI95), height = 0.25) +
  facet_grid(. ~ hemi, scales = "free_x") +
  scale_alpha_manual(values = c(1, 0.5)) +
  labs(x = "Estimate Values") + 
  labs(title = "Developmental Trajectories of the Thalamic Nuclei") +
  theme(legend.position = "bottom",
    legend.text = element_text(family = "Arial", size = 10, face = "plain"),  
    legend.title = element_text(family = "Arial",size = 11, face = "plain"), 
    axis.title.x = element_text(family = "Arial",face = "plain", size = 12),  
    axis.text.x = element_text(family = "Arial",size = 10),  
    axis.title.y = element_text(family = "Arial",face = "plain", size = 12), 
    axis.text.y = element_text(family = "Arial",size = 10), 
    strip.text = element_text(family = "Arial",face = "plain", size = 11),
    plot.title = element_text(family = "Arial",face = "bold", size = 14)  
  ) +
  panel_border()

ggsave(p1, file = "PlotModels/fullModel.png", width = 10, height = 8)
# plot_grid(p0,p1, nrow=2)


# check some specific ones...

## MDl
## where we see:
# * whole thalamus right + effect
# * negtive effect of age, sexM, age*sexM
vol="Right.MDl"
plot_tmp<-ggplot(data=data,aes_string(x="age",y=vol,color="sex")) +
  geom_point() +
  geom_smooth(method='lm',formula= y~x)




