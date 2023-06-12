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
estimates0_df <- read_csv("models/tables/thalamic_involvement_language_dev/estimates0_AI.table.csv")

##full model output
estimates_df <- read_csv("models/tables/thalamic_involvement_language_dev/estimates_AI.table.csv")


#--------------------------------------
# forest plot for intercept
#--------------------------------------
#get the 95% CI (upper and lower CI) of estimates at the baseline
estimates0_df <- estimates0_df %>%
  mutate(lCI95 = Estimate - qt(0.975, df) * Std..Error, #97.5th percentile of the t-distribution with df degrees of freedom
         uCI95 = Estimate + qt(0.975, df) * Std..Error) #it gives the t-value for a two-tailed test with a 95% confidence level

# estimate0_df<-estimates0_df %>% mutate(iv=gsub("Left.|Right.","",.))  # either whole_thalamus, (intercept)

# estimates0_df <- separate(estimates0_df, dv, into = c("hemi", "nucleus"), sep = "\\.")



p0<-estimates0_df %>% 
  ggplot(data=., aes(y=dv,x=Estimate)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh(aes(xmin = lCI95, xmax = uCI95), height = 0.25) +
  facet_grid(.~iv,scales="free_x") +
  # scale_alpha_manual(values=c(0.2,1)) +
  theme(legend.position="bottom") +
  panel_border() +
  NULL

ggsave(p0, file = "PlotModels/baselineModel_AI.png", width = 8, height = 7)

#get the 95% CI (upper and lower CI) of estimates at the full model
estimates_df <- estimates_df %>%
  mutate(lCI95 = Estimate - qt(0.975, df) * Std..Error, #97.5th percentile of the t-distribution with df degrees of freedom
         uCI95 = Estimate + qt(0.975, df) * Std..Error) 

# estimates_df <- separate(estimates_df, dv, into = c("hemi", "nucleus"), sep = "\\.")
# estimates_df <- separate(estimates_df, iv, into = c("hemisphere", "region"), sep = "\\.")
# estimates_df$region[is.na(estimates_df$region)]<-"Intercept"


# estimates_df <- estimates_df  %>% mutate(hemisphere=gsub("Left|Right","Whole.thalamus",hemisphere))

estimates_df$Pfdr<-p.adjust(estimates_df$Pr...t..,method = "fdr")
estimates_df$Pbonf<-p.adjust(estimates_df$Pr...t..,method = "bonferroni")

# to visualize only significant
estimates_df <- estimates_df %>% mutate(sig=if_else(Pbonf<0.05,"*","n.s."))

write.csv(estimates_df, file=paste0(out_dir,"tables/thalamic_involvement_language_dev/AI.estimates_table_bonferroni.csv"), row.names = FALSE)

## edit to plot the final results
# change nuclei in iv column with "Thalamic Nuclei"
nuclei <- estimates_df$iv[grep("AI", estimates_df$iv)]
nuclei <- estimates_df$iv[grep("Whole|(Intercept)|age:sexM|age|sexM", estimates_df$iv,invert = TRUE)] %>% unique()

estimates_df$iv %in% nuclei
estimates_df$iv[estimates_df$iv %in% nuclei]
estimates_df$iv[estimates_df$iv %in% nuclei] <- "Thalamic_Nuclei"

estimates_df$iv=="AI.Whole_thalamus" & estimates_df$model=="AI.Whole_thalamus"
estimates_df$iv[estimates_df$iv=="AI.Whole_thalamus" & estimates_df$model=="AI.Whole_thalamus"] <- "Thalamic_Nuclei"
# estimates_df <- separate(estimates_df, model, into = c("Index", "nucleus"), sep = "\\.")

# drop ivs we are not interested in (i.e., Intercept, age, age*sex, sexM)
excluded_ivs <- c("(Intercept)", "age", "sexM", "age:sexM")
estimates_df <- estimates_df[!estimates_df$iv %in% excluded_ivs, ]

# # collapse Left and Right.Whole_thalamus in the plot
# estimates_df <- separate(estimates_df, iv, into = c("hemisphere", "region"), sep = "\\.")
# estimates_df$region[is.na(estimates_df$region)]<-"Thalamic_Nuclei"

# reorder the nuclei: whole_thal, first-order and higher-order relay thal_nuclei
nucleus_order <- c("AI.PuM","AI.PuL","AI.PuI","AI.PuA","AI.totalPul",
                   "AI.MDm","AI.MDl","AI.totalMD",
                   "AI.AV",
                   "AI.MGN","AI.LGN",
                   "AI.Whole_thalamus")

# Reorder the "nucleus" variable 
estimates_df$model <- factor(estimates_df$model, levels = nucleus_order)

# edit names (i.e., reading and vocabulary names to plot "reading" & "vocabulary"; etc)

estimates_df$dv[estimates_df$dv == "nihtbx_picvocab_uncorrected"] <- "Vocabulary Scores"
estimates_df$dv[estimates_df$dv == "nihtbx_reading_uncorrected"] <- "Reading Scores"

estimates_df$iv[estimates_df$iv == "AI.Whole_thalamus"] <- "AI Total Thalamus"
estimates_df$iv[estimates_df$iv == "Thalamic_Nuclei"] <- "Thalamic Nuclei"

# edit the names of the y-label text
estimates_df$model <- factor(estimates_df$model,levels = c("AI.PuM","AI.PuL","AI.PuI","AI.PuA","AI.totalPul","AI.MDm","AI.MDl","AI.totalMD","AI.AV","AI.MGN","AI.LGN","AI.Whole_thalamus"),
                          labels = c("AI PuM","AI PuL","AI PuI","AI PuA","AI Total Pulvinar","AI MDm","AI MDl","AI Total MD","AI AV","AI MGN","AI LGN","AI Total Thalamus"))


# # facet grid order
facet_grid <- c("Thalamic Nuclei", "AI Total Thalamus") #not working though

p1 <- estimates_df %>% filter(iv=="Thalamic Nuclei") %>%
  ggplot(data = ., aes(y = model, x = Estimate, alpha = sig, shape = sig)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh(aes(xmin = lCI95, xmax = uCI95), height = 0.25) +
  facet_grid(.~dv, scales = "free_x") +
  scale_alpha_manual(values = c(1, 0.3)) +
  labs(title="Asymmetry Indexes of the Thalamic Nuclei Involvement in Language Development", 
       x = "Estimate Values",
       y = "Nuclei of Interest") +
  theme(legend.position = "bottom",
        legend.text = element_text(family = "Arial", size = 10, face = "plain"),  
        legend.title = element_text(family = "Arial",size = 11, face = "plain"), 
        axis.title.x = element_text(family = "Arial",face = "plain", size = 12),  
        axis.text.x = element_text(family = "Arial",size = 10),  
        axis.title.y = element_text(family = "Arial",face = "plain", size = 12), 
        axis.text.y = element_text(family = "Arial",size = 10), 
        strip.text = element_text(family = "Arial",face = "plain", size = 11),
        plot.title = element_text(family = "Arial",face = "bold", size = 14)) + 
  panel_border() +
  NULL

ggsave(p1, file = "PlotModels/fullModel/AI.full_model_thalamic_involvement_language_dev.png", width = 10, height = 8)
# plot_grid(p0,p1, nrow=2)


# check some specific ones...

## MDl
## where we see:
# * whole thalamus right + effect
# * negative effect of age, sexM, age*sexM
vol="Right.MDl"
plot_tmp<-ggplot(data=data,aes_string(x="age",y=vol,color="sex")) +
  geom_point() +
  geom_smooth(method='lm',formula= y~x)