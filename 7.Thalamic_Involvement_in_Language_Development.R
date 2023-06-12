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
if (dir.exists("models/tables/thalamic_involvement_language_dev/")==FALSE) { dir.create("models/tables/thalamic_involvement_language_dev/") }
#--------------------------------------
# read clean data:
data<-readRDS("dataZ_AMEPA_cleaned.rds")
#--------------------------------------
# define variables
# independent variables
thal_cols <- readRDS("thal_AMEPA.rds")[-1] %>% colnames()
# thal_cols<-colnames(data)[grep("Left|Right",colnames(data))]
# should exclude the small nuclei...

# define columns of interest
thal_cols_of_interest<-thal_cols[grep("MD|Pu|LGN|MGN|AV",thal_cols)] # AG: me he dado cuenta que no cogia los "Pulvinar" porque habia "Pul", Ya he cambiado a "Pu"
thal_cols2adjust<-thal_cols[grep("Whole",thal_cols)]
thal_cols_of_interest2 <- c(thal_cols_of_interest,thal_cols2adjust)
# define random factors
random<-colnames(data)[grep("family|subject$|device",colnames(data))]
random_part<-paste(paste0("(1|",random,")"),collapse=" + ")
# define independent variables
ivs<-c("age","sex")
ivs2<-c(thal_cols_of_interest2)

# define dependent variables
dvs <- data[grep("nihtbx_picvocab_uncorrected|nihtbx_reading_uncorrected", names(data))]
dvs2<-colnames(dvs)

# Define independent variables
# ivs <- c("age", "sex", thal_cols_of_interest)  # Add the nuclei as independent variables

# TO RUN for each dv
for (dv in dvs2) {
  
  for (iv in ivs2){
  cat('\n------------------------------------\n')
  cat("dv=", dv, "\nivs=", ivs, "\nivs2=", ivs2, "\nrandom=", random_part)
  cat('\n------------------------------------\n')
  
  # Define whole thal to control for
  hemi <- strsplit(iv, "\\.") %>% sapply("[[", 1)  # Extract hemisphere from name
  global <- paste0(hemi, ".Whole_thalamus")
  
  # Define if it's a global thalamic measure or not, the baseline model will be different
  if (dv != global) {
    model_b <- paste(dv, " ~ ", random_part, "  + ", global)
  } else {
    model_b <- paste(dv, " ~ ", random_part)
  }
  
  # Baseline model
  m0 <- do.call("lmer", list(as.formula(model_b), data = data))
  
  # Full model and final model
  model_full <- paste0(model_b, " + ", paste(ivs, collapse = "*"), " + ", paste(iv, collapse = "+"))
  m1 <- do.call("lmer", list(as.formula(model_full), data = data))
  
  # Extract estimates from m1
  table_m1 <- summary(m1) %>% coef() %>% data.frame()
  
  # Use step to select model
  step_m <- step(m1, reduce.random = FALSE)
  step_m %>% print()
  
  m_final <- get_model(step_m)  # Get final model after step-wise regression
  
  # Check final model
  anova(m_final)
  summary(m_final)
  
  # Check model assumptions, diagnostic plots
  # Assumptions
  m_plot_diag <- plot_model(m_final, type = "diag")
  m_plot_diag[5] <- m_plot_diag[[2]][1]
  m_plot_diag[6] <- m_plot_diag[[2]][2]
  m_plot_diag[7] <- m_plot_diag[[2]][3]
  cowplot::plot_grid(plotlist = m_plot_diag[c(1, 3:7)], nrow = 2) %>%
    ggsave(., file = paste(out_dir, "/model_diagnosis/thalamic_involvement_language_dev/", iv, "_m_final", "_sjplot", ".png", sep = ""),
           width = 15, height = 6)
  
  #--------------------------------------
  # plot main effects
  ## these take some time to plot... so could comment when testing
  plot_model(m_final,type="std") %>% 
    ggsave(.,
           file=paste(out_dir,
                      "/thalamic_involvement_language_dev/",dv, "_", iv,"_m_final","_std_sjplot",".png",sep=""))
  #
  plot_model(m1,terms=c(iv, "age","sex"),type="pred") %>%
    ggsave(.,
           file=paste(out_dir,
                      "/thalamic_involvement_language_dev/",dv, "_", iv,"_m_final_int","_std_sjplot",".png",sep=""))
  
  #--------------------------------------
  # summary table, parameter table - to get t-values
  sum_table<-summary(m1) %>% coef()  %>% data.frame()
  sum_table0<-summary(m0) %>% coef()  %>% data.frame()
  
  # ANOVA table, for all parameters
  anova_table<-step_m$fixed %>% as.data.frame()
  # effect size, r2 for the fixed effects(r2m)
  effect_table<-data.frame(model=formula(m_final) %>% as.character() %>% paste0(.,collapse="") %>% gsub("^~","",.) %>% gsub(dv,paste0(dv," ~ "),.),
                           R2m=r.squaredGLMM(m_final)[1],
                           R2c=(r.squaredGLMM(m_final)[2]))
  # save summary tables
  write.csv(anova_table,file=paste0(out_dir,"tables/thalamic_involvement_language_dev/",dv,"_final_","anova_table.csv"),row.names = TRUE)
  write.csv(effect_table,file=paste0(out_dir,"tables/thalamic_involvement_language_dev/",dv,"_final_","effect_r2m_table.csv"),row.names = FALSE)
  write.csv(table_m1, file=paste0(out_dir,"tables/thalamic_involvement_language_dev/",dv,"_m1_","summary_table.csv"), row.names = TRUE)
  
  
  
  # add a column to indicate which model it corresponds to
  # define nucleus of the model
  model <- iv
  
  # Variables to discard
  variables_to_discard <- c("(Intercept)", "Whole_thalamus", "age", "sexM", "age:sexM")
  # Select only the variables not in variables_to_discard
  
  selected_variable <- subset(model, !model %in% variables_to_discard)
  # estimates0_df <- estimates0_df %>% mutate(Model = selected_variable) #doesn't work, it doesn't assign it as I wish
  # estimates_df <-estimates_df %>% mutate(model = selected_variable)
  
  # combine all estimates
  sum_table<-sum_table %>% mutate(dv=dv,iv=rownames(.))
  sum_table<-sum_table %>% mutate(model=selected_variable)
  sum_table0<-sum_table0 %>% mutate(dv=dv,iv=rownames(.))
  if(exists("estimates_df")){
    estimates_df<-merge(estimates_df,sum_table,all=TRUE)
    estimates0_df<-merge(estimates0_df,sum_table0,all=TRUE)
    
  } else {
    estimates_df<-sum_table
    estimates0_df<-sum_table0
  }
  
  write.csv(estimates_df, file=paste0(out_dir,"tables/thalamic_involvement_language_dev/estimates_table.csv"), row.names = FALSE)
  write.csv(estimates0_df, file=paste0(out_dir,"tables/thalamic_involvement_language_dev/estimates0_table.csv"), row.names = FALSE)
  #--------------------------------------
  # clean intermediate
  rm(m0,m1,step_m,m_final)
  # rm(m_plot_diag)
  #--------------------------------------
  }
}

