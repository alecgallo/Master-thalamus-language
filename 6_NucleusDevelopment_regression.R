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
# define random factors
random<-colnames(data)[grep("family|subject$|device",colnames(data))]
random_part<-paste(paste0("(1|",random,")"),collapse=" + ")
# define independent variables
ivs<-c("age","sex")

#--------------------------------------
# TO RUN for each dv
# still work in progress but something to start from
#--------------------------------------
for (dv in c(thal_cols_of_interest,thal_cols2adjust)){ 
  cat('\n------------------------------------\n')
  cat("dv=", dv,"\nivs=",ivs,"\nrandom=",random_part)
  cat('\n------------------------------------\n')
  
  # define whole thal to control for
  hemi=strsplit(dv,"\\.") %>% sapply("[[",1) # extract hemisphere from name
  global=paste0(hemi,".Whole_thalamus")
  
  # define if it's a global thalamic measure or not, the baseline model will be different
  if(dv!=global) {
    model_b<-paste(dv," ~ ", random_part, "  + ", global)
  } else {
    model_b<-paste(dv," ~ ", random_part)
  }
  #--------------------------------------
  # baseline model
  #--------------------------------------
  m0<-do.call("lmer",list(as.formula(model_b),data=data))
  #--------------------------------------
  # full model and final model
  #--------------------------------------
  model_full<-paste0(model_b, "+ ", paste(ivs, collapse="*"))
  m1<-do.call("lmer",list(as.formula(model_full),data=data))
  
  # extract estimates from m1 
  table_m1 <- summary(m1) %>% coef() %>% data.frame()
  
  # use step to select model
  step_m<-step(m1,reduce.random=FALSE)
  step_m %>% print()
  
  m_final<-get_model(step_m) # get final model after step-wise regression
  # anova(m0,m_final)
  #--------------------------------------
  # check final model
  anova(m_final)
  summary(m_final)
  
  # check model assumptions, diagnostic plots
  # # assumptions
  m_plot_diag<-plot_model(m_final,type = "diag")
  m_plot_diag[5]<-m_plot_diag[[2]][1]
  m_plot_diag[6]<-m_plot_diag[[2]][2]
  m_plot_diag[7]<-m_plot_diag[[2]][3]
  cowplot::plot_grid(plotlist=m_plot_diag[c(1,3:7)],nrow=2) %>%
    ggsave(.,
           file=paste(out_dir,
                      "/model_diagnosis/",dv,"_m_final","_sjplot",".png",sep=""),
           width=15,height=6
    )
  #--------------------------------------
  # plot main effects
  ## these take some time to plot... so could comment when testing
  plot_model(m_final,type="std") %>% 
    ggsave(.,
           file=paste(out_dir,
                      "/",dv,"_m_final","_std_sjplot",".png",sep=""))
  #
  plot_model(m1,terms=c("age","sex"),type="pred") %>%
    ggsave(.,
           file=paste(out_dir,
                      "/",dv,"_m_final_int","_std_sjplot",".png",sep=""))
  
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
  write.csv(anova_table,file=paste0(out_dir,"/tables/",dv,"_final_","anova_table.csv"),row.names = TRUE)
  write.csv(effect_table,file=paste0(out_dir,"/tables/",dv,"_final_","effect_r2m_table.csv"),row.names = FALSE)
  write.csv(table_m1, file=paste0(out_dir,"/tables/",dv,"_m1_","summary_table.csv"), row.names = TRUE)
  
  # combine all estimates
  sum_table<-sum_table %>% mutate(dv=dv,iv=rownames(.))
  sum_table0<-sum_table0 %>% mutate(dv=dv,iv=rownames(.))
  if(exists("estimates_df")){
    estimates_df<-merge(estimates_df,sum_table,all=TRUE)
    estimates0_df<-merge(estimates0_df,sum_table0,all=TRUE)
    
  } else {
    estimates_df<-sum_table
    estimates0_df<-sum_table0
  }
  write.csv(estimates_df, file=paste0(out_dir,"/tables/estimates_table.csv"), row.names = FALSE)
  write.csv(estimates0_df, file=paste0(out_dir,"/tables/estimates0_table.csv"), row.names = FALSE)
  #--------------------------------------
  # clean intermediate
  rm(m0,m1,step_m,m_final)
  # rm(m_plot_diag)
  #--------------------------------------
  
}


# p1<-plot_model(m1,terms=c("age","sex"),type="pred", 
#                xlab = "Age", 
#                ylab = "Left Total MD (volume)",
#                title="Left Total MD's Developmental Trajectory") #left.totalMD
# 
# ggsave(.,
#        file=paste(out_dir,
#                   "/pulvinars/",dv,"_m_final_int","_std_sjplot",".png",sep=""))
# 
# 
# p2<-plot_model(m1,terms=c("age","sex"),type="pred") #right.totalMD
# p3<-plot_model(m1,terms=c("age","sex"),type="pred") #left.MDl
# p4<-plot_model(m1,terms=c("age","sex"),type="pred") #right.MDl
# p5<-plot_model(m1,terms=c("age","sex"),type="pred") #left.MDm
# p6<-plot_model(m1,terms=c("age","sex"),type="pred") #right.MDm
# 
# 
# p<-plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
# 
# 
# p0<-plot_model(m1,terms=c("age","sex"),type="pred",
#                xlab = "Age", 
#                ylab = "Left Total MD (volume)",
#                title="Left Total MD's Developmental Trajectory") #left.totalMD)
# 
# text(x = min(p0$model$age), y = max(p0$model$left.totalMD),
#      labels = "A", pos = 2, cex = 1.5)
# 
# 
# 
# p1 <- plot_model(m1, terms = c("age", "sex"), type = "pred")
# p1 <- p1 + xlab("Age") + ylab("Right Total Pulvinar (volume)") + 
#   ggtitle("A. Developmental Trajectory of the Right Total Pulvinar") + theme(legend.position = "none")
#           
#           
# p2 <- plot_model(m1, terms = c("age", "sex"), type = "pred")          
# p2 <- p2 + xlab("Age") + ylab("Right Medial Pulvinar (volume)") + 
#   ggtitle("B. Developmental Trajectory of Right Medial Pulvinar")        
#           
# library(gridExtra)
# library(grid)
# 
# combined_plot <- grid.arrange(p1, p2, nrow = 1)
# 
# letter_A <- textGrob("A", x = 0, y = 1, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
# letter_B <- textGrob("B", x = 0, y = 1, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
# 
# combined_plot <- arrangeGrob(
#   combined_plot,
#   top = gTree(children = gList(letter_A)),
#   bottom = gTree(children = gList(letter_B)),
#   widths = unit.c(unit(0.4, "npc"), unit(0.6, "npc")),  # Adjust the width ratio of p1 and p2
#   heights = unit(0.6, "npc")  # Set the height of the combined plot
# )
# grid.newpage()
# grid.draw(combined_plot)



