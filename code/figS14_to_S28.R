# Purpose: Barplot showing the number and percentages of resistant/sensitive residues
# Author:  Li Li
# Date:    Oct 27, 2025

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")

library(gridExtra)
library(ggpubr)
library(tidyverse)
library(plyr)


dataWk80 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9c.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) 

#week 80 - week 104
dataPostWk80 <- read.csv(file.path(datDir, "d_wk80_wk104_survival_dataset_sieve.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) 



dataWk80_cases <- filter(dataWk80, hiv1event == 1)
dataPostWk80_cases <- filter(dataPostWk80, hiv1event == 1)

escapeDesc <- read.csv(file.path(datDir, "VRC01escape.DescFile.csv"))
escapeDesc <- filter(escapeDesc, pub_id%in% c(dataWk80_cases$pub_id, dataPostWk80_cases$pub_id))
escapeDesc$cohort <- ifelse(escapeDesc$pub_id %in% dataWk80_cases$pub_id, "Weeks 0 to 80", "Weeks 80 to 104")

features <- colnames(escapeDesc)[grepl("score", colnames(escapeDesc))]
features_desc <- c("60","197","234","276","278","279","280","281","317","364","365","369","371","372","425","456",         
                   "459","460","461","463","465","471","474","616","618",
                   "Length of V1V2","No. of PNGs in V1V2","Length of V5","No. of PNGs in V5")
names(features_desc) <- features
tabledf <- tibble("trial" = character(),"cohort" = character(), "tx" = character(),   "feature" = character(), "feature_d" = character(),
                  "num_total" = numeric(), "num_sens" = numeric(), "num_resistant" = numeric(), "num_neutral" = numeric())

plotdf <- tibble("trial" = character(),"cohort" = character(), "tx" = character(),   "feature" = character(), "feature_d" = character(),
                 "id" = numeric(),"value" = numeric())
#calculate the number/percentage of endpoints with resistant/neutral/sensitive features
for(trial in c("703", "704", "703and704")){
  if(trial == "703"){
    df0 <- filter(escapeDesc, Protocol == "V703")
  }else if (trial == "704"){
    df0 <- filter(escapeDesc, Protocol == "V704")
  }else{
    df0 <- escapeDesc
  }
  for(c in c("Weeks 0 to 80", "Weeks 80 to 104")){
    for(x in c("C3", "T1", "T2")){
      df1 <- filter(df0, Tx == x & cohort == c)
      
      for(f in features){
        dfsummary <- table(df1[, f])
        tabledf <- add_row(.data = tabledf, "trial" = trial, "cohort" = c, "tx" = ifelse(x == "C3", "Placebo", 
                                                                                         ifelse(x == "T1", "VRC01, low dose", "VRC01, high dose")), "feature" = f, "feature_d" = features_desc[f],
                           "num_total" = length(df1$pub_id), "num_sens" = ifelse(is.na(dfsummary["-1"]), 0,dfsummary["-1"]) , 
                           "num_resistant" = ifelse(is.na(dfsummary["1"]), 0, dfsummary["1"]), 
                           "num_neutral" = ifelse(is.na(dfsummary["0"]), 0, dfsummary["0"]))
        
        n_sens <- ifelse(is.na(dfsummary["-1"]), 0,dfsummary["-1"])
        n_res <- ifelse(is.na(dfsummary["1"]), 0, dfsummary["1"])
        n_neu <- ifelse(is.na(dfsummary["0"]), 0, dfsummary["0"])
        if(n_res>0){
          plotdf <- add_row(.data = plotdf, "trial" = trial, "cohort" = c, "tx" = ifelse(x == "C3", "Placebo", 
                                                                                         ifelse(x == "T1", "VRC01, low dose", "VRC01, high dose")), "feature" = f, "feature_d" = features_desc[f],
                            "id" = seq(1, n_res, 1), "value" = 1)
        }
        if(n_neu >0){
          plotdf <- add_row(.data = plotdf, "trial" = trial, "cohort" = c, "tx" = ifelse(x == "C3", "Placebo", 
                                                                                         ifelse(x == "T1","VRC01, low dose", "VRC01, high dose")), "feature" = f, "feature_d" = features_desc[f],
                            "id" = seq((n_res+1), n_res+n_neu, 1), "value" = 0)
        }
        if(n_sens >0){
          plotdf <- add_row(.data = plotdf, "trial" = trial, "cohort" = c, "tx" = ifelse(x == "C3", "Placebo", 
                                                                                         ifelse(x == "T1", "VRC01, low dose", "VRC01, high dose")), "feature" = f, "feature_d" = features_desc[f],
                            "id" = seq((n_res+n_neu+1), n_res+n_neu+n_sens, 1), "value" = -1)
        }
        
        
      }
      
    }
  }
  
  
}

for(trialx in c("703", "704", "703and704")){
  plotdf_trialx<- filter(plotdf, trial == trialx)
  plotdf_trialx$z <- 1
  plotdf_trialx$group <- paste(plotdf_trialx$cohort, "\n",plotdf_trialx$tx, sep = "")
  tabledf_703and704 <- filter(tabledf, trial == trialx)
  if(trialx == "703" ){
    plotdf_trialx$width_value <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", 0.8, 0.96)
    plotdf_trialx$ymin1 <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id-0.35, plotdf_trialx$id-0.45)
    plotdf_trialx$ymax1<- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id+0.35, plotdf_trialx$id+0.45)
    
    
  }else if (trialx == "703and704"){
    plotdf_trialx$width_value <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", 0.8, 0.92)
    
    plotdf_trialx$ymin1 <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id-0.2, plotdf_trialx$id-0.43)
    plotdf_trialx$ymax1<- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id+0.2, plotdf_trialx$id+0.43)
  }else{
    plotdf_trialx$width_value <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", 0.8, 0.92)
    plotdf_trialx$ymin1 <- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id-0.35, plotdf_trialx$id-0.45)
    plotdf_trialx$ymax1<- ifelse(plotdf_trialx$cohort == "Weeks 0 to 80", plotdf_trialx$id+0.35, plotdf_trialx$id+0.45)
    
  }
  plotdf_trialx <- left_join(plotdf_trialx, tabledf_703and704, by = c("trial", "cohort", "tx", "feature", "feature_d"))
  
  
  plotdf_trialx$xmin1 <- 0
  plotdf_trialx$xmax1<- 1
  
  #ordering the features by Week 0-80 Placebo frequency of Sensitive
  features_g1 <- c("60","197","234","276","278","279","280","281","317","364","365","369","371","372","425","456",         
                   "459","460","461","463","465","471","474","616","618")
  
  for(plot_groupi in 1:4){
    if(plot_groupi == 1){
      tempdf <- filter(plotdf_trialx, feature_d %in% features_g1[1:7])
      tempdf$feature_d <- factor(tempdf$feature_d, levels = features_g1)
      
    }else if (plot_groupi == 2){
      tempdf <- filter(plotdf_trialx, feature_d %in% features_g1[8:13])
      tempdf$feature_d <- factor(tempdf$feature_d, levels = features_g1)
      
    }else if (plot_groupi == 3){
      tempdf <- filter(plotdf_trialx, feature_d %in% features_g1[14:19])
      tempdf$feature_d <- factor(tempdf$feature_d, levels = features_g1)
      
    }else{
      tempdf <- filter(plotdf_trialx, feature_d %in% features_g1[20:25])
      tempdf$feature_d <- factor(tempdf$feature_d, levels = features_g1)
      
    }
    
    p1 <- list()
    i = 1
    for(g in c("Weeks 0 to 80\nPlacebo", "Weeks 0 to 80\nVRC01, low dose",  "Weeks 0 to 80\nVRC01, high dose",
               "Weeks 80 to 104\nPlacebo", "Weeks 80 to 104\nVRC01, low dose",  "Weeks 80 to 104\nVRC01, high dose")){
      p1[[i]] <- ggplot(data = filter(tempdf, group == g))+
        geom_col(aes(x = id, y = z, fill = as.character(value)), width = 1, color = "white", linewidth = 0.1)+
        scale_y_continuous(breaks = NULL, minor_breaks = NULL, expand = expansion(0))+
        scale_x_discrete(labels = NULL, expand = expansion(0))+
        scale_fill_manual(breaks = c("1", "0", "-1"), values = c("#D7191C", "gray50","#2C7BB6"),
                          labels = c("Resistant", "Neutral", "Sensitive"))+
        coord_flip()+
        facet_grid(rows = vars(feature_d), scales = "free", 
                   switch = "y")+
        ggtitle(g)+
        ylab("") +
        xlab("")
      
      if(i==1){
        p1[[i]] <-  p1[[i]] + 
          theme(axis.text.x = element_blank(),
                plot.margin = margin(t = 0.1, r = 0.05, b = 0.1, l = 0.05, unit = "cm"),
                legend.position = "top",
                legend.direction = "horizontal",
                legend.title = element_blank(),
                panel.spacing = unit(0.05, "lines"),
                plot.title = element_text (size = 10, hjust = 0.5),
                axis.title.x =  element_text (size = 22),
                axis.title.y =  element_text (size = 22, vjust = 0.5),
                strip.text.x = element_text (size = 12, vjust = 0.5),
                strip.text.y = element_text (size = 12, vjust = 0.5)
          )
      }else{
        p1[[i]] <-  p1[[i]] + 
          theme(axis.text.x = element_blank(),
                plot.margin = margin(t = 0.1, r = 0.05, b = 0.1, l = -0.8, unit = "cm"),
                legend.position = "bottom",
                legend.direction = "horizontal",
                panel.spacing = unit(0.05, "lines"),
                plot.title = element_text (size = 10, hjust = 0.5),
                axis.title.x =  element_text (size = 22),
                axis.title.y =  element_text (size = 22, vjust = 0.5),
                strip.text.x = element_blank(),
                strip.text.y = element_blank ()
          )
      }
      
      i = i+1
    }
    
    multiPanelPlot1 <- ggpubr::ggarrange(plotlist = p1[1:6],widths = c(1.55, 1, 1, 1, 1, 1),
                                         nrow = 1, align='h',
                                         common.legend = TRUE,
                                         legend = "bottom")
    
    
    multiPanelPlot1 <- annotate_figure(multiPanelPlot1,,
                                       left=text_grob("Feature", rot=90, size=16))
    
    
    ggsave(file.path(figDir, paste0("trial",trialx, "AAresidueplot_", plot_groupi, ".pdf")),
           plot=multiPanelPlot1, width=8, height=11)
  }
  

  p3 <- list()
  i = 1
  for(g in c("Weeks 0 to 80\nPlacebo", "Weeks 0 to 80\nVRC01, low dose",  "Weeks 0 to 80\nVRC01, high dose",
             "Weeks 80 to 104\nPlacebo", "Weeks 80 to 104\nVRC01, low dose",  "Weeks 80 to 104\nVRC01, high dose")){
    p3[[i]] <- ggplot(data = filter(plotdf_trialx , !feature_d %in% features_g1 & group == g))+
      scale_x_discrete(labels = NULL, expand = expansion(0))+
      geom_col(aes(x = id, y = z, fill = as.character(value)), width = 1, color = "white", linewidth = 0.1)+
      scale_y_continuous(breaks = NULL, minor_breaks = NULL, expand = expansion(0))+
      scale_fill_manual(breaks = c("1", "0", "-1"), values = c("#D7191C", "gray50","#2C7BB6"),
                        labels = c("Resistant", "Neutral", "Sensitive"))+
      coord_flip()+
      facet_grid(rows = vars(feature_d), scales = "free", 
                 switch = "y")+
      ggtitle(g)+
      ylab("") +
      xlab("")
    
    if(i==1){
      p3[[i]] <-  p3[[i]] + 
        theme(axis.text.x = element_blank(),
              plot.margin = margin(t = 0.1, r = 0.05, b = 0.1, l = 0.05, unit = "cm"),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.spacing = unit(0.05, "lines"),
              plot.title = element_text (size = 10, hjust = 0.5),
              axis.title.x =  element_text (size = 22),
              axis.title.y =  element_text (size = 22, vjust = 0.5),
              strip.text.x = element_text (size = 12, vjust = 0.5),
              strip.text.y = element_text (size = 12, vjust = 0.5)
        )
    }else{
      p3[[i]] <-  p3[[i]] + 
        theme(axis.text.x = element_blank(),
              plot.margin = margin(t = 0.1, r = 0.05, b = 0.1, l = -0.8, unit = "cm"),
              legend.position = "bottom",
              legend.direction = "horizontal",
              panel.spacing = unit(0.05, "lines"),
              plot.title = element_text (size = 10, hjust = 0.5),
              axis.title.x =  element_text (size = 22),
              axis.title.y =  element_text (size = 22, vjust = 0.5),
              strip.text.x = element_blank(),
              strip.text.y = element_blank ()
        )
    }
    
    i = i+1
  }
  
  multiPanelPlot2 <- ggpubr::ggarrange(plotlist = p3[1:6],widths = c(1.55, 1, 1, 1, 1, 1),
                                       nrow = 1, align='h',
                                       common.legend = TRUE,
                                       legend = "bottom")
  
  
  multiPanelPlot2 <- annotate_figure(multiPanelPlot2,,
                                     left=text_grob("Feature", rot=90, size=16))
  
  
  ggsave(file.path(figDir, paste0("trial",trialx, "AAfeatureplot.pdf")),
         plot=multiPanelPlot2, width=8, height=8)
  
}