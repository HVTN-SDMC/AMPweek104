#Purpose: descriptive figures for quantitative marks
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

library(tidyverse)

logit <- function(p){
  return(log(p/(1-p)))
}

dataWk80 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9c.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt80ls
dataWk80$gmt80ls[dataWk80$gmt80ls == ">100"] <- "100"
dataWk80$gmt80ls <- as.numeric(dataWk80$gmt80ls)
#log10 transformation 
dataWk80$gmt80ls <- log10(dataWk80$gmt80ls)




#week 80 - week 104
dataPostWk80 <- read.csv(file.path(datDir, "d_wk80_wk104_survival_dataset_sieve.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt80ls
dataPostWk80$gmt80ls[dataPostWk80$gmt80ls == ">100"] <- "100"
dataPostWk80$gmt80ls <- as.numeric(dataPostWk80$gmt80ls)
#log10 transformation 
dataPostWk80$gmt80ls <- log10(dataPostWk80$gmt80ls)

#week 104
dataWk104 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9_wk104_v2.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt80ls
dataWk104$gmt80ls[dataWk104$gmt80ls == ">100"] <- "100"
dataWk104$gmt80ls <- as.numeric(dataWk104$gmt80ls)
#log10 transformation 
dataWk104$gmt80ls <- log10(dataWk104$gmt80ls)


escapeMark <- read.csv(file.path(datDir, "VRC01escape.DescFile.csv"))
dataWk104 <- left_join(dataWk104, escapeMark, by = "pub_id")
dataWk104$ntx <- ifelse(dataWk104$tx == "C3", "Placebo", ifelse(dataWk104$tx == "T1", "VRC01, low dose", "VRC01, high dose"))

dataWk104$epitope.dist.subtype.ls <- ifelse(dataWk104$protocol == "HVTN 703", dataWk104$epitope.dist.c.ls, dataWk104$epitope.dist.b.ls)
dataWk104$resmark_sens_tbl2 <- 27 - dataWk104$resmark_sens_tbl2
#boxplots and violin plots comparing continuous marks
quantMarks <- c("parscore1.ls","parscore2.ls", "gmt80ls", "epitope.dist.subtype.ls", 
                "hdist.zspace.sites.preselect.all.ls", 
                "hdist.zspace.sites.binding.all.ls", "resmark_sens_tbl2", "resmark_comb_tbl2")

# parscore1.xx: logit predicted probability IC80 >= 1 ug/ml
# parscore2.xx: log base 10 predicted IC80 (ug/ml)
# gmt80ls: IC80ls
# epitope.dist.any.ls: epitope distance to subtype-agnostic reference
# hdist.zspace.sites.preselect.all.ls: PC-weighted HD in the neutralization-associated set
# hdist.zspace.sites.binding.all.ls: PC-weighted HD in the VRC01/CD4 binding set

markFileString <- c("logitPredProbResIC80_ls", "log10predIC80_ls", "IC80_ls", "epitopeDist_ls", 
                    "hdist_ls", "hdist_CD4binding_ls", "resmark_sens_tbl2", "resmark_comb_tbl2")

yLabels <- c(expression("Predicted Probability of" ~ IC[80] > 1 ~ mu * "g/ml"), expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)"),
             expression(IC[80] ~ "(" * mu * "g/ml)"), "VRC01 Epitope Distance to\n Subtype-Specific Reference", 
             "PC-Weighted Hamming Distance in\n27 Positions Predictive of Neutralization",
             "PC-Weighted Hamming Distance in\n50 VRC01 or CD4 Binding Positions",
             "Max(27) - No. of VRC01 Sensitive Features", "VRC01 Resistance Feature Score")



dataWk104_cases <- filter(dataWk104, hiv1event == 1)
dataWk80_cases <- filter(dataWk80, hiv1event == 1)
dataPostWk80_cases <- filter(dataPostWk80, hiv1event == 1)

dataWk104_cases <- filter(dataWk104_cases, pub_id%in% c(dataWk80_cases$pub_id, dataPostWk80_cases$pub_id))
dataWk104_cases$followup <- ifelse(dataWk104_cases$pub_id %in% dataPostWk80_cases$pub_id, "Weeks 80 to 104", "Weeks 0 to 80")
dataWk104_cases$plotg <- paste(dataWk104_cases$followup,"\n", dataWk104_cases$ntx, sep = "")
dataWk104_cases$plotg <- factor(dataWk104_cases$plotg, levels = c("Weeks 0 to 80\nPlacebo",
                         "Weeks 0 to 80\nVRC01, low dose", "Weeks 0 to 80\nVRC01, high dose", 
                         "Weeks 80 to 104\nPlacebo",
                         "Weeks 80 to 104\nVRC01, low dose", "Weeks 80 to 104\nVRC01, high dose"))

for (i in 1:length(quantMarks)){
  if(quantMarks[i] == "gmt80ls" | quantMarks[i] == "parscore2.ls"){
    ytickLab = c(0.07, 0.3, 1, 3, 5, 10, 50, 100)
    ytickAt = log10(ytickLab)
    markRng <- range(dataWk104_cases[, quantMarks[i]], na.rm = TRUE)
    ylim = c(markRng[1] , markRng[2])
    dataWk104_cases$mark <- dataWk104_cases[, quantMarks[i]]
    
  }else if (quantMarks[i] == "parscore1.ls"){
    ytickAt = c(-1, -0.5, 0, 0.5, 1, 2,  3)
    ytickLab = round(exp(ytickAt)/(1+exp(ytickAt)),1)
    markRng <- range(dataWk104_cases[, quantMarks[i]], na.rm = TRUE)
    ylim = exp(c(min(markRng[1]) , markRng[2]))/(1+exp(c(min(markRng[1]) , markRng[2])))
    dataWk104_cases$mark <- dataWk104_cases[, quantMarks[i]]
    
  }else{
    dataWk104_cases$mark <- dataWk104_cases[, quantMarks[i]]
    markRng <- range(dataWk104_cases[, quantMarks[i]], na.rm = TRUE)
    ytickAt = NULL
    ytickLab = NULL
    ylim = NULL
    
  }
  
  for (trial in c("703", "704", "703and704")){
    if(trial == "703"){
      df <- filter(dataWk104_cases, protocol == "HVTN 703")
    }else if (trial == "704"){
      df <- filter(dataWk104_cases, protocol == "HVTN 704")
    }else{
      df <- dataWk104_cases
    }
    set.seed(1)
    p <- ggplot(data = df)+
      geom_violin (aes(x = plotg, y = mark, color = tx))+
      geom_jitter(aes(x = plotg, y = mark, color = tx),width = 0.2,stroke=1.3, fill = "white", shape=21, 
                  size = 4, na.rm = TRUE, height=0, alpha = 0.9)+
      geom_boxplot(aes(x = plotg, y = mark, color = tx), width = 0.15, outlier.shape=NA, fill = NA, fatten = 5)+
      ylab(yLabels[i])+
      xlab("")
    if(!is.null(ylim)){
      p <- p +  scale_y_continuous(limits = markRng, breaks = ytickAt, labels = ytickLab, minor_breaks = NULL)
    }else{
      p <- p +  scale_y_continuous(breaks = scales::breaks_pretty(n = 5))
    }
    p <- p+ 
      scale_color_manual(breaks = c("C3", "T1", "T2"), values = c("blue", "red3","red3")) +
      guides(color = "none") + 
      theme_bw()+
      theme(legend.title = element_blank(),
            plot.margin = unit(c(0.5,0.5,0.5,1), "cm"), #t, r, b, l
            legend.position = c(0.8, 0.95),
            legend.direction = "horizontal",
            legend.text = element_text (size = 16),
            title = element_text(size = 14),
            axis.title.x =  element_text (size = 22),
            axis.title.y =  element_text (size = 22, vjust = 0.5),
            axis.text.x =  element_text (size = 20),
            axis.text.y =  element_text (size = 20)) 
    ggsave(file.path(figDir, paste0(trial,"_desc_", markFileString[i], "_VRC01_placebo.pdf")),
           plot=p, width=16, height=7)
  }
  
}
