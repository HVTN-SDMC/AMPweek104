# PE by IC50 using the nonparametric method

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")

library(sievePH)
library(tidyverse)

logit <- function(p){
  return(log(p/(1-p)))
}

source(file.path(repoDir, "code/macro/plot.summary.sievePH.R"))
source(file.path(repoDir, "code/macro/plot.summary.sievePH_wk80Wk104.R"))

#week 80 
dataWk80 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9c.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt50ls
dataWk80$gmt50ls[dataWk80$gmt50ls == ">100"] <- "100"
dataWk80$gmt50ls <- as.numeric(dataWk80$gmt50ls)
#log10 transformation 
dataWk80$gmt50ls <- log10(dataWk80$gmt50ls)


#week 80 - week 104
dataPostWk80 <- read.csv(file.path(datDir,  "d_wk80_wk104_survival_dataset_sieve.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt50ls
dataPostWk80$gmt50ls[dataPostWk80$gmt50ls == ">100"] <- "100"
dataPostWk80$gmt50ls <- as.numeric(dataPostWk80$gmt50ls)
#log10 transformation 
dataPostWk80$gmt50ls <- log10(dataPostWk80$gmt50ls)
#PubIDs for the VRC01 participants in figure 6B with IC80 = 100 and the placebo recipient with an IC80 ~7
dataPostWk80$pub_id[!is.na(dataPostWk80$gmt50ls) & dataPostWk80$gmt50ls == log10(100) & dataPostWk80$tx != "C3"]
dataPostWk80$pub_id[!is.na(dataPostWk80$gmt50ls) & dataPostWk80$gmt50ls > log10(5) & dataPostWk80$tx == "C3"]


#week 104
dataWk104 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9_wk104_v2.csv")) %>%
  # two 703 ppts with missing sequences have also a missing time-to-event
  filter(!is.na(hiv1fpday)) %>%
  # stratification variable 
  mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                            protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                            protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                            protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))

#pre-process gmt50ls
dataWk104$gmt50ls[dataWk104$gmt50ls == ">100"] <- "100"
dataWk104$gmt50ls <- as.numeric(dataWk104$gmt50ls)
#log10 transformation 
dataWk104$gmt50ls <- log10(dataWk104$gmt50ls)


mark = "gmt50ls"
xLim <- range(c(dataPostWk80[, mark], dataWk80[, mark]), na.rm=TRUE)
dataPostWk80Sub <- dplyr::select(dataPostWk80, all_of(c("tx_pool", "hiv1fpday", "hiv1event", mark,  "stratVar")))
dataWk80Sub <- dplyr::select(dataWk80, all_of(c("tx_pool", "hiv1fpday", "hiv1event", mark,  "stratVar")))
dataWk104Sub <-  dplyr::select(dataWk104, all_of(c("tx_pool", "hiv1fpday", "hiv1event", mark,  "stratVar")))
colnames(dataPostWk80Sub) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
colnames(dataWk80Sub) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
colnames(dataWk104Sub) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")

dataPostWk80Sub$tx <- ifelse(as.character(dataPostWk80Sub$tx)=="C3", 0, 1)
dataWk80Sub$tx <- ifelse(as.character(dataWk80Sub$tx)=="C3", 0, 1)
dataWk104Sub$tx <- ifelse(as.character(dataWk104Sub$tx)=="C3", 0, 1)

# convert mark values for non-primary endpoints through tau to NA
dataPostWk80Sub$mark <- ifelse(dataPostWk80Sub$eventInd==0, NA, dataPostWk80Sub$mark)
dataWk80Sub$mark <- ifelse(dataWk80Sub$eventInd==0, NA, dataWk80Sub$mark)
dataWk104Sub$mark <- ifelse(dataWk104Sub$eventInd==0, NA, dataWk104Sub$mark)

# complete-case analysis, i.e., discard cases with a missing mark
dataPostWk80Sub <- subset(dataPostWk80Sub, !(eventInd==1 & is.na(mark)))
dataWk80Sub <- subset(dataWk80Sub, !(eventInd==1 & is.na(mark)))
dataWk104Sub <- subset(dataWk104Sub, !(eventInd==1 & is.na(mark)))


dataPostWk80Sub$treatment <- ifelse(dataPostWk80Sub$tx==0, "Placebo", "VRC01")
dataWk80Sub$treatment <- ifelse(dataWk80Sub$tx==0, "Placebo", "VRC01")


#refit week 80, post week 80 week 104 using the nonparametric method
loadresult <- FALSE
defaultbandwidth <- FALSE
if(loadresult){
  if(defaultbandwidth){
    fitWk80 <- readRDS(file.path(outDir, "nonparametricfitWk80_ic50_defaultbw.rds"))
    fitWk104 <- readRDS(file.path(outDir, "nonparametricfitWk104_ic50_defaultbw.rds"))
    fitPostWk80 <- readRDS(file.path(outDir, "nonparametricfitPostWk80_ic50_defaultbw.rds"))
    
  }else{
    fitWk80 <- readRDS(file.path(outDir, "nonparametricfitWk80_ic50_fixedbandwitdh.rds"))
    fitWk104 <- readRDS(file.path(outDir, "nonparametricfitWk104_ic50_fixedbandwitdh.rds"))
    fitPostWk80 <- readRDS(file.path(outDir, "nonparametricfitPostWk80_ic50_fixedbandwitdh.rds"))
    
  }
  
  
}else{
  dataWk80Sub$eventTime2 <- pmax(dataWk80Sub$eventTime, 1) / 365.25
  dataPostWk80Sub$eventTime2 <- pmax(dataPostWk80Sub$eventTime, 1) / 365.25
  dataWk104Sub$eventTime2 <- pmax(dataWk104Sub$eventTime, 1) / 365.25
  
  if(defaultbandwidth){
    hbandic50 <- 9*sqrt(var((dataWk80Sub$mark - min(dataWk80Sub$mark, na.rm = TRUE))/(max(dataWk80Sub$mark, na.rm = TRUE) 
                                                                                      - min(dataWk80Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataWk80Sub$eventInd)^{-1/3}
    fitWk80 <- with(dataWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 85.9/52, tband = 0.5,
                                                hband = hbandic50, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    hbandic50 <- 9*sqrt(var((dataWk104Sub$mark - min(dataWk104Sub$mark, na.rm = TRUE))/(max(dataWk104Sub$mark, na.rm = TRUE) 
                                                                                        - min(dataWk104Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataWk104Sub$eventInd)^{-1/3}
    fitWk104 <- with(dataWk104Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 108.6/52, tband = 0.5,
                                                  hband = hbandic50, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    hbandic50 <- 9*sqrt(var((dataPostWk80Sub$mark - min(dataPostWk80Sub$mark, na.rm = TRUE))/(max(dataPostWk80Sub$mark, na.rm = TRUE) 
                                                                                              - min(dataPostWk80Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataPostWk80Sub$eventInd)^{-1/3}
    fitPostWk80 <- with(dataPostWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 24/52, tband = 0.5,
                                                        hband = hbandic50, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    saveRDS(fitWk80, file.path(outDir, "nonparametricfitWk80_ic50_defaultbw.rds"))
    saveRDS(fitWk104, file.path(outDir, "nonparametricfitWk104_ic50_defaultbw.rds"))
    saveRDS(fitPostWk80, file.path(outDir, "nonparametricfitPostWk80_ic50_defaultbw.rds"))
  }else{
    fitWk80 <- with(dataWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 85.9/52, tband = 0.5,
                                                hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    
    fitWk104 <- with(dataWk104Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 108.6/52, tband = 0.5,
                                                  hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    fitPostWk80 <- with(dataPostWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 24/52, tband = 0.5,
                                                        hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    saveRDS(fitWk80, file.path(outDir, "nonparametricfitWk80_ic50_fixedbandwitdh.rds"))
    saveRDS(fitWk104, file.path(outDir, "nonparametricfitWk104_ic50_fixedbandwitdh.rds"))
    saveRDS(fitPostWk80, file.path(outDir, "nonparametricfitPostWk80_ic50_fixedbandwitdh.rds"))
    
  }
}

sfitWk80 <- summary(fitWk80, sieveAlternative = "oneSided")
sfitWk104 <- summary(fitWk104, sieveAlternative = "oneSided")
sfitPostWk80 <- summary(fitPostWk80, sieveAlternative = "oneSided")


#plot Wk80 and Wk104 together
sfitWk80$cohort <- "Week 80"
sfitWk104$cohort <- "Week 104"
sfit = list(sfitWk80, sfitWk104)
dataWk104Sub$treatment <- ifelse(dataWk104Sub$tx==0, "Placebo", "VRC01")
xtickLab = c(0.07, 0.3, 1, 3, 5, 10, 25, 50, 100)
xtickAt = log10(xtickLab)

p <- ggplotSievePH_wk80wk104(sfit, data1 = dataWk80Sub[, c("treatment", "mark", "tx")], 
                             data2 = dataWk104Sub[, c("treatment", "mark", "tx")],  
                             xlab=expression(IC[50] ~ "(" * mu * "g/ml)"),
                             breaks.x = xtickAt, labels.x = xtickLab, title="")

ggsave(file.path(figDir, paste0("703and704_sievePH_PEbyIC50_ls_VRC01_placebo_Wk80Wk104nonparametric", ifelse(defaultbandwidth, "_dbw",""), ".pdf")),
       plot=p, width=1 * 7, height=0.9 * 8)




pPostWk80 <- ggplotSievePH_onecohort (sfitPostWk80, dataPostWk80Sub[, c("treatment", "mark", "tx")],
                                      xlab=expression(IC[50] ~ "(" * mu * "g/ml)"),
                                      breaks.x = xtickAt, labels.x = xtickLab, title="", color_line = "darkgreen",
                                      height.panels = c(1, 2), cohort = "Weeks 80 to 104", legendLabels = "Weeks 80 to 104",
                                      limits.x = range(dataWk104Sub$mark, na.rm =TRUE))


ggsave(file.path(figDir, paste0("703and704_sievePH_PEbyIC50_ls_VRC01_placebo_postWk80nonparametric", ifelse(defaultbandwidth, "_dbw",""), ".pdf")),
       plot=pPostWk80, width=1 * 7, height=0.9 * 8)


