# Purpose: Primary endpoint analysis; estimation of hazard ratio-based PE by quantitaive marks 

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

quantMarks <- c("parscore1.ls","parscore2.ls", "gmt80ls", "epitope.dist.subtype.ls", 
          "hdist.zspace.sites.preselect.all.ls", 
          "hdist.zspace.sites.binding.all.ls")
# parscore1.xx: logit predicted probability IC80 >= 1 ug/ml
# parscore2.xx: log base 10 predicted IC80 (ug/ml)
# gmt80ls: IC80ls
# epitope.dist.any.ls: epitope distance to subtype-agnostic reference
# hdist.zspace.sites.preselect.all.ls: PC-weighted HD in the neutralization-associated set
# hdist.zspace.sites.binding.all.ls: PC-weighted HD in the VRC01/CD4 binding set

# tags that are part of the output PDF file names
trialFileString <- "704and703"
doseFileString <- "VRC01pooled"
markFileString <- c("logitPredProbResIC80_ls", "log10predIC80_ls", "IC80_ls", "epitopeDist_ls", "hdist_ls", "hdist_CD4binding_ls")

# tags that are part of the plot titles
doseTitleString <- c("PE")
markType <- rep(c("ProbResIC80", "IC80"), each=3)

# plot labels
VRC01lab <- c("VRC01")
xLabels <- c(expression("Predicted Probability of" ~ IC[80] > 1 ~ mu * "g/ml"), expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)"),
            expression(IC[80] ~ "(" * mu * "g/ml)"), "VRC01 Epitope Distance to\n Subtype-Matched Reference", 
            "PC-Weighted Hamming Distance in\n27 Positions Predictive of Neutralization",
            "PC-Weighted Hamming Distance in\n50 VRC01 or CD4 Binding Positions")


#week 80 
dataWk80 <- read.csv(file.path(datDir, "amp_sieve_marks_wk80.csv")) %>%
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
dataWk80$epitope.dist.subtype.ls <- ifelse(dataWk80$protocol == "HVTN 703", dataWk80$epitope.dist.c.ls, dataWk80$epitope.dist.b.ls)




#week 80 - week 104
dataPostWk80 <- read.csv(file.path(datDir, "amp_sieve_marks_wk80to104.csv")) %>%
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
dataPostWk80$gmt80ls <- log10(dataPostWk80$gmt80ls)
dataPostWk80$pub_id[!is.na(dataPostWk80$gmt80ls) & dataPostWk80$gmt80ls == log10(100) & dataPostWk80$tx != "C3"]
dataPostWk80$pub_id[!is.na(dataPostWk80$gmt80ls) & dataPostWk80$gmt80ls > log10(5) & dataPostWk80$tx == "C3"]
dataPostWk80$epitope.dist.subtype.ls <- ifelse(dataPostWk80$protocol == "HVTN 703", dataPostWk80$epitope.dist.c.ls, dataPostWk80$epitope.dist.b.ls)


#week 104
dataWk104 <- read.csv(file.path(datDir, "amp_sieve_marks_wk104.csv")) %>%
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

escapeMark <- read.csv(file.path(datDir, "VRC01escape_DescFile.csv"))
dataWk104$ntx <- ifelse(dataWk104$tx == "C3", "Placebo", ifelse(dataWk104$tx == "T1", "VRC01 10 mg/kg", "VRC01 30 mg/kg"))
dataWk104$epitope.dist.subtype.ls <- ifelse(dataWk104$protocol == "HVTN 703", dataWk104$epitope.dist.c.ls, dataWk104$epitope.dist.b.ls)
dataWk104 <- left_join(dataWk104, escapeMark, by = "pub_id")


# For each quantitative mark, use Juraska et al 2013 method to fit the model and plot the PE by mark curve
i = 1
for (mark in quantMarks){

  xLim <- range(c(dataPostWk80[, mark], dataWk80[, mark]), na.rm=TRUE)
  
  dataPostWk80Sub <- dplyr::select(dataPostWk80, all_of(c("tx_pool", "hiv1fpday", "hiv1event", mark,  "stratVar")))
  dataWk80Sub <- dplyr::select(dataWk80, all_of(c("tx_pool", "hiv1fpday", "hiv1event", mark,  "stratVar")))
  
  colnames(dataPostWk80Sub) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
  colnames(dataWk80Sub) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
  
  dataPostWk80Sub$tx <- ifelse(as.character(dataPostWk80Sub$tx)=="C3", 0, 1)
  dataWk80Sub$tx <- ifelse(as.character(dataWk80Sub$tx)=="C3", 0, 1)
  
  # convert mark values for non-primary endpoints through tau to NA
  dataPostWk80Sub$mark <- ifelse(dataPostWk80Sub$eventInd==0, NA, dataPostWk80Sub$mark)
  dataWk80Sub$mark <- ifelse(dataWk80Sub$eventInd==0, NA, dataWk80Sub$mark)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  dataPostWk80Sub <- subset(dataPostWk80Sub, !(eventInd==1 & is.na(mark)))
  dataWk80Sub <- subset(dataWk80Sub, !(eventInd==1 & is.na(mark)))
  
  # fit the mark-specific HR model
  markRng <- range(dataPostWk80Sub$mark, na.rm=TRUE)
  markGrid <- seq(markRng[1], markRng[2], length.out=100)
  fitPostWk80 <- with(dataPostWk80Sub, sievePH(eventTime, eventInd, mark, tx, stratVar))
  sfitPostWk80 <- summary(fitPostWk80, markGrid=markGrid, sieveAlternative="oneSided")
  
  markRng <- range(dataWk80Sub$mark, na.rm=TRUE)
  markGrid <- seq(markRng[1], markRng[2], length.out=100)
  fitWk80 <- with(dataWk80Sub, sievePH(eventTime, eventInd, mark, tx, stratVar))
  sfitWk80 <- summary(fitWk80, markGrid=markGrid, sieveAlternative="oneSided")
  
 
  dataPostWk80Sub$treatment <- ifelse(dataPostWk80Sub$tx==0, "Placebo", "VRC01")
  dataWk80Sub$treatment <- ifelse(dataWk80Sub$tx==0, "Placebo", "VRC01")
  title <- ifelse(mark == "gmt80ls", "Most Resistant Founder", "Predicted Most Resistant Founder")
  p.PostWk80 <- sfitPostWk80$pWald.HRconstant.1sided.HRincrease
  fmt.p.PostWk80 <- ifelse(p.PostWk80<0.001, "< 0.001", paste0("= ", format(p.PostWk80, digits=2, nsmall=2)))
  p.Wk80 <- sfitWk80$pWald.HRconstant.1sided.HRincrease
  fmt.p.Wk80 <- ifelse(p.Wk80<0.001, "< 0.001", paste0("= ", format(p.Wk80, digits=2, nsmall=2)))
  
  subtitle <- paste0("Week 80 One-Sided Unadjusted Sieve P ", fmt.p.Wk80,
                     "\nPost-Week 80 One-Sided Unadjusted Sieve P ", fmt.p.PostWk80)
  
  sfitWk80$cohort <- "Week 80"
  sfitPostWk80$cohort <- "Post-Week 80"
  sfit = list(sfitWk80, sfitPostWk80)
  
  if(mark == "gmt80ls" | mark == "parscore2.ls"){
    xtickLab = c(0.07, 0.3, 1, 3, 5, 10, 50, 100)
    xtickAt = log10(xtickLab)
    xlim = c(min(markRng[1]) , markRng[2])
    
    p <- ggplotSievePH(sfit, data1 = dataWk80Sub, data2 = dataPostWk80Sub,  
                  xlab=xLabels[i],
                  breaks.x = xtickAt, labels.x = xtickLab, title=subtitle)
    
  }else if (mark == "parscore1.ls"){
   
    xtickAt = c(-1, -0.5, 0, 0.5, 1, 2,  3)
    xtickLab = round(exp(xtickAt)/(1+exp(xtickAt)),1)
    xlim = exp(c(min(markRng[1]) , markRng[2]))/(1+exp(c(min(markRng[1]) , markRng[2])))
    
    p <- ggplotSievePH(sfit, data1 = dataWk80Sub, data2 = dataPostWk80Sub,  
                       xlab=xLabels[i],
                       breaks.x = xtickAt, labels.x = xtickLab, title=subtitle)
    
    
  }else{
    p <- ggplotSievePH(sfit, data1 = dataWk80Sub, data2 = dataPostWk80Sub,  
                       xlab=xLabels[i], title=subtitle)
  }
  
#margin goes top, right, bottom, left  
  ggsave(file.path(figDir, paste0("703and704_sievePH_PEby", markFileString[i], "_VRC01_placebo.pdf")),
         plot=p, width=1 * 7, height=0.9 * 8)
  i = i + 1
  
}




#replot sievePH IC80 up to 100; combine week 80 and post week 80
mark = "gmt80ls"
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

# fit the mark-specific HR model
markRng <- range(dataPostWk80Sub$mark, na.rm=TRUE)
markGrid <- seq(markRng[1], markRng[2], length.out=100)
fitPostWk80 <- with(dataPostWk80Sub, sievePH(eventTime, eventInd, mark, tx, stratVar))
sfitPostWk80 <- summary(fitPostWk80, markGrid=markGrid, sieveAlternative="oneSided")

markRng <- range(dataWk80Sub$mark, na.rm=TRUE)
markGrid <- seq(markRng[1], markRng[2], length.out=100)
fitWk80 <- with(dataWk80Sub, sievePH(eventTime, eventInd, mark, tx, stratVar))
sfitWk80 <- summary(fitWk80, markGrid=markGrid, sieveAlternative="oneSided")


dataPostWk80Sub$treatment <- ifelse(dataPostWk80Sub$tx==0, "Placebo", "VRC01")
dataWk80Sub$treatment <- ifelse(dataWk80Sub$tx==0, "Placebo", "VRC01")
title <- ifelse(mark == "gmt80ls", "Most Resistant Founder", "Predicted Most Resistant Founder")
p.PostWk80 <- sfitPostWk80$pWald.HRconstant.1sided.HRincrease
fmt.p.PostWk80 <- ifelse(p.PostWk80<0.001, "< 0.001", paste0("= ", format(p.PostWk80, digits=2, nsmall=2)))
p.Wk80 <- sfitWk80$pWald.HRconstant.1sided.HRincrease
fmt.p.Wk80 <- ifelse(p.Wk80<0.001, "< 0.001", paste0("= ", format(p.Wk80, digits=2, nsmall=2)))

subtitle <- paste0("Week 80 One-Sided Unadjusted Sieve P ", fmt.p.Wk80,
                   "\nPost-Week 80 One-Sided Unadjusted Sieve P ", fmt.p.PostWk80)

sfitWk80$cohort <- "Week 80"
sfitPostWk80$cohort <- "Post-Week 80"
sfit = list(sfitWk80, sfitPostWk80)

xtickLab = c(0.07, 0.3, 1, 3, 5, 10, 25, 50, 100)
xtickAt = log10(xtickLab)
xlim = c(min(markRng[1]) , markRng[2])

p <- ggplotSievePH(sfit, data1 = dataWk80Sub, data2 = dataPostWk80Sub,  
                   xlab=expression(IC[80] ~ "(" * mu * "g/ml)"),
                   breaks.x = xtickAt, labels.x = xtickLab, title=subtitle)
ggsave(file.path(figDir, paste0("703and704_sievePH_PEbyIC80_ls_VRC01_placebo_v2.pdf")),
       plot=p, width=1 * 7, height=0.9 * 8)



#refit week 80, post week 80 week 104 using the nonparametric method
loadresult <- FALSE
defaultbandwidth <- FALSE
if(loadresult){
  if(defaultbandwidth){
    fitWk80 <- readRDS(file.path(outDir, "nonparametricfitWk80_defaultbw.rds"))
    fitWk104 <- readRDS(file.path(outDir, "nonparametricfitWk104_defaultbw.rds"))
    fitPostWk80 <- readRDS(file.path(outDir, "nonparametricfitPostWk80_defaultbw.rds"))
    
  }else{
    fitWk80 <- readRDS(file.path(outDir, "nonparametricfitWk80.rds"))
    fitWk104 <- readRDS(file.path(outDir, "nonparametricfitWk104.rds"))
    fitPostWk80 <- readRDS(file.path(outDir, "nonparametricfitPostWk80.rds"))
  }
  
  
}else{
  dataWk80Sub$eventTime2 <- pmax(dataWk80Sub$eventTime, 1) / 365.25
  dataPostWk80Sub$eventTime2 <- pmax(dataPostWk80Sub$eventTime, 1) / 365.25
  dataWk104Sub$eventTime2 <- pmax(dataWk104Sub$eventTime, 1) / 365.25
  
  if(defaultbandwidth){
    hbandic80 <- 5*sqrt(var((dataWk80Sub$mark - min(dataWk80Sub$mark, na.rm = TRUE))/(max(dataWk80Sub$mark, na.rm = TRUE) 
                  - min(dataWk80Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataWk80Sub$eventInd)^{-1/3}
    fitWk80 <- with(dataWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 85.9/52, tband = 0.5,
                                                hband = hbandic80, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    hbandic80 <- 5*sqrt(var((dataWk104Sub$mark - min(dataWk104Sub$mark, na.rm = TRUE))/(max(dataWk104Sub$mark, na.rm = TRUE) 
                 - min(dataWk104Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataWk104Sub$eventInd)^{-1/3}
    fitWk104 <- with(dataWk104Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 108.6/52, tband = 0.5,
                                                  hband = hbandic80, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    hbandic80 <- 5*sqrt(var((dataPostWk80Sub$mark - min(dataPostWk80Sub$mark, na.rm = TRUE))/(max(dataPostWk80Sub$mark, na.rm = TRUE) 
                 - min(dataPostWk80Sub$mark, na.rm = TRUE)), na.rm = TRUE))*sum(dataPostWk80Sub$eventInd)^{-1/3}
    fitPostWk80 <- with(dataPostWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 24/52, tband = 0.5,
                                                        hband = hbandic80, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    saveRDS(fitWk80, file.path(outDir, "nonparametricfitWk80_defaultbw.rds"))
    saveRDS(fitWk104, file.path(outDir, "nonparametricfitWk104_defaultbw.rds"))
    saveRDS(fitPostWk80, file.path(outDir, "nonparametricfitPostWk80_defaultbw.rds"))
  }else{
    fitWk80 <- with(dataWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 85.9/52, tband = 0.5,
                                                hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    
    fitWk104 <- with(dataWk104Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 108.6/52, tband = 0.5,
                                                  hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    fitPostWk80 <- with(dataPostWk80Sub, kernel_sievePH(eventTime2, eventInd, mark, tx, tau = 24/52, tband = 0.5,
                                                        hband = 0.42, nvgrid = 100, a = 0.01, b = 1, nboot = NULL, seed = 1))
    
    saveRDS(fitWk80, file.path(outDir, "nonparametricfitWk80.rds"))
    saveRDS(fitWk104, file.path(outDir, "nonparametricfitWk104.rds"))
    saveRDS(fitPostWk80, file.path(outDir, "nonparametricfitPostWk80.rds"))
    
  }
}

sfitWk80 <- summary(fitWk80, sieveAlternative = "oneSided")
sfitWk104 <- summary(fitWk104, sieveAlternative = "oneSided")
sfitPostWk80 <- summary(fitPostWk80, sieveAlternative = "oneSided")

p.PostWk80 <- sfitPostWk80$HRconstant.1sided$Tint2minc[2]
fmt.p.PostWk80 <- ifelse(p.PostWk80<0.001, "< 0.001", paste0("= ", format(p.PostWk80, digits=2, nsmall=2)))
p.Wk80 <- sfitWk80$HRconstant.1sided$Tint2minc[2]
fmt.p.Wk80 <- ifelse(p.Wk80<0.001, "< 0.001", paste0("= ", format(p.Wk80, digits=2, nsmall=2)))

subtitle <- paste0("Week 80 One-Sided Unadjusted Sieve P ", fmt.p.Wk80,
                   "\nPost-Week 80 One-Sided Unadjusted Sieve P ", fmt.p.PostWk80)

sfitWk80$cohort <- "Week 80"
sfitPostWk80$cohort <- "Post-Week 80"
sfit = list(sfitWk80, sfitPostWk80)
dataPostWk80Sub$treatment <- ifelse(dataPostWk80Sub$tx==0, "Placebo", "VRC01")
dataWk80Sub$treatment <- ifelse(dataWk80Sub$tx==0, "Placebo", "VRC01")
xtickLab = c(0.07, 0.3, 1, 3, 5, 10, 25, 50, 100)
xtickAt = log10(xtickLab)
xlim = c(min(markRng[1]) , markRng[2])


#plot Wk80 and Wk104 together
sfitWk80$cohort <- "Week 80"
sfitWk104$cohort <- "Week 104"
sfit = list(sfitWk80, sfitWk104)
dataWk104Sub$treatment <- ifelse(dataWk104Sub$tx==0, "Placebo", "VRC01")

p <- ggplotSievePH_wk80wk104(sfit, data1 = dataWk80Sub[, c("treatment", "mark", "tx")], 
                             data2 = dataWk104Sub[, c("treatment", "mark", "tx")],  
                   xlab=expression(IC[80] ~ "(" * mu * "g/ml)"),
                   breaks.x = xtickAt, labels.x = xtickLab, title="")

ggsave(file.path(figDir, paste0("703and704_sievePH_PEbyIC80_ls_VRC01_placebo_Wk80Wk104nonparametric", ifelse(defaultbandwidth, "_dbw",""), ".pdf")),
       plot=p, width=1 * 7, height=0.9 * 8)



#plot postWk80 seperately
sfitPostWk80 <- summary(fitPostWk80, sieveAlternative = "oneSided")
pPostWk80 <- ggplotSievePH_onecohort (sfitPostWk80, dataPostWk80Sub[, c("treatment", "mark", "tx")],
                                   xlab=expression(IC[80] ~ "(" * mu * "g/ml)"),
                                   breaks.x = xtickAt, labels.x = xtickLab, title="", color_line = "darkgreen",
                                   height.panels = c(1, 2), cohort = "Weeks 80 to 104", legendLabels = "Weeks 80 to 104",
                                   limits.x = range(dataWk104Sub$mark, na.rm =TRUE))


ggsave(file.path(figDir, paste0("703and704_sievePH_PEbyIC80_ls_VRC01_placebo_postWk80nonparametric", ifelse(defaultbandwidth, "_dbw",""), ".pdf")),
       plot=pPostWk80, width=1 * 7, height=0.9 * 8)

