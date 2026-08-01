#Instantaneous PE for sensitive and resistant viruses, VRC01 pooled vs. control, Week 104 analysis

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")
source(file.path(repoDir, "code/macro/instantaneousPEutils.R"))
library(tidyverse)
library(survival)
logit <- function(p){
  return(log(p/(1-p)))
}


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
dataWk104$gmt80ls2 <- dataWk104$gmt80ls
dataWk104$gmt80ls2[!is.na(dataWk104$gmt80ls) & 10^dataWk104$gmt80ls>=100] <- log10(100)

# four intermediate-to-resistant viruses, all in the VRC01 groups, that were diagnosed at the Week 96 or 104 blood draw
df <- filter(dataWk104, tx_pool!="C3" & (hiv1fpday>700) & !is.na(gmt80ls) & gmt80ls > 1)
df <- dplyr::select(df, all_of(c("protocol","tx", "hiv1fpday", "hiv1event", "gmt80ls",  "stratVar")))


#For sensitive viruses
dataWk104Sub <-  dplyr::select(dataWk104, all_of(c("tx_pool", "hiv1fpday", "hiv1event", "gmt80ls2",  "stratVar")))
colnames(dataWk104Sub) <- c("treatment", "time", "status", "mark", "stratVar")
dataWk104Sub <- subset(dataWk104Sub, !(status==1 & is.na(mark)))

#censor IC80>1 cases
dataWk104Sub$status [dataWk104Sub$status == 1 & dataWk104Sub$mark >= log10(1)] <- 0
dataWk104Sub$treatment <- ifelse(as.character(dataWk104Sub$treatment)=="C3", 1, 2)
dataWk104Sub$time <- dataWk104Sub$time/7

paste.p <- function(p){
  if(p<0.001){return(" < 0.001")
  }else{return(paste(" = ",as.character(format(p,digits=2))))}}

phTestPvalue <- cox.zph(coxph(Surv(time, status) ~ treatment, data = dataWk104Sub ))$table["treatment", "p"]
legendPhTest <- paste("PH Test P-value", paste.p(phTestPvalue))


pdf(file=file.path(figDir, "instPEsensitiveVirus_VRC01pooledVsControl_fixedBw.pdf"), 
    width=1.05 * 5, height=1.05 * 5)
fit1 <- plotSmoothHazTE (dataWk104Sub,  summaryType = "efficacy", 
                 title = expression("Dose-Pooled VRC01 vs. Placebo\n", IC[80] < 1 ~ "" * mu * "g/ml") , "Weeks since Enrollment", "Instantaneous Hazard-Ratio PE (%)", 
                 xtickAt = c(0, 16, 32, 48, 64, 80, 96, 108), xtickLabel = c(0, 16, 32, 48, 64, 80, 96, 108), 
                 xMax=108, showLegend=FALSE, phTestPvalue  = paste.p(phTestPvalue),
                 nBoot=NULL, loadFile=NULL, saveFile="instantPEsensitiveVirus_Wk104_VRC01pooledvsControl_fixedBw.rds", saveDir=outDir)

dev.off()

pdf(file=file.path(figDir, "instPEsensitiveVirus_VRC01pooledVsControl_optimBw.pdf"), 
    width=1.05 * 5, height=1.05 * 5)
fit2 <- plotSmoothHazTE (dataWk104Sub,  summaryType = "efficacy", 
                 title = expression("Dose-Pooled VRC01 vs. Placebo\n", IC[80] < 1 ~ "" * mu * "g/ml") , "Weeks since Enrollment", "Instantaneous Hazard-Ratio PE (%)", 
                 xtickAt = c(0, 16, 32, 48, 64, 80, 96, 108), xtickLabel = c(0, 16, 32, 48, 64, 80, 96, 108), 
                 xMax=108, showLegend=FALSE, phTestPvalue  = paste.p(phTestPvalue),
                 nBoot=500, loadFile=NULL, saveFile="instantPEsensitiveVirus_Wk104_VRC01pooledvsControl_optimBw.rds", saveDir=outDir)

dev.off()



#For resistant viruses
dataWk104Sub <-  dplyr::select(dataWk104, all_of(c("tx_pool", "hiv1fpday", "hiv1event", "gmt80ls2",  "stratVar")))
colnames(dataWk104Sub) <- c("treatment", "time", "status", "mark", "stratVar")
dataWk104Sub <- subset(dataWk104Sub, !(status==1 & is.na(mark)))

#censor IC80>1 cases
dataWk104Sub$status [dataWk104Sub$status == 1 & dataWk104Sub$mark < log10(1)] <- 0

dataWk104Sub$treatment <- ifelse(as.character(dataWk104Sub$treatment)=="C3", 1, 2)
dataWk104Sub$time <- dataWk104Sub$time/7

phTestPvalue <- cox.zph(coxph(Surv(time, status) ~ treatment, data = dataWk104Sub ))$table["treatment", "p"]

#For resistant viruses
pdf(file=file.path(figDir, "instPEresistantVirus_VRC01pooledVsControl_fixedBw.pdf"), 
    width=1.05 * 5, height=1.05 * 5)
fit4 <- plotSmoothHazTE (dataWk104Sub,  summaryType = "efficacy", 
                 title = expression("Dose-Pooled VRC01 vs. Placebo\n", IC[80] > 1 ~ "" * mu * "g/ml") , "Weeks since Enrollment", "Instantaneous Hazard-Ratio PE (%)", 
                 xtickAt = c(0, 16, 32, 48, 64, 80, 96, 108), xtickLabel = c(0, 16, 32, 48, 64, 80, 96, 108), 
                 xMax=108, showLegend=FALSE, phTestPvalue  = paste.p(phTestPvalue),
                 nBoot=NULL, loadFile=NULL, saveFile="instantPEresistantVirus_Wk104_VRC01pooledvsControl_fixedBw.rds", saveDir=outDir)

dev.off()

pdf(file=file.path(figDir, "instPEresistantVirus_VRC01pooledVsControl_optimBw.pdf"), 
    width=1.05 * 5, height=1.05 * 5)
fit5 <- plotSmoothHazTE (dataWk104Sub,  summaryType = "efficacy", 
                 title = expression("Dose-Pooled VRC01 vs. Placebo\n", IC[80] > 1 ~ "" * mu * "g/ml") , "Weeks since Enrollment", "Instantaneous Hazard-Ratio PE (%)", 
                 xtickAt = c(0, 16, 32, 48, 64, 80, 96, 108), xtickLabel = c(0, 16, 32, 48, 64, 80, 96, 108), 
                 xMax=108, showLegend=FALSE, phTestPvalue  = paste.p(phTestPvalue),
                 nBoot=500, loadFile=NULL, saveFile="instantPEresistantVirus_Wk104_VRC01pooledvsControl_optimBw.rds", saveDir=outDir)

dev.off()



#bandwidth
fit1[1, ]
fit2[1, ]

