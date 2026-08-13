library(tidyverse)
library(here)
here::i_am("README.md")
repoDir <- here::here()

outDatDir <- file.path(repoDir, "data")
datDir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data_final" # file.path(repoDir, "data")
macroDir <- file.path(repoDir, 'code/macro')

# load macros
source(file.path(macroDir, "cuminc_functions.R"))

# input files
dataFile <- file.path(datDir, "amp_survival_wk104.csv")
sDataFile <- file.path(datDir, 'amp_sieve_marks_wk104.csv')

# output file names
CIR_pool_80.csvFile <- file.path(outDatDir, "amp_cir_wk104_pool_cmpriskIC80ls_2cat_trunc.csv")

# specify variables names from dataset

# name of variable containing groupings to be compared
grpVar_ind  <- "tx"
grpVar_pool <- "tx_pool" 

# name of variable containing follow-up time information
timeVar <- "hiv1survday"

# name of variable containing info on which type of event has occurred
# It's best to make this into a factor after it's read in, being careful to set the
# eventType that corresponds to 'censoring' to be the first level of the factor.
# the second 'v2' version combines sensitive and moderately sensitive into a single level
eventTypeVar80 <- "hiv1event80ls" 
eventTypeVar80_v2 <- "hiv1event80ls_v2" 

# idVar - name of variable containing a unique identifier
idVar <- "pub_id"

# reference level of your group variable (e.g reference trt group) 
refLvl <- "C3"

# comparison level of your group variable
cmpLvl_pool <- c("T1+T2")
cmpLvl_ind <- c("T1", "T2")

# read in time-to-event dataset and sieve data file
dat <- read.csv( dataFile, stringsAsFactors = FALSE ) %>%
  select(protocol, pub_id, tx, tx_pool, hiv1survday)
sdat = read.csv( sDataFile ) %>%
  select(pub_id, hiv1event, gmt80ls) 

tmp = merge(dat, sdat) %>%
  mutate(gmt80ls = as.numeric(sub('>','',gmt80ls)))

tmp$hiv1event80ls_v2 = case_when(
  tmp$hiv1event==0 ~ 0,
  tmp$gmt80ls <= 1 ~ 1,
  tmp$gmt80ls > 1 ~ 2,
  TRUE ~ -1
)

dat = tmp


# exclude cases with missing breakthrough data
dat <- subset(dat, hiv1event80ls_v2!=-1)

table(dat$tx, dat$hiv1event80ls_v2, useNA="always")

## Create eventTypeVar
dat[[ eventTypeVar80_v2 ]] <- 
  factor( dat[[ eventTypeVar80_v2 ]], 
          levels= c(0, 1, 2), 
          labels = c("Censored", "Sensitive", "Resistant") )

table(dat$tx, dat$hiv1event80ls_v2, useNA="always")

## get cumulative incidence via competing risks for pooled VRC01. 
## Note our input dataset is already censored how we need it to be.

cuminc_pool_80 <- 
  cmpRisks(data=dat, 
           futimeVar=timeVar, 
           groupVar= grpVar_pool, 
           eventVar= eventTypeVar80_v2)

cuminc_ind_80 <- 
  cmpRisks(data=dat, 
           futimeVar=timeVar, 
           groupVar= grpVar_ind,
           eventVar= eventTypeVar80_v2)

## get CIR/PE estimates for each cohort
CIR_pool_80 <- EffCIR(cuminc_pool_80, refLvl= refLvl, cmpLvl= cmpLvl_pool)
CIR_ind_80  <- EffCIR(cuminc_ind_80,  refLvl = refLvl, cmpLvl=cmpLvl_ind)

# aggregate results for nicer output


# then for IC80
CIR.summary.80 <- data.frame(comparison = c(rep("VRC01 Pooled vs. Control", 2), 
                                            rep("VRC01 30 mg/kg vs. Control", 2),
                                            rep("VRC01 10 mg/kg vs. Control", 2)),
                             eventType = rep(c("Sensitive", "Resistant"), 3),
                             tau_wks = c(rep(max(dat[dat$hiv1event==0,]$hiv1survday), 6)/7),
                             n_cases_tau = NA,
                             eff = NA,
                             lo.eff  = NA,
                             up.eff = NA)
CIR.summary.80$temp1 <- c(sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Sensitive")),
                          sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Resistant")),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("T2", "C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("T2", "C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("T1", "C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("T1", "C3"))))
CIR.summary.80$temp2 <- c(sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_pool_80$data[["tx_pool"]] %in% c("T1+T2"))),
                          sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_pool_80$data[["tx_pool"]] %in% c("T1+T2"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("T2"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("T2"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("T1"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("T1"))))
CIR.summary.80$temp3 <- c(sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_pool_80$data[["tx_pool"]] %in% c("C3"))),
                          sum(c(CIR_pool_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_pool_80$data[["tx_pool"]] %in% c("C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Sensitive" & CIR_ind_80$data[["tx"]] %in% c("C3"))),
                          sum(c(CIR_ind_80$data[["hiv1event80ls_v2"]]=="Resistant" & CIR_ind_80$data[["tx"]] %in% c("C3"))))
CIR.summary.80$n_cases_tau <- with(CIR.summary.80, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
CIR.summary.80$temp1 <- CIR.summary.80$temp2 <- CIR.summary.80$temp3 <- NULL
CIR.summary.80$eff <- c(CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Sensitive"),]$eff,
                        CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Resistant"),]$eff,
                        CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$eff,
                        CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$eff,
                        CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$eff,
                        CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$eff)
CIR.summary.80$lo.eff <- c(CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Sensitive"),]$lo.eff,
                           CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Resistant"),]$lo.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$lo.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$lo.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$lo.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$lo.eff)
CIR.summary.80$up.eff <- c(CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Sensitive"),]$up.eff,
                           CIR_pool_80$eff[CIR_pool_80$eff$eventType %in% c("Resistant"),]$up.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$up.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T2 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$up.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Sensitive"),]$up.eff,
                           CIR_ind_80$eff[CIR_ind_80$eff$comparison %in% c("T1 vs. C3") & CIR_ind_80$eff$eventType %in% c("Resistant"),]$up.eff)

## output CIR/PE estimates
write.csv(CIR_pool_80$CIR, file=CIR_pool_80.csvFile, na="", row.names=FALSE, quote=FALSE)

q(save = "no")

