library(here)
here::i_am("README.md")
repoDir <- here::here()
macroDir <- file.path(repoDir, "code/macro")
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data_final" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")

# load macros
source(file.path(macroDir, "phReg.R")) 
source(file.path(macroDir, "cuminc_functions.R"))

# input file name
dataFile <- file.path(dat2Dir, "v703_survival.csv")

# output file names
CIR.csvFile_summary_trunc <- file.path(tabDir, "v703_cir_efficacy_wk104_trunc.csv")

# specify variables names from input dataset

# groupings to be compared
grpVar_pool  <- "rx_pool"
grpVar_ind <- "rx_code"

# follow-up time information
timeVar <- "fudayswk104"
timeVar_trunc <- "fudayswk104_trunc"

# event indicator
eventIndVar <- "statuswk104"

# unique identifier
idVar <- "pub_id"

# the strata variable
strataVar <- "rx_code"

# reference level of your group variable
refLvl <- "C3"

# comparison level of your group variable
cmpLvl_pool <- c("T1+T2")
cmpLvl_ind <- c("T1", "T2")

# source input data and subset to get records only for MITT participants
dat <- read.csv( dataFile, stringsAsFactors = FALSE )
mitt <- subset(dat, subset=(efficacy_flag == 1), select = c(idVar, grpVar_pool, grpVar_ind, timeVar, eventIndVar, strataVar) )

# create strata weights
    
strataWts <- c(C3 = 1, T2=0.5, T1=0.5)

# get pooled and individual cumulative incidence estimates
# uncensored, tau censored and truncated
mitt.cuminc_pool <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar, 
    eventVar = eventIndVar, 
    groupVar = grpVar_pool, 
    strataVar = strataVar,
    idVar = idVar,
    stratifiedEstimator = TRUE,
    strataWeights = strataWts)

mitt.cuminc_pool_tau <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar, 
    eventVar = eventIndVar, 
    groupVar = grpVar_pool, 
    strataVar = strataVar,
    idVar = idVar,
    censor = list( minAtRisk=150),
    stratifiedEstimator = TRUE,
    strataWeights = strataWts)

# truncate the followup time at the censored tau time and generate the
# truncated cumulative estimate values
stopifnot(mitt.cuminc_pool_tau$censorInfo$censorTime==760)
mitt$fudayswk104_trunc = pmin(mitt.cuminc_pool_tau$censorInfo$censorTime,
                              mitt$fudayswk104)

mitt.cuminc_pool_trunc <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar_trunc, 
    eventVar = eventIndVar, 
    groupVar = grpVar_pool, 
    strataVar = strataVar,
    idVar = idVar,
    stratifiedEstimator = TRUE,
    strataWeights = strataWts)


mitt.cuminc_ind <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar, 
    eventVar = eventIndVar, 
    groupVar = grpVar_ind, 
    idVar = idVar)

mitt.cuminc_ind_tau <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar, 
    eventVar = eventIndVar, 
    groupVar = grpVar_ind, 
    idVar = idVar,
    censor = list( minAtRisk=150))

mitt.cuminc_ind_trunc <- 
  naCumInc( 
    data = mitt,
    futimeVar = timeVar_trunc, 
    eventVar = eventIndVar, 
    groupVar = grpVar_ind, 
    idVar = idVar)

# get pooled and individual CIR/PE estimates 
mitt.CIR_pool_trunc <- EffCIR( mitt.cuminc_pool, refLvl = refLvl, cmpLvl=cmpLvl_pool, nullHypEff=0)

mitt.CIR_pool_trunc$tests$cmpLvl <- mitt.CIR_pool_trunc$eff$cmpLvl
mitt.CIR_pool_trunc$tests$refLvl <- mitt.CIR_pool_trunc$eff$refLvl
mitt.CIR_pool_trunc$tests$comparison <- paste0(mitt.CIR_pool_trunc$eff$cmpLvl, " vs. ", mitt.CIR_pool_trunc$eff$refLvl)
mitt.CIR_pool_trunc$tests <- mitt.CIR_pool_trunc$tests[,c(7,8,9,1,2,3,4,5,6)]

mitt.CIR_ind_trunc  <- EffCIR( mitt.cuminc_ind,  refLvl = refLvl, cmpLvl=cmpLvl_ind, nullHypEff=0)

mitt.CIR_ind_trunc$tests$cmpLvl <- mitt.CIR_ind_trunc$eff$cmpLvl
mitt.CIR_ind_trunc$tests$refLvl <- mitt.CIR_ind_trunc$eff$refLvl
mitt.CIR_ind_trunc$tests$comparison <- paste0(mitt.CIR_ind_trunc$eff$cmpLvl, " vs. ", mitt.CIR_ind_trunc$eff$refLvl)
mitt.CIR_ind_trunc$tests <- mitt.CIR_ind_trunc$tests[,c(7,8,9,1,2,3,4,5,6)]


# aggregate results for nicer output 

# CIR/PE
CIR.summary_trunc <- data.frame(comparison = c("Pooled VRC01 vs. Control", "VRC01 30 mg/kg vs. Control", "VRC01 10 mg/kg vs. Control"),
                                tau_wks = NA,
                                n_cases_total = NA,
                                pe = NA,
                                bounds = NA,
                                pvalue = NA,
                                adj.pvalue = NA)
CIR.summary_trunc$temp1 <- c(sum(mitt.CIR_pool_trunc$data[["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T2", "C3"),][["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T1", "C3"),][["orig_statuswk104"]]))
CIR.summary_trunc$temp2 <- c(sum(mitt.CIR_pool_trunc$data[mitt.CIR_pool_trunc$data$rx_pool %in% c("T1+T2"),][["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T2"),][["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T1"),][["orig_statuswk104"]]))
CIR.summary_trunc$temp3 <- c(sum(mitt.CIR_pool_trunc$data[mitt.CIR_pool_trunc$data$rx_pool %in% c("C3"),][["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("C3"),][["orig_statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("C3"),][["orig_statuswk104"]]))
CIR.summary_trunc$n_cases_total <- with(CIR.summary_trunc, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
CIR.summary_trunc$temp1 <- CIR.summary_trunc$temp2 <- CIR.summary_trunc$temp3 <- NULL
CIR.summary_trunc$tau_wks <- c(mitt.cuminc_pool_tau$censorInfo$censorTime,
                               mitt.cuminc_ind_tau$censorInfo$censorTime,
                               mitt.cuminc_ind_tau$censorInfo$censorTime)/7
CIR.summary_trunc$temp1 <- c(sum(mitt.CIR_pool_trunc$data[["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T2", "C3"),][["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T1", "C3"),][["statuswk104"]]))
CIR.summary_trunc$temp2 <- c(sum(mitt.CIR_pool_trunc$data[mitt.CIR_pool_trunc$data$rx_pool %in% c("T1+T2"),][["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T2"),][["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("T1"),][["statuswk104"]]))
CIR.summary_trunc$temp3 <- c(sum(mitt.CIR_pool_trunc$data[mitt.CIR_pool_trunc$data$rx_pool %in% c("C3"),][["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("C3"),][["statuswk104"]]),
                             sum(mitt.CIR_ind_trunc$data[mitt.CIR_ind_trunc$data$rx_code %in% c("C3"),][["statuswk104"]]))
CIR.summary_trunc$n_cases_total <- with(CIR.summary_trunc, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
CIR.summary_trunc$temp1 <- CIR.summary_trunc$temp2 <- CIR.summary_trunc$temp3 <- NULL
CIR.summary_trunc$pe <- paste0(round(c(mitt.CIR_pool_trunc$eff$eff,
                                       mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T2 vs. C3"),]$eff,
                                       mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T1 vs. C3"),]$eff)*100, 2), "%")
CIR.summary_trunc$temp1 <- round(c(mitt.CIR_pool_trunc$eff$lo.eff,
                                   mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T2 vs. C3"),]$lo.eff,
                                   mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T1 vs. C3"),]$lo.eff)*100, 2)
CIR.summary_trunc$temp2 <- round(c(mitt.CIR_pool_trunc$eff$up.eff,
                                   mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T2 vs. C3"),]$up.eff,
                                   mitt.CIR_ind_trunc$eff[mitt.CIR_ind_trunc$eff$comparison %in% c("T1 vs. C3"),]$up.eff)*100, 2)
CIR.summary_trunc$bounds <- with(CIR.summary_trunc, paste0("(", temp1, "% to ", temp2, "%)"))
CIR.summary_trunc$temp1 <- CIR.summary_trunc$temp2 <- NULL
CIR.summary_trunc$pvalue <- round(c(mitt.CIR_pool_trunc$tests$pvalue,
                                    mitt.CIR_ind_trunc$tests[mitt.CIR_ind_trunc$tests$comparison %in% c("T2 vs. C3"),]$pvalue,
                                    mitt.CIR_ind_trunc$tests[mitt.CIR_ind_trunc$tests$comparison %in% c("T1 vs. C3"),]$pvalue), 4)
# for holm-bonferroni adjustment must be careful about order of pvalues
pvalue.in <- mitt.CIR_ind_trunc$tests$pvalue #T1 then T2
pvalue.out <- round(p.adjust(pvalue.in, method="holm"), 4) #T1 then T2
CIR.summary_trunc$adj.pvalue <- c(" ", pvalue.out[2], pvalue.out[1])

# output PE estimate on truncated data
write.csv( CIR.summary_trunc, file=CIR.csvFile_summary_trunc, na="", row.names=FALSE, quote=FALSE)


q(save = "no")
