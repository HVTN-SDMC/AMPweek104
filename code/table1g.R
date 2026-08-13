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
dataFile <- file.path(dat2Dir, "v704_survival.csv")

# output file names
CIR.csvFile_summary <- file.path(tabDir, "v704_cir_efficacy_wk80.csv")

# specify variables names from input dataset

# groupings to be compared
grpVar_pool  <- "rx_pool"
grpVar_ind <- "rx_code"

# follow-up time information
timeVar <- "fudayswk80"

# event indicator
eventIndVar <- "statuswk80"

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
mitt.cuminc_pool <- 
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

mitt.cuminc_ind <- 
  naCumInc( 
     data = mitt,
     futimeVar = timeVar, 
     eventVar = eventIndVar, 
     groupVar = grpVar_ind, 
     idVar = idVar,
     censor = list( minAtRisk=150))

# get pooled and individual CIR/PE estimates 
mitt.CIR_pool <- EffCIR( mitt.cuminc_pool, refLvl = refLvl, cmpLvl=cmpLvl_pool, nullHypEff=0)

  mitt.CIR_pool$tests$cmpLvl <- mitt.CIR_pool$eff$cmpLvl
  mitt.CIR_pool$tests$refLvl <- mitt.CIR_pool$eff$refLvl
  mitt.CIR_pool$tests$comparison <- paste0(mitt.CIR_pool$eff$cmpLvl, " vs. ", mitt.CIR_pool$eff$refLvl)
  mitt.CIR_pool$tests <- mitt.CIR_pool$tests[,c(7,8,9,1,2,3,4,5,6)]

mitt.CIR_ind  <- EffCIR( mitt.cuminc_ind,  refLvl = refLvl, cmpLvl=cmpLvl_ind, nullHypEff=0)

  mitt.CIR_ind$tests$cmpLvl <- mitt.CIR_ind$eff$cmpLvl
  mitt.CIR_ind$tests$refLvl <- mitt.CIR_ind$eff$refLvl
  mitt.CIR_ind$tests$comparison <- paste0(mitt.CIR_ind$eff$cmpLvl, " vs. ", mitt.CIR_ind$eff$refLvl)
  mitt.CIR_ind$tests <- mitt.CIR_ind$tests[,c(7,8,9,1,2,3,4,5,6)]
  
# get pooled and individual cox estimates 
mitt.cox_pool <- 
  phReg( 
     mitt, 
     futimeVar = timeVar, 
     eventVar  = eventIndVar,
     groupVar  = grpVar_pool, 
     refLvl    = refLvl, 
     cmpLvl    = cmpLvl_pool,
     cox.zph   = TRUE,
     nullHypEff = 0)

  mitt.cox_pool$tests$cmpLvl <- mitt.cox_pool$eff$cmpLvl
  mitt.cox_pool$tests$refLvl <- mitt.cox_pool$eff$refLvl
  mitt.cox_pool$tests$comparison <- paste0(mitt.cox_pool$eff$cmpLvl, " vs. ", mitt.cox_pool$eff$refLvl)
  mitt.cox_pool$tests <- mitt.cox_pool$tests[,c(7,8,9,1,2,3,4,5,6)]

mitt.cox_ind_B <- 
  phReg( 
     mitt, 
     futimeVar = timeVar, 
     eventVar  = eventIndVar,
     groupVar  = grpVar_ind, 
     refLvl    = refLvl, 
     cmpLvl    = cmpLvl_ind[2],
     cox.zph   = TRUE,
     nullHypEff = 0)

  mitt.cox_ind_B$tests$cmpLvl <- mitt.cox_ind_B$eff$cmpLvl
  mitt.cox_ind_B$tests$refLvl <- mitt.cox_ind_B$eff$refLvl
  mitt.cox_ind_B$tests$comparison <- paste0(mitt.cox_ind_B$eff$cmpLvl, " vs. ", mitt.cox_ind_B$eff$refLvl)
  mitt.cox_ind_B$tests <- mitt.cox_ind_B$tests[,c(7,8,9,1,2,3,4,5,6)]
  
mitt.cox_ind_C <- 
  phReg( 
     mitt, 
     futimeVar = timeVar, 
     eventVar  = eventIndVar,
     groupVar  = grpVar_ind, 
     refLvl    = refLvl, 
     cmpLvl    = cmpLvl_ind[1],
     cox.zph   = TRUE,
     nullHypEff = 0)

  mitt.cox_ind_C$tests$cmpLvl <- mitt.cox_ind_C$eff$cmpLvl
  mitt.cox_ind_C$tests$refLvl <- mitt.cox_ind_C$eff$refLvl
  mitt.cox_ind_C$tests$comparison <- paste0(mitt.cox_ind_C$eff$cmpLvl, " vs. ", mitt.cox_ind_C$eff$refLvl)
  mitt.cox_ind_C$tests <- mitt.cox_ind_C$tests[,c(7,8,9,1,2,3,4,5,6)]

# aggregate results for nicer output (only need score test for cox model)
  
  # CIR/PE
  CIR.summary <- data.frame(comparison = c("Pooled VRC01 vs. Control", "VRC01 30 mg/kg vs. Control", "VRC01 10 mg/kg vs. Control"),
                            n_cases_total = NA,
                            tau_wks = NA,
                            n_cases_tau = NA,
                            pe = NA,
                            bounds = NA,
                            pvalue = NA,
                            adj.pvalue = NA)
  CIR.summary$temp1 <- c(sum(mitt.CIR_pool$data[["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2", "C3"),][["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1", "C3"),][["orig_statuswk80"]]))
  CIR.summary$temp2 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("T1+T2"),][["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2"),][["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1"),][["orig_statuswk80"]]))
  CIR.summary$temp3 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("C3"),][["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["orig_statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["orig_statuswk80"]]))
  CIR.summary$n_cases_total <- with(CIR.summary, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
  CIR.summary$temp1 <- CIR.summary$temp2 <- CIR.summary$temp3 <- NULL
  CIR.summary$tau_wks <- c(mitt.CIR_pool$censorInfo$censorTime, 
                           mitt.CIR_ind$censorInfo$censorTime,
                           mitt.CIR_ind$censorInfo$censorTime)/7
  CIR.summary$temp1 <- c(sum(mitt.CIR_pool$data[["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2", "C3"),][["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1", "C3"),][["statuswk80"]]))
  CIR.summary$temp2 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("T1+T2"),][["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2"),][["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1"),][["statuswk80"]]))
  CIR.summary$temp3 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("C3"),][["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["statuswk80"]]),
                         sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["statuswk80"]]))
  CIR.summary$n_cases_tau <- with(CIR.summary, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
  CIR.summary$temp1 <- CIR.summary$temp2 <- CIR.summary$temp3 <- NULL
  CIR.summary$pe <- paste0(round(c(mitt.CIR_pool$eff$eff,
                                   mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T2 vs. C3"),]$eff,
                                   mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T1 vs. C3"),]$eff)*100, 2), "%")
  CIR.summary$temp1 <- round(c(mitt.CIR_pool$eff$lo.eff,
                               mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T2 vs. C3"),]$lo.eff,
                               mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T1 vs. C3"),]$lo.eff)*100, 2)
  CIR.summary$temp2 <- round(c(mitt.CIR_pool$eff$up.eff,
                               mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T2 vs. C3"),]$up.eff,
                               mitt.CIR_ind$eff[mitt.CIR_ind$eff$comparison %in% c("T1 vs. C3"),]$up.eff)*100, 2)
  CIR.summary$bounds <- with(CIR.summary, paste0("(", temp1, "% to ", temp2, "%)"))
  CIR.summary$temp1 <- CIR.summary$temp2 <- NULL
  CIR.summary$pvalue <- round(c(mitt.CIR_pool$tests$pvalue,
                                mitt.CIR_ind$tests[mitt.CIR_ind$tests$comparison %in% c("T2 vs. C3"),]$pvalue,
                                mitt.CIR_ind$tests[mitt.CIR_ind$tests$comparison %in% c("T1 vs. C3"),]$pvalue), 4)
  # for holm-bonferroni adjustment must be careful about order of pvalues
  pvalue.in <- mitt.CIR_ind$tests$pvalue #T1 then T2
  pvalue.out <- round(p.adjust(pvalue.in, method="holm"), 4) #T1 then T2
  CIR.summary$adj.pvalue <- c(" ", pvalue.out[2], pvalue.out[1])

  # Cox PE
  cox.summary <- CIR.summary[,c("comparison", "n_cases_total")]
  cox.summary$pvalue.zph <- round(c(mitt.cox_pool$zph$p[1], mitt.cox_ind_B$zph$p[1], mitt.cox_ind_C$zph$p[1]), 4)
  cox.summary$pe <- paste0(round(c(mitt.cox_pool$eff$eff, mitt.cox_ind_B$eff$eff, mitt.cox_ind_C$eff$eff)*100, 2), "%")
  cox.summary$temp1 <- paste0(round(c(mitt.cox_pool$eff$lo.eff, mitt.cox_ind_B$eff$lo.eff, mitt.cox_ind_C$eff$lo.eff)*100, 2))
  cox.summary$temp2 <- paste0(round(c(mitt.cox_pool$eff$up.eff, mitt.cox_ind_B$eff$up.eff, mitt.cox_ind_C$eff$up.eff)*100, 2))
  cox.summary$bounds <- with(cox.summary, paste0("(", temp1, "% to ", temp2, "%)"))
  cox.summary$temp1 <- cox.summary$temp2 <- NULL
  cox.summary$pvalue <- c(subset(mitt.cox_pool$tests, test %in% c("Score"))$pvalue, 
                          subset(mitt.cox_ind_B$tests, test %in% c("Score"))$pvalue,
                          subset(mitt.cox_ind_C$tests, test %in% c("Score"))$pvalue)
  # for holm-bonferroni adjustment must be careful about order of pvalues
  pvalue.in <- cox.summary$pvalue[2:3] #T2 then T1
  pvalue.out <- round(p.adjust(pvalue.in, method="holm"), 4) #T2 then T1
  cox.summary$pvalue <- round(cox.summary$pvalue, 4)
  cox.summary$adj.pvalue <- c(" ", pvalue.out[1], pvalue.out[2])
  
write.csv( CIR.summary, file=CIR.csvFile_summary, na="", row.names=FALSE, quote=FALSE)

q(save = "no")
