library(here)
here::i_am("README.md")
repoDir <- here::here()
macroDir <- file.path(repoDir, "code/macro")
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")

# load macros
source(file.path(macroDir, "phReg.R")) 
source(file.path(macroDir, "cuminc_functions.R"))

# input file name
dataFile <- file.path(datDir, "amp_survival_postwk80.csv")

# output file names
CIR.csvFile_summary_703 <- file.path(tabDir, "v703_cir_efficacy_postwk80_trunc.csv")
CIR.csvFile_summary_704 <- file.path(tabDir, "v704_cir_efficacy_postwk80_trunc.csv")

# specify variables names from input dataset

# groupings to be compared
grpVar_pool  <- "rx_pool"
grpVar_ind <- "rx_code"

# follow-up time information
timeVar <- "fudays_postwk80"

# event indicator
eventIndVar <- "status_postwk80"

# unique identifier
idVar <- "pub_id"

# the strata variable
strataVar <- "rx_code"

# reference level of your group variable
refLvl <- "C3"

# comparison level of your group variable
cmpLvl_pool <- c("T1+T2")
cmpLvl_ind <- c("T1", "T2")

for( trial_name in c('v703', 'v704') ) {

  # set output file
  CIR.csvFile_summary = ifelse(trial_name=='v703',
                               CIR.csvFile_summary_703,
                               CIR.csvFile_summary_704)

  # source input data and subset to get records only for MITT participants
  if( trial_name=='v703') {
    dat <- read.csv( dataFile, stringsAsFactors = FALSE )
    dat <- subset(dat, protocol=='HVTN 703')
    mitt <- subset(dat, subset=(efficacy_flag == 1), select = c(idVar, grpVar_pool, grpVar_ind, timeVar, eventIndVar, strataVar))
  } else {
    dat <- read.csv( dataFile, stringsAsFactors = FALSE )
    dat <- subset(dat, protocol=='HVTN 704')
    mitt <- subset(dat, subset=(efficacy_flag == 1), select = c(idVar, grpVar_pool, grpVar_ind, timeVar, eventIndVar, strataVar))
  }
  
  # create strata weights
  
  strataWts <- c(C3 = 1, T2=0.5, T1=0.5)
  
  # get pooled and individual cumulative incidence estimates
  # remove censoring to get estiamtes through the last infection
  mitt.cuminc_pool <- 
    naCumInc( 
      data = mitt,
      futimeVar = timeVar, 
      eventVar = eventIndVar, 
      groupVar = grpVar_pool, 
      strataVar = strataVar,
      idVar = idVar,
      #censor = list( minAtRisk=150),
      stratifiedEstimator = TRUE,
      strataWeights = strataWts)
  
  mitt.cuminc_ind <- 
    naCumInc( 
      data = mitt,
      futimeVar = timeVar, 
      eventVar = eventIndVar, 
      groupVar = grpVar_ind, 
#      censor = list( minAtRisk=150),
      idVar = idVar)

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
    CIR.summary$temp1 <- c(sum(mitt.CIR_pool$data[["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2", "C3"),][["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1", "C3"),][["orig_status_postwk80"]]))
    CIR.summary$temp2 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("T1+T2"),][["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2"),][["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1"),][["orig_status_postwk80"]]))
    CIR.summary$temp3 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("C3"),][["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["orig_status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["orig_status_postwk80"]]))
    CIR.summary$n_cases_total <- with(CIR.summary, paste0(temp1, " (", temp2, " vs. ", temp3, ")"))
    CIR.summary$temp1 <- CIR.summary$temp2 <- CIR.summary$temp3 <- NULL
    CIR.summary$tau_wks <- c(mitt.CIR_pool$censorInfo$censorTime, 
                             mitt.CIR_ind$censorInfo$censorTime,
                             mitt.CIR_ind$censorInfo$censorTime)/7
    CIR.summary$temp1 <- c(sum(mitt.CIR_pool$data[["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2", "C3"),][["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1", "C3"),][["status_postwk80"]]))
    CIR.summary$temp2 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("T1+T2"),][["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T2"),][["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("T1"),][["status_postwk80"]]))
    CIR.summary$temp3 <- c(sum(mitt.CIR_pool$data[mitt.CIR_pool$data$rx_pool %in% c("C3"),][["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["status_postwk80"]]),
                           sum(mitt.CIR_ind$data[mitt.CIR_ind$data$rx_code %in% c("C3"),][["status_postwk80"]]))
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
    

  write.csv( CIR.summary, file=CIR.csvFile_summary, na="", row.names=FALSE, quote=FALSE)
}  

q(save = "no")
