library(here)
here::i_am("README.md")
repoDir <- here::here()

outDatDir <- file.path(repoDir, "data")
datDir <- file.path(repoDir, "data")
macroDir <- file.path(repoDir, 'code/macro')

# load macros
source(file.path(macroDir, "phReg.R")) 
source(file.path(macroDir, "cuminc_functions.R"))

# input file name
dataFile <- file.path(datDir, "amp_survival_postwk80.csv")

# output file names
CIR.csvFile_pool  <- file.path(outDatDir, "amp_cir_pool_postwk80_v2.csv")
CIR.csvFile_ind   <- file.path(outDatDir, "amp_cir_ind_postwk80_v2.csv")

# specify variables names from input dataset

# groupings to be compared
grpVar_pool  <- "rx_pool"
grpVar_ind <- "rx_code"

# follow-up time information
timeVar <- "fudays_postwk80"
timeVar_trunc <- "fudays_postwk80_trunc"

# event indicator
eventIndVar <- "status_postwk80"

# unique identifier
idVar <- "pub_id"

# the strata variable
strataVar_pool <- "prot_rx"
strataVar_ind <- "protocol"

# reference level of your group variable
refLvl <- "C3"

# comparison level of your group variable
cmpLvl_pool <- c("T1+T2")
cmpLvl_ind <- c("T1", "T2")

# source input data and subset to get records only for MITT participants
dat <- read.csv( dataFile, stringsAsFactors = FALSE )
mitt <- subset(dat, subset=(efficacy_flag == 1), select = c(idVar, grpVar_pool, grpVar_ind, timeVar, eventIndVar, strataVar_pool, strataVar_ind))

# create strata weights
strataWts_pool <- c("HVTN 703 C3" = 0.5, "HVTN 704 C3" = 0.5, 
                    "HVTN 703 T2" = 0.25, "HVTN 704 T2" = 0.25, 
                    "HVTN 703 T1" = 0.25, "HVTN 704 T1" = 0.25)
strataWts_ind <- c("HVTN 703" = 0.5, "HVTN 704" = 0.5)

# get pooled and individual cumulative incidence estimates
# no censoring, no truncation
mitt.cuminc_pool <- naCumInc(
  data = mitt,
  futimeVar = timeVar, 
  eventVar = eventIndVar, 
  groupVar = grpVar_pool, 
  strataVar = strataVar_pool,
  idVar = idVar,
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_pool
)

mitt.cuminc_ind <- naCumInc(
  data = mitt,
  futimeVar = timeVar, 
  eventVar = eventIndVar, 
  groupVar = grpVar_ind, 
  strataVar = strataVar_ind,
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_ind
)

# censor at tau, only used to find tau
mitt.cuminc_pool_tau <- naCumInc(
  data = mitt,
  futimeVar = timeVar, 
  eventVar = eventIndVar, 
  groupVar = grpVar_pool, 
  strataVar = strataVar_pool,
  idVar = idVar,
  censor = list( minAtRisk=150),
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_pool
)

mitt.cuminc_ind_tau <- naCumInc(
  data = mitt,
  futimeVar = timeVar, 
  eventVar = eventIndVar, 
  groupVar = grpVar_ind, 
  censor = list(minAtRisk=150),
  strataVar = strataVar_ind,
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_ind
)

tau = mitt.cuminc_pool_tau$censorInfo$censorTime
stopifnot(tau==mitt.cuminc_ind_tau$censorInfo$censorTime)
mitt[[timeVar_trunc]] = pmin(tau, mitt[[timeVar]])

# compute for truncated follow_up times
mitt.cuminc_pool_trunc <- naCumInc(
  data = mitt,
  futimeVar = timeVar_trunc, 
  eventVar = eventIndVar, 
  groupVar = grpVar_pool, 
  strataVar = strataVar_pool,
  idVar = idVar,
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_pool
)

mitt.cuminc_ind_trunc <- naCumInc(
  data = mitt,
  futimeVar = timeVar_trunc, 
  eventVar = eventIndVar, 
  groupVar = grpVar_ind, 
  strataVar = strataVar_ind,
  stratifiedEstimator = TRUE,
  strataWeights = strataWts_ind
)


# get pooled and individual CIR estimates first based on uncensored and untruncated data
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


# get pooled and individual CIR/PE estimates first based on truncated data
mitt.CIR_pool_trunc <- EffCIR( mitt.cuminc_pool_trunc, refLvl = refLvl, cmpLvl=cmpLvl_pool, nullHypEff=0)

  mitt.CIR_pool_trunc$tests$cmpLvl <- mitt.CIR_pool_trunc$eff$cmpLvl
  mitt.CIR_pool_trunc$tests$refLvl <- mitt.CIR_pool_trunc$eff$refLvl
  mitt.CIR_pool_trunc$tests$comparison <- paste0(mitt.CIR_pool_trunc$eff$cmpLvl, " vs. ", mitt.CIR_pool_trunc$eff$refLvl)
  mitt.CIR_pool_trunc$tests <- mitt.CIR_pool_trunc$tests[,c(7,8,9,1,2,3,4,5,6)]

mitt.CIR_ind_trunc  <- EffCIR( mitt.cuminc_ind_trunc,  refLvl = refLvl, cmpLvl=cmpLvl_ind, nullHypEff=0)

  mitt.CIR_ind_trunc$tests$cmpLvl <- mitt.CIR_ind_trunc$eff$cmpLvl
  mitt.CIR_ind_trunc$tests$refLvl <- mitt.CIR_ind_trunc$eff$refLvl
  mitt.CIR_ind_trunc$tests$comparison <- paste0(mitt.CIR_ind_trunc$eff$cmpLvl, " vs. ", mitt.CIR_ind_trunc$eff$refLvl)
  mitt.CIR_ind_trunc$tests <- mitt.CIR_ind_trunc$tests[,c(7,8,9,1,2,3,4,5,6)]


# output CIR, these files estimate throught he last infection date
write.csv( mitt.CIR_pool$CIR, file=CIR.csvFile_pool, na="", row.names=FALSE, quote=FALSE)
write.csv( mitt.CIR_ind$CIR,  file=CIR.csvFile_ind, na="", row.names=FALSE, quote=FALSE)

q(save = "no")
