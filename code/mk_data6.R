library(here)
here::i_am("README.md")
repoDir <- here::here()

outDatDir <- file.path(repoDir, "data")
datDir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data" # file.path(repoDir, "data")
macroDir <- file.path(repoDir, 'code/macro')

# load macros
source(file.path(macroDir, "phReg.R")) 
source(file.path(macroDir, "cuminc_functions.R"))

# input file name
dataFile <- file.path(datDir, "v703_survival.csv")

# output file names
cuminc.csvFile_pool <- file.path(outDatDir, "v703_cuminc_wk104_pool.csv")
cuminc.csvFile_ind  <- file.path(outDatDir, "v703_cuminc_wk104_ind.csv")

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

# output cuminc/cumhaz estimates uncensored
write.csv( mitt.cuminc_pool$cuminc, file=cuminc.csvFile_pool, na="", row.names=FALSE, quote=FALSE)
write.csv( mitt.cuminc_ind$cuminc, file=cuminc.csvFile_ind, na="", row.names=FALSE, quote=FALSE)

q(save = "no")
