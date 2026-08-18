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
dataFile <- file.path(datDir, "amp_survival.csv")

# output file names
CIR.csvFile_pool  <- file.path(outDatDir, "amp_cir_wk104_pool.csv")

# specify variables names from input dataset

# groupings to be compared
grpVar_pool  <- "rx_pool"
grpVar_ind <- "rx_code"

# follow-up time information
timeVar <- "fudayswk104"

# event indicator
eventIndVar <- "statuswk104"

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


# cumulative incidence uncensored
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


# get pooled and individual CIR/PE estimates based on truncated data
mitt.CIR_pool <- EffCIR( mitt.cuminc_pool, refLvl = refLvl, cmpLvl=cmpLvl_pool, nullHypEff=0)

  mitt.CIR_pool$tests$cmpLvl <- mitt.CIR_pool$eff$cmpLvl
  mitt.CIR_pool$tests$refLvl <- mitt.CIR_pool$eff$refLvl
  mitt.CIR_pool$tests$comparison <- paste0(mitt.CIR_pool$eff$cmpLvl, " vs. ", mitt.CIR_pool$eff$refLvl)
  mitt.CIR_pool$tests <- mitt.CIR_pool$tests[,c(7,8,9,1,2,3,4,5,6)]

# output CIR on uncensored data
write.csv( mitt.CIR_pool$CIR, file=CIR.csvFile_pool, na="", row.names=FALSE, quote=FALSE)

q(save = "no")
