# ***************************************************************************************
# Code History 
# Document date of changes and person making edits
# Provide a summary of edits made - both what was changed and why): 
# ----------------------------------------------------------------
#   Original verision created (???)
#
#   < many other modifications > 

#   2019Dec31  Doug Grove  
#     - Modified code to allow multiple levels to be specified in the 'cmpLvl' argument.
#       Allows a user to compute HRs for multi-arm trials within the same cox model.
# 
#     - Modified argument 'cmpLvl' to have default value of NULL, and if the user does
#       not specify the argument (so it remains NULL), the function will interpret this
#       to mean that all levels of groupVar (other than refLvl) should be compared to
#       refLvl (the reference level)
#
#     - Modified the "eff" component of the output object, to include three additional 
#       columns: 'cmpLvl', 'refLvl' and 'comparison', which provide info on which 
#       comparisons the estimates in a row correspond to.  This is required when there 
#       are multiple comparisons being done and is useful in all cases.
#

phReg <-
  function(data,  futimeVar, eventVar,  groupVar,  refLvl, cmpLvl=NULL,
           strataVar= NULL, confLvl=0.95, nullHypEff=NULL, cox.zph=FALSE ) 
{

  # NOTES:
  # ------
  #  1. This is not a general purpose coxph alternative, it's a wrapper that allows 
  #     similar info to be obtained from the cox model as from our use of the 
  #     cumulative-incidence-ratio.  The usage allowed here is fairly restrictive,
  #     please only use it for it's intended purpose.
  #
  #  2. if the data.frame 'data' contains more than two groups (specified via
  #     the values of 'groupVar'), then the data will be subset down to just those
  #     two groups before fitting the cox model

  # Arguments
  # ---------
  # data:  a data.frame containing the variables specified by arguments 
  #        'futimeVar', 'eventVar', 'groupVar' and 'strataVar' (if given)
  #
  # futimeVar: a character string giving the name of a variable in 'data' containing
  #            the follow-up times
  #
  # eventVar: a character string giving the name of a variable in 'data' containing
  #           an integer vector indicating the occurence of an event (1=yes, 0=no)
  #
  # groupVar:  a character string giving the name of a variable in 'data'.
  #     The variable must be categorical and should identify a characteristic of interest.
  #     Typically this would be something like treatment group, but could be anything.
  #
  # refLvl: the "level"/category of the group variable to be used as the reference-level
  #         (i.e. denominator) when computing the hazard ratio.
  #         This argument is required
  #
  # cmpLvl: the level(s) of the group variable to be compared to the reference.
  #         It is used as the numerator of the hazard ratio
  #         If cmpLvl is left as 'NULL', then *all* non-reference levels of groupVar
  #         will be used 
  #
  # strataVar:  an optional character *string* giving the name of the variable in 'data'
  #    variable in 'data' to be stratified on.  Stratification in a cox model implies 
  #    allowing different levels to have different baseline hazards.  If you have multiple
  #    variables you want to stratify on (e.g. sex and dosage) you can do this by
  #    creating a single variable that jointly specifies the levels of all those variables
  #    and use that. You can create such a variable using e.g. the function interaction().
  #
  # confLevel: the confidence level used to constructing the CI for the hazard-ratio
  #
  # nullHypEff:  The value of the null hypotheses that the user wishes to test.  It must
  #    be a numeric vector with values < 1, more than one hypothesis can be specified.
  #    The values provided should specify the null for the *efficacy* parameter (i.e. 
  #    1 - hazard-ratio), not for e.g. the log-hazard-ratio. 
  #    If the argument is not specified by the user, the function returns results for the
  #    nullHypEff=0 (Efficacy = 0 , Hazard-ratio = 1) 
  #    Results are returned in the 'tests' component of the output list 
  #
  #    All tests are two-sided tests of
  #       H0: Efficacy=xxx vs. HA: Efficacy != xxx,
  #    they are NOT testing the directional hypothesis
  #       H0: Efficacy <= xxx  vs. HA: Efficacy > xxx
  #
  # cox.zph: logical values indicating whether or not to run cox.zph to test the
  #    assumption of proportional hazards.  Default value is FALSE.
  #    Results are returned in the 'zph' component of the output list 
  #  

  ## These first 4 checks are done by R for us
  #if ( missing(data) )
  #  stop("You have not provided a value for the argument 'data'\n\n")

  #if ( missing(futimeVar) )
  #  stop("You have not provided a value for the argument 'futimeVar'\n\n")

  #if ( missing(eventVar) )
  #  stop("You have not provided a value for the argument 'eventVar'\n\n")

  #if ( missing(groupVar) )
  #  stop("You have not provided a value for the argument 'groupVar'\n\n")

  ## Do basic check on nullHypEff if it was specified
  if ( !is.null(nullHypEff) ) {
    if (!is.numeric(nullHypEff) || any(nullHypEff >= 1) )
      stop("Argument 'nullHypEff' must be a numeric vector with all values < 1 \n\n")
  }

  Strata <- ( !is.null(strataVar) )


  ## make sure variables specified are in the data.frame provided
  vars<- c(futimeVar, eventVar, groupVar, if (Strata) strataVar)
  dataVars <- names(data)
  for (v in vars) {
    if ( !v %in% dataVars ) 
      stop("The variable ", v, " is not present in data.frame '", 
           as.list(match.call())$data, "' specified by argument 'data'\n\n")
  }

  grpLvls <- unique( data[[groupVar]] )

  ## check that 'refLvl' is among the group levels found in the input object
  if ( !(refLvl %in% grpLvls) )
    stop("The value specified for 'refLvl' does not match any values of variable ", 
         groupVar, "\n", "in the data.frame provided via argument 'data'\n\n")

  ## if cmpLvl is specified, check that the value(s) specified are among the group levels
  if ( !is.null( cmpLvl ) ) {
    notFound <- !( cmpLvl %in% grpLvls )
    if ( any(notFound) )
      stop("The following levels specified via argument 'cmpLvl' were not found in the input:\n",
            paste( cmpLvl[ notFound ], collapse=","), ".\n\n" )
  } else {
    ## if cmpLvl wasn't specified, assign it's default value, i.e.  all non-reference levels
    cmpLvl <- grpLvls[ grpLvls != refLvl ]
  }


  ## subset data to the groups specified by cmpLvl and refLvl and
  ## to the variables we need to fit the model
  d <- data[ data[[groupVar]] %in% c(cmpLvl, refLvl), 
             c(futimeVar, eventVar, groupVar, strataVar) ]


  ## convert the group variable into a factor, and set refLvl to be the first level
  ## (coxph uses the first level as it's reference level for estimation)
  d[[groupVar]] <- factor( d[[groupVar]], levels=c(refLvl, cmpLvl) )


  if ( Strata ) {
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ ", groupVar,
                            " + strata(", strataVar, ")" )
  } else {
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ ", groupVar )
  }


  ## load survival package
  library(survival)

  ## coerce our character string 'formula.char' to an actual "formula"
  SurvFormula <- as.formula( formula.char )

  ## assign the default value to 'nullHypEff' so we can use one set of code to
  ## process all situations (i.e. when user has and has not specified nullHypEff)
  if (is.null(nullHypEff)) 
    nullHypEff <- 0

  ## We need to the model once for each value of vector 'nullHypEff'
  ## Results of these fits will be identical except for the results of the tests specified

  ## create a list to store output of each test into
  testList <- vector("list", length=length(nullHypEff))

  ## loop over tests,
  for (ii in 1:length(nullHypEff)) {

      hyp <- nullHypEff[ii]

      ## Note: the value of 'init' needs to have the same length as 'cmpLvl' (i.e. as the
      ##  number of parameters being estimated), hence the usage of 'rep()' below
      cobj <- coxph( SurvFormula, data=d, init = log(1 - rep(hyp,length(cmpLvl)))  )
      summ <- summary( cobj, conf.int = confLvl)

      ## bind together results from three tests done by default
      tmpTests <- as.data.frame( do.call("rbind", summ[c("waldtest","sctest","logtest")] ) )

      ## re-arrange/re-name results and supplement with name and null hypothesis being tested
      testList[[ ii ]] <-
          data.frame(
              test = c("Wald", "Score", "Likelihood Ratio"),
              nullHypEff = hyp,
              nullHypHR  = 1 - hyp,
              testStatistic = tmpTests$test,
              df = tmpTests$df,
              pvalue = tmpTests$pvalue,
              stringsAsFactors = FALSE )

  }
  rm(tmpTests)


  if ( isTRUE(cox.zph) ) {
    zph <- cox.zph( cobj )
    cat("Results of cox.zph(), test of proportionality of hazards:\n\n")
    cat("\n", "-------------------------------------------------\n\n", sep="")
    print( zph )
    cat("\n", "-------------------------------------------------\n\n", sep="")
  }

  ## extract HazardRatio and confidence Interval, then convert to VE
  ests <- summ$conf.int

  ## get the name of the group(s) used as the numerator in the hazard ratio
  grpLvl <- rownames( ests )
  if ( ! all(grpLvl %in% paste0(groupVar, cmpLvl) ) )
    stop("The name(s) of the output parameter(s) are not what is expected.\n",
         "The output names are ", grpLvl, ", but we are expecting ", 
          paste0(groupVar, cmpLvl), ".\n  Please investigate\n\n")

  ## create output object
  out <- data.frame(
             cmpLvl = cmpLvl,
             refLvl = refLvl,
             comparison = paste0( cmpLvl, " vs. ", refLvl),

             HR    = ests[, "exp(coef)"],
             lo.HR = ests[, substr(colnames(ests),1,5) == "lower"],
             up.HR = ests[, substr(colnames(ests),1,5) == "upper"], 

             eff   = NA,
             lo.eff = NA,
             up.eff = NA,
             confLevel = confLvl )

  ## compute values for eff, lo.eff and up.eff (simpler to do in a 'within' statement)
  out <- within(out, {
             eff    <- 1 - HR
             lo.eff <- 1 - up.HR
             up.eff <- 1 - lo.HR 
         })



  ## To do:  Make a proper output object
  ##    Include components for (?):
  ##    Call, ... 
  outObj <- list(Call = NULL, eff = out, tests=do.call("rbind", testList)) 

  ## if we performed test of proportionality, return primary output
  if ( isTRUE(cox.zph) ) {
    outObj$zph <- as.data.frame( zph$table )
  }

  return( outObj ) 
}

