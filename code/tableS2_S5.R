# Purpose: linear regression for continuous variables including group (Placebo vs. VRC01 pooled) 
#           adjusting for the 4-level region factor (RSA vs. not-RSA in Africa and US+Switzarland vs. 
#           South America for the Americas trial).  
#           For ordinal marks run a proportional odds model again adjusting for region and group.  
#           Report the beta, 95% CI, and p-vlaue for group.  
# Author:  Li Li
# Date:    Oct 29, 2025

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")

library(gridExtra)
library(ggpubr)
library(tidyverse)
library(plyr)
library(MASS)

logit <- function(p){
  return(log(p/(1-p)))
}

dataWk80 <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9c.csv")) %>%
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
dataWk80$epitope.dist.subtype.ls <- ifelse(dataWk80$protocol == "HVTN 703", dataWk80$epitope.dist.c.ls, 
                                               dataWk80$epitope.dist.b.ls)

#week 80 - week 104
dataPostWk80 <- read.csv(file.path(datDir, "d_wk80_wk104_survival_dataset_sieve.csv")) %>%
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
#log10 transformation 
dataPostWk80$gmt80ls <- log10(dataPostWk80$gmt80ls)
dataPostWk80$epitope.dist.subtype.ls <- ifelse(dataPostWk80$protocol == "HVTN 703", dataPostWk80$epitope.dist.c.ls, 
                                               dataPostWk80$epitope.dist.b.ls)


dataWk80_cases <- filter(dataWk80, hiv1event == 1)
dataPostWk80_cases <- filter(dataPostWk80, hiv1event == 1)

#number of PNGs in v5
mean(dataPostWk80_cases$num.pngs.v5.ls[dataPostWk80_cases$tx_pool=="T1+T2"])
mean(dataPostWk80_cases$num.pngs.v5.ls[dataPostWk80_cases$tx_pool=="C3"])


escapeDesc <- read.csv(file.path(datDir, "VRC01escape_DescFile.csv"))
escapeDesc <- filter(escapeDesc, pub_id%in% c(dataWk80_cases$pub_id, dataPostWk80_cases$pub_id))
escapeDesc$cohort <- ifelse(escapeDesc$pub_id %in% dataWk80_cases$pub_id, "Weeks 0 to 80", "Weeks 80 to 104")

escapeFeatures <- c(colnames(escapeDesc)[grepl("score", colnames(escapeDesc))], "resmark_comb_tbl2")
escapeFeatures_desc <- c("60","197","234","276","278","279","280","281","317","364","365","369","371","372","425","456",         
                   "459","460","461","463","465","471","474","616","618",
                   "Length of V1V2","No. of PNGs in V1V2","Length of V5","No. of PNGs in V5", "resmark_comb_tbl2")
names(escapeFeatures_desc) <- escapeFeatures


dataWk80_cases <- left_join(dataWk80_cases, dplyr::select(escapeDesc, all_of(c("pub_id", escapeFeatures))), by = "pub_id")
dataPostWk80_cases <- left_join(dataPostWk80_cases, dplyr::select(escapeDesc, all_of(c("pub_id", escapeFeatures))), by = "pub_id")

#Linear regression for continuous features
quantMarks <- c("parscore1.ls","parscore2.ls", "gmt80ls", "epitope.dist.subtype.ls",  
                "hdist.zspace.sites.binding.all.ls", "resmark_comb_tbl2")

quantMarkslabels <- c("Predicted Odds of IC80 > 1", "Predicted IC80","IC80", "VRC01 Epitope Distance to Subtype-Specific Reference", 
                      "PC-Weighted Hamming Distance in 50 VRC01 or CD4 Binding Positions","VRC01 Resistance Feature Score")
names(quantMarkslabels) <- quantMarks

results_q <- tibble("cohort" = character(), "feature" = character(), "VRC01" = character(), "Placebo" = character(), "beta" = character(), "95% CI" = character(), 
                  "P-value" = numeric(), "P-value2" = character())

results_e <- tibble("cohort" = character(), "feature" = character(), "beta" = character(), "95% CI" = character(), 
                    "P-value" = numeric(), "P-value2" = character(), "n_vaccine(-1)" = character(), "n_vaccine(0)" = character(), "n_vaccine(1)" = character(),
                    "n_placebo(-1)" = character(), "n_placebo(0)" = character(), "n_placebo(1)" = character())

for(cohortx in c("Weeks 0 to 80", "Weeks 80 to 104")){
  if(cohortx == "Weeks 0 to 80"){
    df <- dplyr::select(dataWk80_cases, all_of(c("tx_pool", "stratVar", quantMarks, escapeFeatures)))
  }else{
    df <- dplyr::select(dataPostWk80_cases, all_of(c("tx_pool", "stratVar", quantMarks, escapeFeatures)))
  }
  for(f in quantMarks){
    df$mark <- df[, f]
    fit <- lm(mark ~ tx_pool + stratVar, data = df)
    sfit <- summary(fit)
    betaest <- sfit$coefficients["tx_poolT1+T2","Estimate"]
    betastd <- sfit$coefficients["tx_poolT1+T2","Std. Error"]
    ct <- qt(0.975, df = sfit$df[2])
    ul <- betaest + ct*betastd
    ll <- betaest - ct*betastd
    fit_reduced <- lm(mark ~ stratVar, data = df)
    lr <- anova(fit, fit_reduced)
    p <- lr$`Pr(>F)`[2]
    p2 <- ifelse(p<0.001, "< 0.001", paste0(format(p, digits=2, nsmall=2)))
    
    df2 <- filter(df, !is.na(mark))
    coefest <- sfit$coefficients[,"Estimate"]
    design_matrix <- model.matrix(~ stratVar, data = df2)
    newx1 <- cbind(rep(1, dim(design_matrix)[1]),rep(1, dim(design_matrix)[1]), design_matrix[, 2:4])
    newx0 <- cbind(rep(1, dim(design_matrix)[1]),rep(0, dim(design_matrix)[1]), design_matrix[, 2:4])
    
    mean_vaccine <- mean(newx1%*%coefest, na.rm = TRUE)
    mean_placebo <- mean(newx0%*%coefest, na.rm = TRUE)
    
    #mean_vaccine <- mean(df$mark[df$tx_pool=="T1+T2"], na.rm = TRUE)
    #mean_placebo <- mean(df$mark[df$tx_pool=="C3"], na.rm = TRUE)
    
    if(f %in% c("parscore2.ls", "gmt80ls")){
      results_q <- add_row(.data = results_q, "cohort" = cohortx, "feature" = quantMarkslabels[f],
                           "VRC01" = format(10^mean_vaccine, digits = 2, nsmall = 2), "Placebo" = format(10^mean_placebo, digits = 2, nsmall = 2),
                           "beta" = format(10^betaest, digits=2, nsmall=2),
                           "95% CI" = paste0("(", format(10^ll, digits = 2, nsmall = 2), ", ", 
                                             format(10^ul, digits = 2, nsmall = 2), ")"),
                           "P-value" = p, "P-value2" = p2)
    }else if (f=="parscore1.ls"){
      results_q <- add_row(.data = results_q, "cohort" = cohortx, "feature" = quantMarkslabels[f],
                           "VRC01" = format(exp(mean_vaccine), digits = 2, nsmall = 2), "Placebo" = format(exp(mean_placebo), digits = 2, nsmall = 2),
                           "beta" = format(exp(betaest), digits=2, nsmall=2),
                           "95% CI" = paste0("(", format(exp(ll), digits = 2, nsmall = 2), ", ", 
                                             format(exp(ul), digits = 2, nsmall = 2), ")"),
                           "P-value" = p, "P-value2" = p2)
    }else{
      results_q <- add_row(.data = results_q, "cohort" = cohortx, "feature" = quantMarkslabels[f],
                           "VRC01" = format(mean_vaccine, digits = 2, nsmall = 2), "Placebo" = format(mean_placebo, digits = 2, nsmall = 2),
                           "beta" = format(betaest, digits=2, nsmall=2),
                           "95% CI" = paste0("(", format(ll, digits = 2, nsmall = 2), ", ", 
                                             format(ul, digits = 2, nsmall = 2), ")"),
                           "P-value" = p, "P-value2" = p2)
    }
      
      
    }
    
  }


for(cohortx in c("Weeks 0 to 80", "Weeks 80 to 104")){
  if(cohortx == "Weeks 0 to 80"){
    df <- dplyr::select(dataWk80_cases, all_of(c("tx_pool", "stratVar", quantMarks, escapeFeatures)))
  }else{
    df <- dplyr::select(dataPostWk80_cases, all_of(c("tx_pool", "stratVar", quantMarks, escapeFeatures)))
  }
  for(f in escapeFeatures){

    mark <- as.factor(df[, f])
    tx_pool <- df$tx_pool
    startVar <- df$stratVar
    frequencies <- table(mark)
    
    df0 <- data.frame(tx_pool =  tx_pool, mark = mark)
    frequenciesbyTx <- ddply(df0, .(tx_pool), function(dfx){
      nt <- length(dfx$mark) 
      ans <- c(sum(dfx$mark == "-1"), sum(dfx$mark == "0"), sum(dfx$mark == "1"))
      p <- round(ans/nt*100, 1)
      ans2<- paste(ans, " (", p,"%", ")", sep="")
      return(ans2)
    })
    colnames(frequenciesbyTx) <- c("treatment","-1", "0", "1")
    rownames(frequenciesbyTx) <- c("C3", "T1+T2")
    
    
    if(length(names(frequencies)) == 1){
      results_e <- add_row(.data = results_e, "cohort" = cohortx, "feature" = f,
                         "beta" = NA,"95% CI" = NA,"P-value" = NA, "P-value2" = NA,
                         "n_vaccine(-1)" = frequenciesbyTx["T1+T2","-1"], "n_vaccine(0)" = frequenciesbyTx["T1+T2","0"], "n_vaccine(1)" = frequenciesbyTx["T1+T2","1"],
                         "n_placebo(-1)" = frequenciesbyTx["C3","-1"], "n_placebo(0)" = frequenciesbyTx["C3","0"], "n_placebo(1)" = frequenciesbyTx["C3","1"])
    }else if (length(names(frequencies)) == 2){
      if(min(frequencies) <5){
        results_e <- add_row(.data = results_e, "cohort" = cohortx, "feature" = f,
                           "beta" = NA,"95% CI" = NA,"P-value" = NA, "P-value2" = NA,
                           "n_vaccine(-1)" = frequenciesbyTx["T1+T2","-1"], "n_vaccine(0)" = frequenciesbyTx["T1+T2","0"], "n_vaccine(1)" = frequenciesbyTx["T1+T2","1"],
                           "n_placebo(-1)" = frequenciesbyTx["C3","-1"], "n_placebo(0)" = frequenciesbyTx["C3","0"], "n_placebo(1)" = frequenciesbyTx["C3","1"])
      }else{
        fit <- tryCatch(glm(mark~ tx_pool + startVar, family = binomial(link = "logit" )))
        sfit <- summary(fit)
        betaest <- sfit$coefficients["tx_poolT1+T2","Estimate"]
        betastd <- sfit$coefficients["tx_poolT1+T2","Std. Error"]
        ct <- qt(0.975, df = fit$df.residual)
        ul <- betaest + ct*betastd
        ll <- betaest - ct*betastd
        
        fit_reduced <- glm(mark~ startVar, family = binomial(link = "logit" ))
        lr <- anova(fit_reduced, fit)
        p <- lr$`Pr(>Chi)`[2]
        
        p2 <- ifelse(p<0.001, "< 0.001", paste0(format(p, digits=2, nsmall=2)))
        results_e <- add_row(.data = results_e, "cohort" = cohortx, "feature" = f,
                           "beta" = format(exp(betaest), digits=2, nsmall=2),
                           "95% CI" = paste0("(", format(exp(ll), digits = 2, nsmall = 2), ", ", 
                                             format(exp(ul), digits = 2, nsmall = 2), ")"),
                           "P-value" = p, "P-value2" = p2,
                           "n_vaccine(-1)" = frequenciesbyTx["T1+T2","-1"], "n_vaccine(0)" = frequenciesbyTx["T1+T2","0"], "n_vaccine(1)" = frequenciesbyTx["T1+T2","1"],
                           "n_placebo(-1)" = frequenciesbyTx["C3","-1"], "n_placebo(0)" = frequenciesbyTx["C3","0"], "n_placebo(1)" = frequenciesbyTx["C3","1"])
        
      }
      }else{
      if(sum(frequencies>=5) <2){
        results_e <- add_row(.data = results_e, "cohort" = cohortx, "feature" = f,
                           "beta" = NA,"95% CI" = NA,"P-value" = NA, "P-value2" = NA,
                           "n_vaccine(-1)" = frequenciesbyTx["T1+T2","-1"], "n_vaccine(0)" = frequenciesbyTx["T1+T2","0"], "n_vaccine(1)" = frequenciesbyTx["T1+T2","1"],
                           "n_placebo(-1)" = frequenciesbyTx["C3","-1"], "n_placebo(0)" = frequenciesbyTx["C3","0"], "n_placebo(1)" = frequenciesbyTx["C3","1"])
        
      }else{
        fit <- tryCatch(polr(mark~ tx_pool + startVar, method = "logistic" ))
        sfit <- summary(fit)
        betaest <- sfit$coefficients["tx_poolT1+T2","Value"]
        betastd <- sfit$coefficients["tx_poolT1+T2","Std. Error"]
        ct <- qt(0.975, df = fit$df.residual)
        ul <- betaest + ct*betastd
        ll <- betaest - ct*betastd
        
        fit_reduced <- polr(mark~ startVar, method = "logistic" )
        lr <- anova(fit, fit_reduced)
        p <- lr$`Pr(Chi)`[2]
        #p <- 2*(pt(-abs(sfit$coefficients["tx_poolT1+T2","t value"]), df = fit$df.residual))
        p2 <- ifelse(p<0.001, "< 0.001", paste0(format(p, digits=2, nsmall=2)))
        results_e <- add_row(.data = results_e, "cohort" = cohortx, "feature" = f,
                           "beta" = format(exp(betaest), digits=2, nsmall=2),
                           "95% CI" = paste0("(", format(exp(ll), digits = 2, nsmall = 2), ", ", 
                                             format(exp(ul), digits = 2, nsmall = 2), ")"),
                           "P-value" = p, "P-value2" = p2,
                           "n_vaccine(-1)" = frequenciesbyTx["T1+T2","-1"], "n_vaccine(0)" = frequenciesbyTx["T1+T2","0"], "n_vaccine(1)" = frequenciesbyTx["T1+T2","1"],
                           "n_placebo(-1)" = frequenciesbyTx["C3","-1"], "n_placebo(0)" = frequenciesbyTx["C3","0"], "n_placebo(1)" = frequenciesbyTx["C3","1"])
        
      }  
      }
    
  }
}


write.csv(results_e, file.path(tabDir, "vrc01_vs_placebo_escape_feature_regression.csv"), row.names = FALSE)
write.csv(results_q, file.path(tabDir, "vrc01_vs_placebo_quantitative_feature_regression.csv"), row.names = FALSE)

