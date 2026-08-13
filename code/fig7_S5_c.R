# PE by PNGS 230

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")


library(survival)
library(tidyverse)

source(file.path(repoDir, "code/macro/DPEutils.R"))
source(file.path(repoDir, "code/macro/lunnMcneil.R"))
source(file.path(repoDir,"code/macro/forest.R"))


dataPostWk80 <-  read.csv(file.path(datDir, "amp_sieve_marks_wk80to104.csv"))
dataPostWk80.trial <- trial_dose_data (dataPostWk80, "704and703", dose = "T1+T2")

dataWk80 <-  read.csv(file.path(datDir, "amp_sieve_marks_wk80.csv"))
dataWk80.trial <- trial_dose_data (dataWk80, "704and703", dose = "T1+T2")

features.to.analyze.x <- "hxb2.230.pngs.ls"


doseDesc <- "Dose-pooled"
trialDesc <- "Pooled AMP Trials"

# Grab the column containing feature values from the master.trial table
.feature <- features.to.analyze.x
.xPostWk80 <- dataPostWk80.trial[ , .feature]
.xWk80 <- dataWk80.trial[ , .feature]

### Run DVE
.dPE.resultsPostWk80 <- DVE.f( dataPostWk80.trial$hiv1fpday, dataPostWk80.trial$hiv1event, .xPostWk80, 
                               dataPostWk80.trial$txLabel, dataPostWk80.trial$stratVar )
.dPE.resultsWk80 <- DVE.f( dataWk80.trial$hiv1fpday, dataWk80.trial$hiv1event, .xWk80, 
                               dataWk80.trial$txLabel, dataWk80.trial$stratVar )

### List of columns in output 
results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                       "DPE.p.value", 
                       paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                       paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;

results <- data.frame( feature = features.to.analyze.x, feature.name = features.to.analyze.x)

results[ , 3:( length( results.colnames) ) ] <- NA;
colnames(results) <- results.colnames

results[1,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                             paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                             paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.resultsWk80 , 
                                                                                                                                    "T1+T2", .xWk80,
                                                                                                                                    dataWk80.trial)


results[2,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
               paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
               paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.resultsPostWk80 , 
                                                                                                                      "T1+T2", .xPostWk80,
                                                                                                                       dataPostWk80.trial)
#Create a table for forest 
PEtable <- tibble("cohort" = character(), "Haplotype" = character(), "n" = character(),
                  "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                  "P" = character(), "DiffP" = character())


result.f <- results[1, ]
PEtable <- add_row(.data = PEtable ,"cohort" = "Week 80", "Haplotype" = "", "n" = "",
                   "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                   "DiffP" = format.p(result.f$DPE.p.value, 2))          
Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( results$mark.name[1], split = "\\.")[[1]][1:3], collapse="."))
PEtable <- add_row(.data = PEtable, cohort = "", "Haplotype" = Haplotype.x$Haplotype1, 
                   "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                   "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                   "mean" = as.numeric(result.f$PE.type.1.estimate),
                   "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                   "upper" = as.numeric(result.f$PE.type.1.CI.high),
                   "P" = format.p(result.f$PE.type.1.p.value,2),
                   "DiffP" = "")
PEtable <- add_row(.data = PEtable, cohort = "", "Haplotype" = Haplotype.x$Haplotype0, 
                   "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                   "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                   "mean" = as.numeric(result.f$PE.type.0.estimate),
                   "lower" = as.numeric(result.f$PE.type.0.CI.low),
                   "upper" = as.numeric(result.f$PE.type.0.CI.high),
                   "P" = format.p(result.f$PE.type.0.p.value,2),
                   "DiffP" = "")


result.f <- results[2, ]
PEtable <- add_row(.data = PEtable ,"cohort" = "Post-Week 80", "Haplotype" = "", "n" = "",
                   "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                   "DiffP" = format.p(result.f$DPE.p.value, 2)) 
Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( results$mark.name[1], split = "\\.")[[1]][1:3], collapse="."))
PEtable <- add_row(.data = PEtable, cohort = "", "Haplotype" = Haplotype.x$Haplotype1, 
                   "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                   "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                   "mean" = as.numeric(result.f$PE.type.1.estimate),
                   "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                   "upper" = as.numeric(result.f$PE.type.1.CI.high),
                   "P" = format.p(result.f$PE.type.1.p.value,2),
                   "DiffP" = "")
PEtable <- add_row(.data = PEtable, cohort = "", "Haplotype" = Haplotype.x$Haplotype0, 
                   "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                   "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                   "mean" = as.numeric(result.f$PE.type.0.estimate),
                   "lower" = as.numeric(result.f$PE.type.0.CI.low),
                   "upper" = as.numeric(result.f$PE.type.0.CI.high),
                   "P" = format.p(result.f$PE.type.0.p.value,2),
                   "DiffP" = "")




PEtable.forestplot.withSieveT5( PEtable,
                                xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                figDir, paste0("703and704_sievePH_PEbyglycosite230PNGS_ls.pdf"),
                                header = c("  ","Env\n230-232","No. Cases (VRC01 vs. P)\n(Incidence per 100 PYRs)",
                                           "PE (%) (95% CI)","mean", "lower", "upper","P-value","Diff PE\nP-value"))




#look IC80 for 230 PNG
df <- filter(dataPostWk80.trial, hiv1event==1 & hxb2.230.pngs.ls == 1)
df <- dplyr::select(df, all_of(c("protocol", "tx", "gmt80ls")))
write.csv(df, file.path(datDir, "ic80_PNG230_postWk80.csv"))

#plot IC80 by PNGS status, cohort, and treatment
dfWk80 <- dplyr::select(filter(dataWk80, hiv1event==1), all_of(c("protocol", "tx_pool", "gmt80ls", "hxb2.230.pngs.ls")))
dfWk80$gmt80ls[dfWk80$gmt80ls == ">100"] <- "100"
dfWk80$gmt80ls <- as.numeric(dfWk80$gmt80ls)
dfWk80$cohort <- "Weeks 0 to 80"

dfPostWk80 <- dplyr::select(filter(dataPostWk80, hiv1event==1), all_of(c("protocol", "tx_pool", "gmt80ls", "hxb2.230.pngs.ls")))
dfPostWk80$gmt80ls[dfPostWk80$gmt80ls == ">100"] <- "100"
dfPostWk80$gmt80ls <- as.numeric(dfPostWk80$gmt80ls)
dfPostWk80$cohort <- "Weeks 80 to 104"

df <- rbind(dfWk80, dfPostWk80)
df$tx <- ifelse(df$tx_pool=="T1+T2", "VRC01\nDose-pooled", "Placebo")
df$gmt80ls <- log10(df$gmt80ls)
df$group <- ifelse(df$hxb2.230.pngs.ls==1, "PNGS at 230-232", "No PNGS at 230-232")
df$group <- factor(df$group, levels = c("PNGS at 230-232", "No PNGS at 230-232"))

set.seed(8)
p <- ggplot(filter(df, !is.na(group))) +
     geom_violin(aes(x = tx, y = gmt80ls, colour = tx)) +
     geom_boxplot(aes(x = tx, y = gmt80ls, fill = tx, colour = tx), width=0.15, lwd=1.25, alpha = 0.3, outlier.shape=NA)+
     geom_jitter(aes(x = tx, y = gmt80ls, color = tx), width = 0.3, height = 0, size = 2)+
     facet_grid(cols = vars(group), rows = vars(cohort))+
     scale_colour_manual(breaks = c("Placebo", "VRC01\nDose-pooled"), values = c("blue", "red3"))+
     scale_fill_manual(breaks = c("Placebo", "VRC01\nDose-pooled"), values = c("blue", "red3"))+
     scale_y_continuous(breaks = log10(c(0.1, 0.3, 1, 3, 10, 30, 100)), labels = c("0.1", "0.3", "1", "3", "10", "30", expression("">="100")))+
     theme_bw() +
     ylab(expression("Measured IC80 ("*mu*"g/ml)"))+
     xlab("")+
     theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=23),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r=0, b = 0, l = -2)))


ggsave(file.path(figDir, paste0("IC80byPNGS230.pdf")),
       plot=p, width=8, height=8)

#number of participants with observed IC80
library(plyr)
ans <- ddply(df, .(cohort, tx), function(x) c(dim(x)[1], sum(!is.na(x$gmt80ls)), 
                                              sum(!is.na(x$hxb2.230.pngs.ls)),
                                              sum(!is.na(x$gmt80ls) & !is.na(x$hxb2.230.pngs.ls))))
colnames(ans) <- c("cohort", "treatment", "No. of Cases", "No. of Cases with Measured IC80", "No. of Cases with ls Sequence",
                   "No. of Cases with Measured IC80 and ls sequence")

write.csv(ans, file.path(tabDir, "IC80_pngs230_count.csv"), row.names = FALSE)


#plot IC50 by PNGS status, cohort, and treatment
dfWk80 <- dplyr::select(filter(dataWk80, hiv1event==1), all_of(c("protocol", "tx_pool", "gmt50ls", "hxb2.230.pngs.ls")))
dfWk80$gmt50ls[dfWk80$gmt50ls == ">100"] <- "100"
dfWk80$gmt50ls <- as.numeric(dfWk80$gmt50ls)
dfWk80$cohort <- "Weeks 0 to 80"

dfPostWk80 <- dplyr::select(filter(dataPostWk80, hiv1event==1), all_of(c("protocol", "tx_pool", "gmt50ls", "hxb2.230.pngs.ls")))
dfPostWk80$gmt50ls[dfPostWk80$gmt50ls == ">100"] <- "100"
dfPostWk80$gmt50ls <- as.numeric(dfPostWk80$gmt50ls)
dfPostWk80$cohort <- "Weeks 80 to 104"

df <- rbind(dfWk80, dfPostWk80)
df$tx <- ifelse(df$tx_pool=="T1+T2", "VRC01\nDose-pooled", "Placebo")
df$gmt50ls <- log10(df$gmt50ls)
df$group <- ifelse(df$hxb2.230.pngs.ls==1, "PNGS at 230-232", "No PNGS at 230-232")
df$group <- factor(df$group, levels = c("PNGS at 230-232", "No PNGS at 230-232"))

set.seed(8)
p <- ggplot(filter(df, !is.na(group))) +
  geom_violin(aes(x = tx, y = gmt50ls, colour = tx)) +
  geom_boxplot(aes(x = tx, y = gmt50ls, fill = tx, colour = tx), width=0.15, lwd=1.25, alpha = 0.3, outlier.shape=NA)+
  geom_jitter(aes(x = tx, y = gmt50ls, color = tx), width = 0.3, height = 0, size = 2)+
  facet_grid(cols = vars(group), rows = vars(cohort))+
  scale_colour_manual(breaks = c("Placebo", "VRC01\nDose-pooled"), values = c("blue", "red3"))+
  scale_fill_manual(breaks = c("Placebo", "VRC01\nDose-pooled"), values = c("blue", "red3"))+
  scale_y_continuous(breaks = log10(c(0.1, 0.3, 1, 3, 10, 30, 100)), labels = c("0.1", "0.3", "1", "3", "10", "30", expression("">="100")))+
  theme_bw() +
  ylab(expression("Measured IC50 ("*mu*"g/ml)"))+
  xlab("")+
  theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=23),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r=0, b = 0, l = -2)))


ggsave(file.path(figDir, paste0("IC50byPNGS230.pdf")),
       plot=p, width=8, height=8)
