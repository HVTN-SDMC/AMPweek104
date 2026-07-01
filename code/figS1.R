#Simulations for model A-D

rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")

source(file.path(repoDir, "code/macro/simulationUtils.R"))

library(survival)
library(cmprsk)
library(ggplot2)
library(plyr)

seed = 5265
set.seed(seed)
n = 1000000

load = FALSE
for(scenario in c("homogeneous", "moderate", "extreme")){
  for(model in c("A","B","C","D","mixCD1", "mixCD2")){
    print(c(scenario, model))
    if(load){
      df_placebo <- readRDS(file.path(outDir, paste0("simulated_data_", scenario, "_placebo", model, ".rds")))
      df_VRC01 <- readRDS(file.path(outDir, paste0("simulated_data_", scenario, "_vrc01", model, ".rds")))
    } else {
      df_placebo <- simulate_trial_exposure_mix (n = n, scenario = scenario, arm = "placebo", seed = 1)
      df_VRC01 <- simulate_trial_exposure_mix (n = 2*n, scenario = scenario, arm = "vrc01", model = model, seed = 1)
      saveRDS(df_placebo, file = file.path(outDir, paste0("simulated_data_", scenario, "_placebo", model, ".rds")))
      saveRDS(df_VRC01, file = file.path(outDir, paste0("simulated_data_", scenario, "_vrc01", model, ".rds")))
    }
    Y = c(df_placebo$Y, df_VRC01$Y)
    Tr  = c(rep(0, n), rep(1, 2*n))
    mark = c(df_placebo$D_tilde, df_VRC01$D_tilde) + 1
    C = c(df_placebo$C, df_VRC01$C)
    event = ifelse(C==1, mark, 0)
    data <- data.frame("Y" = Y, "Tr"  = factor(Tr, 0:1, labels = c("Placebo", "VRC01")), 
                       "event" = factor(event, 0:2, labels = c("censor", "sensitive", "resistant")))
    fit <- cuminc(data$Y, data$event, data$Tr)
    time <- seq(0.1, 104, 0.1)
    est <- data.frame(t(timepoints(fit, time)$est))
    PE_s <- 1 - est$VRC01.sensitive/est$Placebo.sensitive
    PE_v <- 1 - est$VRC01.resistant/est$Placebo.resistant
    
    titleLabel <- ifelse(scenario == "homogeneous", "Homogeneous Risk", ifelse(scenario == "moderate", "Intermediate Heterogeneity", 
                                                                               "Extreme Heterogeneity"))
    
    cumincdf <- data.frame("time_" = Y, "status_" = ifelse(event ==0, 0, 1), Tr  = factor(Tr, 0:1, labels = c("Placebo", "VRC01")))
    cumincp <- cumIncF(cumincdf[cumincdf$Tr == "Placebo",], timeGrid = time)
    cumincv <- cumIncF(cumincdf[cumincdf$Tr == "VRC01",], timeGrid = time)
    cumincplotdf <- data.frame(rbind(cumincp, cumincv), c(rep("Placebo", length(time)), rep("VRC01", length(time))))
    colnames(cumincplotdf) <- c("time", "cuminc", "treatment")
    tau <-104
    PEatTau <- 1 - cumIncF(cumincdf[cumincdf$Tr == "VRC01",], timeGrid = tau)$cumInc/cumIncF(cumincdf[cumincdf$Tr == "Placebo",], timeGrid = tau)$cumInc
    PEatTau <- paste0(round(PEatTau*100, 1), "%")
    PE <- 1 - cumIncF(cumincdf[cumincdf$Tr == "VRC01",], timeGrid = seq(0.1, 104, 0.1))$cumInc/cumIncF(cumincdf[cumincdf$Tr == "Placebo",], timeGrid = seq(0.1, 104, 0.1))$cumInc
    plotdf <- data.frame(t = c(time, time, time), PE = c(PE_s, PE_v, PE), 
                         type = c(rep("Sensitive", length(time)), rep("Resistant", length(time)), rep("Overall", length(time))))
    
    
    p1 <- ggplot(plotdf, aes(x = t, y = PE, color = type)) + 
      geom_line(aes(linetype = type),linewidth = 2)+
      xlab("Time in weeks") + 
      ylab("Prevention Efficacy (%)") + 
      scale_x_continuous(breaks = seq(0, 14, 1)*8, minor_breaks = NULL, limits = c(0, 104))+
      scale_y_continuous(breaks = seq(-0.4, 1, 0.2), labels = seq(-0.4, 1, 0.2)*100, minor_breaks = NULL, limits = c(-0.5, 1.1))+
      geom_hline(aes(yintercept = 0), linetype = "dashed", size  = 0.7) +
      scale_colour_manual(breaks = c("Sensitive", "Resistant", "Overall"), values = c("#0AB7C9", "sienna3", "black"), labels = c("Sensitive","Intermediate-to-Resistant",  "Overall")) +
      scale_linetype_manual(breaks = c("Sensitive", "Resistant", "Overall"), values = c("solid", "22", "11"), labels = c( "Sensitive","Intermediate-to-Resistant", "Overall")) +
      
      theme_bw()+
      ggtitle(titleLabel) +
      theme(legend.title = element_blank(),
            legend.key.width = unit(1, 'cm'),
            legend.position = c(0.3, 0.2),
            legend.direction = "vertical",
            legend.background = element_blank(),
            legend.text = element_text (size = 14),
            title = element_text(size = 14),
            strip.text.y =  element_text (size = 16),
            strip.text.x =  element_text (size = 16),
            axis.title.x =  element_text (size = 16),
            axis.title.y =  element_text (size = 16),
            axis.text.x =  element_text (size = 16),
            axis.text.y =  element_text (size = 16))
    
    
    p2 <- ggplot() + 
      geom_line( aes(x = time, y = cuminc, color = treatment), data = cumincplotdf, linewidth = 2)+
      xlab("Time in weeks") + 
      ylab("Cumulative Incidence (%)") + 
      scale_x_continuous(breaks = seq(0, 14, 1)*8, minor_breaks = NULL, limits = c(0, 104))+
      scale_y_continuous(breaks = seq(0, 0.1, 0.01), labels = seq(0, 0.1, 0.01)*100, minor_breaks = NULL, limits = c(0, 0.08))+
      scale_colour_manual(breaks = c("VRC01", "Placebo"), values = c("#F8766D", "#619CFF")) +
      annotate(geom ="text", x = 30, y = 0.05, label= paste0("Overall PE = ", PEatTau," at Week ", tau), size = 5)+
      theme_bw()+
      ggtitle(titleLabel) +
      theme(legend.title = element_blank(),
            legend.position = c(0.7, 0.9),
            legend.direction = "horizontal",
            legend.background = element_blank(),
            legend.text = element_text (size = 14),
            title = element_text(size = 14),
            strip.text.y =  element_text (size = 16),
            strip.text.x =  element_text (size = 16),
            axis.title.x =  element_text (size = 16),
            axis.title.y =  element_text (size = 16),
            axis.text.x =  element_text (size = 16),
            axis.text.y =  element_text (size = 16))
    
    ggsave(file.path(figDir, paste0("PE_plot_", scenario, "_", model, ".pdf")), p1, width = 6, height = 5)
    ggsave(file.path(figDir, paste0("cuminc_plot_", scenario, "_", model, ".pdf")), p2, width = 6, height = 5)
    
  }
}

