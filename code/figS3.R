library(dplyr)
library(ggplot2)
library(grid)

# load macros
source("macro/confbandsurv_2020Jan.R") 

# input directory and file names
adataDir <- "/Volumes/trials/vaccine/p704/analysis/efficacy/code/masking//adata"
adataFile1 <- "amp_survival_postwk80.csv"
adataFile2 <- "amp_cir_pool_postwk80_trunc.csv"
adataFile3 <- "amp_cir_ind_postwk80_trunc.csv"

# output directory and file names
csvDir <- "../output/tables/"
csvFile <- c("amp_confbandsurv_efficacy_pool_postwk80_trunc.csv",
             "amp_confbandsurv_efficacy_high_postwk80_trunc.csv",
             "amp_confbandsurv_efficacy_low_postwk80_trunc.csv")
csvFileSave <- file.path(csvDir, csvFile)

pdfDir <- "../output/figures/"
pdfFile <- c("amp_efficacy_pool_postwk80_trunc_supp.pdf", 
             "amp_efficacy_high_postwk80_trunc_supp.pdf",
             "amp_efficacy_low_postwk80_trunc_supp.pdf")
pdfFileSave <- file.path(pdfDir, pdfFile)

# source input data and divide 'surv' into three datasets to obtain simultaneous
# bounds for each of the comparisons of interest 
surv <- read.csv(file.path(adataDir, adataFile1), stringsAsFactors = FALSE) %>%
  mutate(fudays_postwk80 = pmin(24*7, fudays_postwk80)) # use truncated survival 
cir_pool <- read.csv(file.path(adataDir, adataFile2), stringsAsFactors = FALSE)
cir_ind <- read.csv(file.path(adataDir, adataFile3), stringsAsFactors = FALSE)

    surv_pool <- 
      surv %>%
      mutate(timeD = fudays_postwk80,
             delta = status_postwk80, 
             Rx = ifelse(rx_pool=="T1+T2", 2, 1)) %>%
      select(timeD, delta, Rx)
    
    surv_high <- 
      surv %>%
      filter(prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2", "HVTN 703 C3", "HVTN 704 C3")) %>%
      mutate(timeD = fudays_postwk80,
             delta = status_postwk80,  
             Rx = ifelse(prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2"), 2, 1)) %>%
      select(timeD, delta, Rx)
    
    surv_low <- 
      surv %>%
      filter(prot_rx %in% c("HVTN 703 T1", "HVTN 704 T1", "HVTN 703 C3", "HVTN 704 C3")) %>%
      mutate(timeD = fudays_postwk80,
             delta = status_postwk80, 
             Rx = ifelse(prot_rx %in% c("HVTN 703 T1", "HVTN 704 T1"), 2, 1)) %>%
      select(timeD, delta, Rx)

# set max time through which to estimate simultaneous bounds based on 150 at risk rule
# this is where the survival truncation occurs
tau <- 24*7

# obtain survival results for each of three analysis datasets
run_Confbandsurv <- function(dat){
  Confbandsurv(dat$timeD, dat$delta, dat$Rx, N = 500, t1=0, t2=tau, ngrid=500, contrast="VE")
}

result_surv_pool <- data.frame(run_Confbandsurv(surv_pool))
result_surv_high <- data.frame(run_Confbandsurv(surv_high))
result_surv_low <- data.frame(run_Confbandsurv(surv_low))

# convert time to weeks
daysPerWeek <- 7
result_surv_pool$timesunique <- result_surv_pool$timesunique/daysPerWeek
result_surv_high$timesunique <- result_surv_high$timesunique/daysPerWeek
result_surv_low$timesunique <- result_surv_low$timesunique/daysPerWeek

# export confbandsurv data
write.csv(result_surv_pool, csvFileSave[1], row.names=FALSE)
write.csv(result_surv_high, csvFileSave[2], row.names=FALSE)
write.csv(result_surv_low,  csvFileSave[3], row.names=FALSE)

# combine 'CIR' data with 95% simultaneous bounds from confbandsurv data for plotting plotting, 
# taking steps to ensure plots end at tau (our censoring time point) rather than the last event time
tau_wks <- tau/7

    ## first for pooled VRC01 vs. control
    eff_pool <-
      cir_pool %>%
      add_row(time = tau) %>%
      mutate(time_wks = time/7,
             eff = if_else(time==tau, lag(eff), eff),
             lo.eff = if_else(time==tau, lag(lo.eff), lo.eff),
             up.eff = if_else(time==tau, lag(up.eff), up.eff)) %>%
      select(time_wks, eff, lo.eff, up.eff)

    ## then for high VRC01 vs. control
    eff_high <-
      cir_ind %>%
      filter(comparison=="T2 vs. C3") %>%
      add_row(time = tau) %>%
      mutate(time_wks = time/7,
             eff = if_else(time==tau, lag(eff), eff),
             lo.eff = if_else(time==tau, lag(lo.eff), lo.eff),
             up.eff = if_else(time==tau, lag(up.eff), up.eff)) %>%
      select(time_wks, eff, lo.eff, up.eff)
    
    ## lastly for low VRC01 vs. control
    eff_low <-
      cir_ind %>%
      filter(comparison=="T1 vs. C3") %>%
      add_row(time = tau) %>%
      mutate(time_wks = time/7,
             eff = if_else(time==tau, lag(eff), eff),
             lo.eff = if_else(time==tau, lag(lo.eff), lo.eff),
             up.eff = if_else(time==tau, lag(up.eff), up.eff)) %>%
      select(time_wks, eff, lo.eff, up.eff)
    


# update font sizes for plots w/o footnotes
font_size=24

# plot Pooled VRC01 vs. Control w/o footnotes
p_pool <-
  
  ggplot(data=eff_pool, aes(x=time_wks, y=eff)) +
  scale_x_continuous(name="Weeks since Week 80",
                     limits=c(0, 24),
                     breaks=c(0,4,8,12,16,20,24),
                     labels=c(0,4,8,12,16,20,24)) +
  scale_y_continuous(name="Prevention Efficacy (%)",
                     limits=c(-2, 1),
                     breaks=seq(-2, 1, 0.25),
                     labels=seq(-200, 100, 25)) +
  geom_step(linetype="solid", lwd=1) +
  geom_step(aes(x=time_wks, y=lo.eff), linetype="dashed", lwd=1) +
  geom_step(aes(x=time_wks, y=up.eff), linetype="dashed", lwd=1) +
  # geom_line(data=result_surv_pool, aes(x=timesunique, y=lowband95), linetype="dotted", lwd=1) +
  # geom_line(data=result_surv_pool, aes(x=timesunique, y=upband95), linetype="dotted", lwd=1) +
  geom_step(y=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=20)),
        axis.title.y = element_text(size=font_size, margin=margin(r=20)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), #trbl
        legend.position = 'inside',
        legend.position.inside = c(0.2,0.9),
        legend.text=element_text(size=font_size)) +
  labs(title="Dose Pooled vs. Placebo") 

ggsave(pdfFileSave[1], plot = p_pool, width=11, height=8.5)

# plot High VRC01 vs. Control w/o footnotes
p_high <-
  
  ggplot(data=eff_high, aes(x=time_wks, y=eff)) +
  scale_x_continuous(name="Weeks since Week 80",
                     limits=c(0, 24),
                     breaks=c(0,4,8,12,16,20,24),
                     labels=c(0,4,8,12,16,20,24)) +
  scale_y_continuous(name="Prevention Efficacy (%)",
                     limits=c(-2, 1),
                     breaks=seq(-2, 1, 0.25),
                     labels=seq(-200, 100, 25)) +
  geom_step(linetype="solid", lwd=1) +
  geom_step(aes(x=time_wks, y=lo.eff), linetype="dashed", lwd=1) +
  geom_step(aes(x=time_wks, y=up.eff), linetype="dashed", lwd=1) +
  # geom_line(data=result_surv_high, aes(x=timesunique, y=lowband95), linetype="dotted", lwd=1) +
  # geom_line(data=result_surv_high, aes(x=timesunique, y=upband95), linetype="dotted", lwd=1) +
  geom_step(y=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=20)),
        axis.title.y = element_text(size=font_size, margin=margin(r=20)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm")#, #trbl
        #legend.position = 'inside',
        #legend.position.inside = c(0.2, 0.9),
        #legend.text=element_text(size=font_size)
        ) +
  labs(title="High Dose vs. Placebo") 

ggsave(pdfFileSave[2], plot = p_high, width=11, height=8.5)

# plot Low VRC01 vs. Control w/o footnotes
p_low <-
  
  ggplot(data=eff_low, aes(x=time_wks, y=eff)) +
  scale_x_continuous(name="Weeks since Week 80",
                     limits=c(0, 24),
                     breaks=c(0,4,8,12,16,20,24),
                     labels=c(0,4,8,12,16,20,24)) +
  scale_y_continuous(name="Prevention Efficacy (%)",
                     limits=c(-2, 1),
                     breaks=seq(-2, 1, 0.25),
                     labels=seq(-200, 100, 25)) +
  geom_step(linetype="solid", lwd=1) +
  geom_step(aes(x=time_wks, y=lo.eff), linetype="dashed", lwd=1) +
  geom_step(aes(x=time_wks, y=up.eff), linetype="dashed", lwd=1) +
  # geom_line(data=result_surv_low, aes(x=timesunique, y=lowband95), linetype="dotted", lwd=1) +
  # geom_line(data=result_surv_low, aes(x=timesunique, y=upband95), linetype="dotted", lwd=1) +
  geom_step(y=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=20)),
        axis.title.y = element_text(size=font_size, margin=margin(r=20)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), #trbl
        legend.position = 'inside',
        legend.position.inside = c(0.2, 0.9),
        legend.text=element_text(size=font_size)) +
  labs(title="Low Dose vs. Placebo") 

ggsave(pdfFileSave[3], plot = p_low, width=11, height=8.5)

q(save = "no")
