library(dplyr)
library(ggplot2)
library(grid)

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data_final" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")

# input directory and file names
adataFile1 <- file.path(dat2Dir, "amp_survival.csv")
adataFile2 <- file.path(datDir, "amp_cir_wk104_pool.csv")

# output directory and file names
pdfFileSave <- c(file.path(figDir, "amp_cuminc_wk104_pool.pdf"))

# source input data
surv <- read.csv(adataFile1, stringsAsFactors = FALSE) %>%
  filter(efficacy_flag == 1)
cir_pool <- read.csv(adataFile2, stringsAsFactors = FALSE)

# reformat 'CIR' data for plotting, taking steps to ensure plots
# end at tau (our censoring time point) rather than the last event time
tau <- max(surv$fudayswk104)
tau_wks <- tau/7

    ## first for pooled VRC01 vs. control
    cuminc_vPool <-
      cir_pool %>%
      add_row(time = tau) %>%
      mutate(level = if_else(time==tau, lag(cmpLvl), cmpLvl),
             cuminc = if_else(time==tau, lag(cuminc.cmp), cuminc.cmp)) %>%
      select(time, level, cuminc)
    
    cuminc_p <-
      cir_pool %>%
      add_row(time = tau) %>%
      mutate(level = if_else(time==tau, lag(refLvl), refLvl),
             cuminc = if_else(time==tau, lag(cuminc.ref), cuminc.ref)) %>%
      select(time, level, cuminc)

    cuminc_pool <- 
      bind_rows(cuminc_vPool, cuminc_p) %>%
      mutate(level = factor(level, levels=c("T1+T2", "C3")),
             time_wks = time/7)

# create footnotes using 'surv' data
wks <- surv$fudayswk104/7
f_times <- c(seq(0, 96, 16), tau_wks)

    ## first for pooled VRC01
    v.pool.times <- wks[surv$rx_pool=="T1+T2"]
    v.pool.event.times <- wks[surv$statuswk104==1 & surv$rx_pool=="T1+T2"]
    nrisk.v.pool  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.pool.times)
    nevent.v.pool <- sapply( f_times, function(T, t) sum( t <= T ), t=v.pool.event.times)
    f_vPool_risk <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nrisk.v.pool)
    f_vPool_event <- data.frame(time_wks = f_times,
                                level = "C3",
                                n = nevent.v.pool)
      
    ## then for individual VRC01
    v.low.times <- wks[surv$prot_rx %in% c("HVTN 703 T1", "HVTN 704 T1")]
    v.low.event.times <- wks[surv$statuswk104==1 & surv$prot_rx %in% c("HVTN 703 T1", "HVTN 704 T1")]
    nrisk.v.low  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.low.times)
    nevent.v.low <- sapply( f_times, function(T, t) sum( t <= T ), t=v.low.event.times)
    f_vLow_risk <- data.frame(time_wks = f_times,
                              level = "C3",
                              n = nrisk.v.low)
    f_vLow_event <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nevent.v.low)
    
    v.high.times <- wks[surv$prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2")]
    v.high.event.times <- wks[surv$statuswk104==1 & surv$prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2")]
    nrisk.v.high  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.high.times)
    nevent.v.high <- sapply( f_times, function(T, t) sum( t <= T ), t=v.high.event.times)
    f_vHigh_risk <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nrisk.v.high)
    f_vHigh_event <- data.frame(time_wks = f_times,
                                level = "C3",
                                n = nevent.v.high)
    
    ## lastly for placebo
    p.times <- wks[surv$prot_rx %in% c("HVTN 703 C3", "HVTN 704 C3")]
    p.event.times <- wks[ surv$statuswk104 == 1  & surv$prot_rx %in% c("HVTN 703 C3", "HVTN 704 C3")]
    nrisk.p  <- sapply( f_times, function(T, t) sum( t >= T ), t=p.times)
    nevent.p <- sapply( f_times, function(T, t) sum( t <= T ), t=p.event.times)
    f_p_risk <- data.frame(time_wks = f_times,
                              level = "C3",
                              n = nrisk.p)
    f_p_event <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nevent.p)



# update font sizes for plots w/o footnotes
font_size=32


# plot Pooled VRC01 vs. Control w/o footnotes
p_pool <-
  
  ggplot(data=cuminc_pool, aes(x=time_wks, y=cuminc, color=level, linetype=level)) +
  scale_x_continuous(name="Weeks since Enrollment",
                     limits=c(0, 116),
                     breaks = seq(0, 112, by=16),
                     labels = seq(0, 112, by=16)) +
  scale_y_continuous(name="Cumulative Incidence (%)",
                     limits=c(0, 0.07),
                     breaks=seq(0, 0.07, 0.01),
                     labels=seq(0, 7, 1)) +
  geom_step(lwd=2) +
  scale_color_manual(values = c("T1+T2" = "red", "C3" = "black"), 
                     name="",
                     labels = c("T1+T2" = "Pooled VRC01", "C3" = "Control")) +
  scale_linetype_manual(values = c("T1+T2" = "solid", "C3" = "dotted"),
                        name="",
                        labels = c("T1+T2" = "Pooled VRC01", "C3" = "Control")) +
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
        legend.position = c(0.35, 0.85),
        legend.text=element_text(size=font_size),
        legend.key.width = unit(2,"cm"))

ggsave(pdfFileSave[1], plot = p_pool, width=8.5, height=8.5)


q(save = "no")
