library(dplyr)
library(ggplot2)
library(grid)

# input directory and file names
adataDir <- "/Volumes/trials/vaccine/p704/analysis/efficacy/code/masking/adata"
adataFile1 <- "amp_survival_postwk80.csv"
adataFile2 <- "amp_cir_pool_postwk80_v2.csv"
adataFile3 <- "amp_cir_ind_postwk80_v2.csv"

# output directory and file names
pdfDir <- "../output/figures"
pdfFile <- "amp_cuminc_postwk80_trunc.pdf"
pdfFileSave <- file.path(pdfDir, pdfFile)

# source input data
surv <- read.csv(file.path(adataDir, adataFile1), stringsAsFactors = FALSE)
cir_pool <- read.csv(file.path(adataDir, adataFile2), stringsAsFactors = FALSE)
cir_ind <- read.csv(file.path(adataDir, adataFile3), stringsAsFactors = FALSE)


# note the infections that are past tau
# reformat 'CIR' data for plotting, taking steps to ensure plots
# end at tau (our censoring time point) rather than the last event time
surv_infected = subset(surv, status_postwk80==1)
tau <- max(surv_infected$fudays_postwk80)
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

    ## then for individual VRC01 vs. control
    cuminc_vInd_low <-
      cir_ind %>%
      filter(cmpLvl == 'T1') %>%
      add_row(time = tau) %>%
      mutate(level = if_else(time==tau, lag(cmpLvl), cmpLvl),
             cuminc = if_else(time==tau, lag(cuminc.cmp), cuminc.cmp)) %>%
      select(time, level, cuminc)
    
    cuminc_vInd_high <-
      cir_ind %>%
      filter(cmpLvl == 'T2') %>%
      add_row(time = tau) %>%
      mutate(level = if_else(time==tau, lag(cmpLvl), cmpLvl),
             cuminc = if_else(time==tau, lag(cuminc.cmp), cuminc.cmp)) %>%
      select(time, level, cuminc)
    
    cuminc_ind <- 
      bind_rows(cuminc_vInd_high, cuminc_vInd_low, cuminc_p) %>%
      mutate(level = factor(level, levels=c("T2", "T1", "C3")),
             time_wks = time/7)
    
cuminc = rbind(cuminc_pool %>% filter(level=='T1+T2'), # otherwise we have C3 twice
               cuminc_ind) %>%
  mutate(level = factor(level, levels = c('C3', 'T1+T2', 'T1', 'T2')))


# create footnotes using 'surv' data
wks <- surv$fudays_postwk80/7
f_times <- c(seq(0, 20, 4), tau_wks)

    ## first for pooled VRC01
    v.pool.times <- wks[surv$rx_pool=="T1+T2"]
    v.pool.event.times <- wks[surv$status_postwk80==1 & surv$rx_pool=="T1+T2"]
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
    v.low.event.times <- wks[surv$status_postwk80==1 & surv$prot_rx %in% c("HVTN 703 T1", "HVTN 704 T1")]
    nrisk.v.low  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.low.times)
    nevent.v.low <- sapply( f_times, function(T, t) sum( t <= T ), t=v.low.event.times)
    f_vLow_risk <- data.frame(time_wks = f_times,
                              level = "C3",
                              n = nrisk.v.low)
    f_vLow_event <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nevent.v.low)
    
    v.high.times <- wks[surv$prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2")]
    v.high.event.times <- wks[surv$status_postwk80==1 & surv$prot_rx %in% c("HVTN 703 T2", "HVTN 704 T2")]
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
    p.event.times <- wks[ surv$status_postwk80 == 1  & surv$prot_rx %in% c("HVTN 703 C3", "HVTN 704 C3")]
    nrisk.p  <- sapply( f_times, function(T, t) sum( t >= T ), t=p.times)
    nevent.p <- sapply( f_times, function(T, t) sum( t <= T ), t=p.event.times)
    f_p_risk <- data.frame(time_wks = f_times,
                              level = "C3",
                              n = nrisk.p)
    f_p_event <- data.frame(time_wks = f_times,
                               level = "C3",
                               n = nevent.p)

# create footnote labels and define font sizes
flabel_risk <- data.frame(time_wks=0, level="C3", n="Number at Risk")
flabel_event <- data.frame(time_wks=0, level="C3", n="Cumulative Number of HIV-1 Endpoints")
flabel_vPool <- data.frame(time_wks=0, level="C3", n="VRC01 pooled")
flabel_vHigh <- data.frame(time_wks=0, level="C3", n="VRC01 high dose")
flabel_vLow <- data.frame(time_wks=0, level="C3", n="VRC01, low dose")
flabel_p <- data.frame(time_wks=0, level="C3", n="Placebo")
font_size=18
foot_font_size=4.5
foot_start_line=7

# plot 
p_foot <-
  
  # main plot
  ggplot(data=cuminc, aes(x=time_wks, y=cuminc, color=level, linetype=level)) +
  scale_x_continuous(name="Weeks since Week 80 Visit",
                     breaks=seq(0, 24, by=4),
                     labels=seq(0, 24, by=4)) +
  scale_y_continuous(name="Cumulative Incidence (%)",
                     breaks=seq(0, 0.04, 0.01),
                     labels=seq(0, 4, 1)) +
  geom_step(lwd=1) +
  scale_color_manual(values = c("T1+T2" = "red", "T2" = "blue", "T1" = "forestgreen", "C3" = "black"), 
                     name="",
                     labels = c("T1+T2" = "VRC01 pooled", "T2" = "VRC01, high dose", "T1"= "VRC01, low dose", "C3" = "Placebo")) +
    scale_linetype_manual(values = c("T1+T2" = "solid", "T2" = "solid", "T1" = "dashed", "C3" = "dotted"),
                        name="",
                        labels = c("T1+T2" = "VRC01 pooled", "T2" = "VRC01, high dose", "T1"= "VRC01, low dose", "C3" = "Placebo")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=10)),
        axis.title.y = element_text(size=font_size, margin=margin(r=10)),
        axis.text.x = element_text(size=font_size, margin=margin(t=5)),
        axis.text.y = element_text(size=font_size, margin=margin(r=5)),
        axis.ticks.length = unit(0.15, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,6.5,3), "cm"), #trbl
        legend.position = 'inside',
        legend.position.inside = c(0.25, 0.9),
        legend.text=element_text(size=font_size),
        legend.key.width = unit(2,"cm")) +
    coord_cartesian(xlim=c(0,25), ylim=c(0,0.04), clip='off') +
    
    annotate("text", x = -5.8, y = -.009, label = "Number at Risk", 
             hjust = 0, vjust = 1, color = "black", fontface = "bold") +
    
    annotate("text", x = -5, y = -.011, label = "Placebo", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_p_risk, aes(label=n, x=time_wks, y=-0.011), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5, y = -.013, label = "VRC01 pooled", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vPool_risk, aes(label=n, x=time_wks, y=-0.013), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5, y = -.015, label = "VRC01, low dose", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vLow_risk, aes(label=n, x=time_wks, y=-0.015), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +

    annotate("text", x = -5, y = -.017, label = "VRC01, high dose", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vHigh_risk, aes(label=n, x=time_wks, y=-0.017), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5.8, y = -.02, label = "Cumulative Number of HIV-1 Endpoints", 
             hjust = 0, vjust = 1, color = "black", fontface = "bold") +
    
    annotate("text", x = -5, y = -.022, label = "Placebo", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_p_event, aes(label=n, x=time_wks, y=-0.022), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5, y = -.024, label = "VRC01 pooled", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vPool_event, aes(label=n, x=time_wks, y=-0.024), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5, y = -.026, label = "VRC01, low dose", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vLow_event, aes(label=n, x=time_wks, y=-0.026), 
              size=foot_font_size, show.legend = FALSE, vjust=1) +
    
    annotate("text", x = -5, y = -.028, label = "VRC01, high dose", 
             hjust = 0, vjust = 1, color = "black") +
    geom_text(data=f_vHigh_event, aes(label=n, x=time_wks, y=-0.028), 
              size=foot_font_size, show.legend = FALSE, vjust=1) 
    


pdf(pdfFileSave, width=11, height=8.5)
grid.draw(p_foot)
dev.off()

q(save = "no")
