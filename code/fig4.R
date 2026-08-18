library(dplyr)
library(ggplot2)
library(grid)

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")


# input directory and file names
SURV703 = file.path(datDir, 'v703_survival.csv')
SURV704 = file.path(datDir, 'v704_survival.csv')
adataFile2 <- file.path(datDir, "v704_cuminc_wk104_pool.csv")
adataFile3 <- file.path(datDir, "v704_cuminc_wk104_ind.csv")
adataFile4 <- file.path(datDir, "v703_cuminc_wk104_pool.csv")
adataFile5 <- file.path(datDir, "v703_cuminc_wk104_ind.csv")

# output directory and file names
pdfFileSave <- c(file.path(figDir, "v704_cuminc_wk104_foot.pdf"), 
                 file.path(figDir, "v703_cuminc_wk104_foot.pdf"))

# source input data
cuminc_pool_704 <- read.csv(adataFile2, stringsAsFactors = FALSE) %>%
  filter(rx_pool == 'T1+T2') %>% # Placebo data is also in the individual file
  rename(level=rx_pool) %>%
  mutate(protocol='HVTN 704') %>%
  dplyr::select(protocol, level, time, cuminc, n.risk, n.event)
cuminc_ind_704 <- read.csv(adataFile3, stringsAsFactors = FALSE) %>%
  rename(level=rx_code) %>%
  mutate(protocol='HVTN 704') %>%
  dplyr::select(protocol, level, time, cuminc, n.risk, n.event)
cuminc_pool_703 <- read.csv(adataFile4, stringsAsFactors = FALSE) %>%
  filter(rx_pool == 'T1+T2') %>% # Placebo data is also in the individual file
  rename(level=rx_pool) %>%
  mutate(protocol='HVTN 703') %>%
  dplyr::select(protocol, level, time, cuminc, n.risk, n.event)
cuminc_ind_703 <- read.csv(adataFile5, stringsAsFactors = FALSE) %>%
  rename(level=rx_code) %>%
  mutate(protocol='HVTN 703') %>%
  dplyr::select(protocol, level, time, cuminc, n.risk, n.event)

# combine data
cuminc = rbind(cuminc_pool_704,
               cuminc_ind_704,
               cuminc_pool_703,
               cuminc_ind_703) %>%
  mutate(time_wks = time/7,
         level = factor(level, 
                        levels = c('C3', 'T1+T2','T1','T2')),
         protocol = factor(protocol, levels = c('HVTN 704', 'HVTN 703'))) %>%
  arrange(protocol, level, time) %>% 
  group_by(protocol, level) %>%
  mutate( n.event.cum = cumsum(ifelse(is.na(n.event), 0, n.event)) ) %>%
  ungroup()

         
# create footnotes 

# define font sizes
font_size=18
foot_font_size=4.5
foot_start_line=7

for( i in 1:2 ) {
  
  if( i==1 ) {
    dat = cuminc %>% filter(protocol=='HVTN 704')
    surv = read.csv(SURV704)
  } else {
    dat = cuminc %>% filter(protocol=='HVTN 703')
    surv = read.csv(SURV703)
  }
  surv = surv %>% filter(surv$efficacy_flag==1)
  
  max_wks = max(surv$fudayswk104[which(surv$statuswk104==1)])/7
  wks <- surv$fudayswk104/7
  f_times <- c(seq(0, 96, 16), max_wks)
  
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
  v.low.times <- wks[surv$rx_code=="T1" ]
  v.low.event.times <- wks[surv$statuswk104==1 & surv$rx_code=="T1"]
  nrisk.v.low  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.low.times)
  nevent.v.low <- sapply( f_times, function(T, t) sum( t <= T ), t=v.low.event.times)
  f_vLow_risk <- data.frame(time_wks = f_times,
                            level = "C3",
                            n = nrisk.v.low)
  f_vLow_event <- data.frame(time_wks = f_times,
                             level = "C3",
                             n = nevent.v.low)
  
  v.high.times <- wks[surv$rx_code=="T2" ]
  v.high.event.times <- wks[surv$statuswk104==1 & surv$rx_code=="T2"]
  nrisk.v.high  <- sapply( f_times, function(T, t) sum( t >= T ), t=v.high.times)
  nevent.v.high <- sapply( f_times, function(T, t) sum( t <= T ), t=v.high.event.times)
  f_vHigh_risk <- data.frame(time_wks = f_times,
                             level = "C3",
                             n = nrisk.v.high)
  f_vHigh_event <- data.frame(time_wks = f_times,
                              level = "C3",
                              n = nevent.v.high)
  
  ## lastly for placebo
  p.times <- wks[surv$rx_code == "C3" ]
  p.event.times <- wks[ surv$statuswk104 == 1  &  surv$rx_code == "C3" ]
  nrisk.p  <- sapply( f_times, function(T, t) sum( t >= T ), t=p.times)
  nevent.p <- sapply( f_times, function(T, t) sum( t <= T ), t=p.event.times)
  f_p_risk <- data.frame(time_wks = f_times,
                         level = "C3",
                         n = nrisk.p)
  f_p_event <- data.frame(time_wks = f_times,
                          level = "C3",
                          n = nevent.p)
  
  # plot Pooled VRC01 vs. Placebo w/ footnotes
  p_pool_foot <-
    
    # main plot
    ggplot(data=dat, aes(x=time_wks, y=cuminc, color=level, linetype=level)) +
    scale_x_continuous(name="Weeks since Enrollment",
                       #limits=c(0, 114),
                       breaks=seq(0, 114, by=8),
                       labels=seq(0, 114, by=8)) +
    scale_y_continuous(name="Cumulative Incidence (%)",
                       #limits=c(0, 0.09),
                       breaks=seq(0, 0.09, 0.01),
                       labels=seq(0, 9, 1)) +
    geom_step(lwd=1) +
    scale_color_manual(values = c("C3" = "black", "T1+T2" = "red", "T2" = "blue", "T1" = "forestgreen"), 
                        name="",
                       labels = c("Placebo", "VRC01 pooled", "VRC01, low dose", "VRC01, high dose")) +
    scale_linetype_manual(values = c("C3" = "dotted", "T1+T2" = "solid", "T2" = "solid", "T1" = "dashed"),
                          name="",
                          labels = c("Placebo", "VRC01 pooled", "VRC01, low dose", "VRC01, high dose")) +
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
          plot.margin = unit(c(1,1,6,3), "cm"), #trbl
          legend.position = c(0.25, 0.9),
          legend.text=element_text(size=font_size),
          legend.key.width = unit(2,"cm")) +
    coord_cartesian(xlim=c(0,114), ylim=c(0, 0.09), clip = "off") +
  
    # number at risk
    geom_text(data=f_p_risk, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+1.5) +
    geom_text(data=f_vPool_risk, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+3) +
    geom_text(data=f_vLow_risk, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+4.5) +
    geom_text(data=f_vHigh_risk, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+6) +
    
  
    # number of events
    geom_text(data=f_p_event, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+9.5) +
    geom_text(data=f_vPool_event, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+11) +
    geom_text(data=f_vLow_event, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+12.5) +
    geom_text(data=f_vHigh_event, aes(label=n, x=time_wks, y=-Inf), 
              size=foot_font_size, show.legend = FALSE, vjust=foot_start_line+14) +
    
    annotate("text", x = -34, y = -.02, label = "Number at Risk", 
             hjust = 0, vjust = 1, color = "black", fontface = "bold") +
    annotate("text", x = -30, y = -.024, label = "Placebo", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.028, label = "VRC01 pooled", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.032, label = "VRC01, low dose", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.036, label = "VRC01, high dose", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -34, y = -.04, label = "Cumulative Number of HIV-1 Endpoints", 
             hjust = 0, vjust = 1, color = "black", fontface = "bold") +
    annotate("text", x = -30, y = -.044, label = "Placebo", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.048, label = "VRC01 pooled", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.052, label = "VRC01, low dose", 
             hjust = 0, vjust = 1, color = "black") +
    annotate("text", x = -30, y = -.056, label = "VRC01, high dose", 
             hjust = 0, vjust = 1, color = "black") 
    
  g_pool_foot=ggplot_gtable(ggplot_build(p_pool_foot))
  
  pdf(pdfFileSave[i], width = 8.5, height = 8.5)
  grid::grid.draw(g_pool_foot)
  dev.off()

}

q(save = "no")







