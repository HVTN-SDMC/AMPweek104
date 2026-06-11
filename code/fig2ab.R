library(dplyr)
library(ggplot2)
library(grid)

# input directory and file names
mainDataDir = "/Volumes/trials/vaccine/p704/analysis/efficacy/adata"
adataDir <- "/Volumes/trials/vaccine/p704/analysis/efficacy/code/masking/adata"
adataFile1 <- "amp_survival_wk104_tau_neut_gpdx.csv"
adataFile3 <- "amp_cir_wk104_pool_cmpriskIC80ls_2cat_trunc.csv"

# output directory and file names
pdfDir <- "../output/figures/"
pdfFile.sens <- "amp_efficacy_wk104_neut_cmprskIC80_poolls_trunc_sens.pdf"
pdfFile.res <- "amp_efficacy_wk104_neut_cmprskIC80_poolls_trunc_res.pdf"

pdfFileSave <- file.path(pdfDir, c(pdfFile.sens, pdfFile.res))

# source input data
surv <- read.csv(file.path(mainDataDir, adataFile1), stringsAsFactors = FALSE)

cir_pool80 <- 
  read.csv(file.path(adataDir, adataFile3), stringsAsFactors = FALSE) %>%
  mutate(eventType = factor(eventType, levels=c("Sensitive", "Resistant")))


# reformat 'CIR' data for plotting, taking steps to ensure plots
# end at tau (our censoring time point) rather than the last event time
tau <- max(surv$hiv1survday)

## then for pooled VRC01 vs. control - IC80
eff_Pool80 <-
  cir_pool80 %>%
  group_by(eventType) %>%
  mutate(time = if_else(time==max(time), tau, time),
         time_wks = time/7) %>%
  select(eventType, time_wks, eff, lo.eff, up.eff)
    


# update font sizes for plots w/o footnotes
font_size=24

# plot Pooled VRC01 vs. Control w/o footnotes - IC80 
# Sensitive
eff.last = tail(subset(eff_Pool80, eventType=='Sensitive'), 1)
pts.last = data.frame( x = c(eff.last$time_wks, eff.last$time_wks),
                       y = c(eff.last$lo.eff, eff.last$up.eff) )

p_pool80_sens <-
  
  # main plot
  ggplot(data=subset(eff_Pool80, eventType=='Sensitive'), aes(x=time_wks, y=eff)) +
  scale_x_continuous(name="Weeks since Enrollment",
                     limits=c(0, 117),
                     breaks = c(seq(0, 96, by=16), 108),
                     labels = c(seq(0, 96, by=16), 108)) +
  scale_y_continuous(name="Prevention Efficacy (%)",
                     breaks=seq(-1, 1, 0.25),
                     labels=seq(-100, 100, 25)) +
  coord_cartesian(ylim=c(-1,1)) + 
  geom_line(data = pts.last, aes(x = x, y = y), col='gray', linewidth = 2.5) +
  geom_point(data = eff.last, aes(x = time_wks, y = eff), size = 3) +
  geom_point(data = eff.last, aes(x = time_wks, y = lo.eff), size = 3) +
  geom_point(data = eff.last, aes(x = time_wks, y = up.eff), size = 3) +
  geom_text(data = eff.last, aes(x = time_wks, y = eff, label = sprintf("%.1f%%", 100 * eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_text(data = eff.last, aes(x = time_wks, y = lo.eff, label = sprintf("%.1f%%", 100 * lo.eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_text(data = eff.last, aes(x = time_wks, y = up.eff, label = sprintf("%.1f%%", 100 * up.eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_step(linetype="solid", lwd=1) +
  geom_step(aes(y=lo.eff), linetype="dashed", lwd=1) +
  geom_step(aes(y=up.eff), linetype="dashed", lwd=1) +
  geom_step(y=0, color="black") +
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
        legend.position = "none",
        legend.text=element_text(size=font_size)) 

# plot Pooled VRC01 vs. Control w/o footnotes - IC80
# Intermediate-Resistant
eff.last = tail(subset(eff_Pool80, eventType=='Resistant'), 1)
pts.last = data.frame( x = c(eff.last$time_wks, eff.last$time_wks),
                       y = c(eff.last$lo.eff, eff.last$up.eff) )
p_pool80_res <-
  
  # main plot
  ggplot(data=subset(eff_Pool80, eventType=='Resistant'), aes(x=time_wks, y=eff)) +
  scale_x_continuous(name="Weeks since Enrollment",
                     limits=c(0, 117),
                     breaks = c(seq(0, 96, by=16), 108),
                     labels = c(seq(0, 96, by=16), 108)) +
  scale_y_continuous(name="Prevention Efficacy (%)",
                     breaks=seq(-1, 1, 0.25),
                     labels=seq(-100, 100, 25)) +
  coord_cartesian(ylim=c(-1,1)) + 
  geom_line(data = pts.last, aes(x = x, y = y), col='gray', linewidth = 2.5) +
  geom_point(data = eff.last, aes(x = time_wks, y = eff), size = 3) +
  geom_point(data = eff.last, aes(x = time_wks, y = lo.eff), size = 3) +
  geom_point(data = eff.last, aes(x = time_wks, y = up.eff), size = 3) +
  geom_text(data = eff.last, aes(x = time_wks, y = eff, label = sprintf("%.1f%%", 100 * eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_text(data = eff.last, aes(x = time_wks, y = lo.eff, label = sprintf("%.1f%%", 100 * lo.eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_text(data = eff.last, aes(x = time_wks, y = up.eff, label = sprintf("%.1f%%", 100 * up.eff)),
            hjust = -0.3, vjust = 0.5, size = 4) +
  geom_step(linetype="solid", lwd=1) +
  geom_step(aes(y=lo.eff), linetype="dashed", lwd=1) +
  geom_step(aes(y=up.eff), linetype="dashed", lwd=1) +
  geom_step(y=0, color="black") +
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
        legend.position = "none",
        legend.text=element_text(size=font_size)) 


ggsave(pdfFileSave[1], plot = grid.draw(p_pool80_sens), width=8.5, height=8.5)
ggsave(pdfFileSave[2], plot = grid.draw(p_pool80_res), width=8.5, height=8.5)


q(save = "no")
