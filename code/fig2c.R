library(dplyr)
library(ggplot2)
library(grid)

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")

# input directory and file names
adataFile <- file.path(datDir, "amp_cir_wk104_pool_trunc.csv")

# output directory and file names
pdfFileSave <- file.path(figDir, "amp_efficacy_wk104_pool_trunc.pdf")

cir_pool <- 
  read.csv(adataFile, stringsAsFactors = FALSE) 

# reformat 'CIR' data for plotting, taking steps to ensure plots
# end at tau (our censoring time point) rather than the last event time
tau <- max(cir_pool$time)

## then for pooled VRC01 vs. control - IC80
eff_Pool <-
  cir_pool %>%
  mutate(time_wks = time/7) %>%
  select(time_wks, eff, lo.eff, up.eff)
    
# update font sizes for plots w/o footnotes
font_size=24

eff.last = tail(subset(eff_Pool), 1)
pts.last = data.frame( x = c(eff.last$time_wks, eff.last$time_wks),
                       y = c(eff.last$lo.eff, eff.last$up.eff) )

p_pool <-
  
  # main plot
  ggplot(data=eff_Pool, aes(x=time_wks, y=eff)) +
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


ggsave(pdfFileSave, plot = grid.draw(p_pool), width=8.5, height=8.5)


q(save = "no")
