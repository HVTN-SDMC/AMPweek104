library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(lubridate)
library(survival)

library(here)
here::i_am("README.md")
repoDir <- here::here()
macroDir <- file.path(repoDir, "code/macro")
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")

# input directory and file names
# adataFile1 <- file.path(dat2Dir, "amp_survival.csv")
adataFile2 <- file.path(dat2Dir, "v704_survival_wk104_neut.csv")
adataFile3 <- file.path(dat2Dir, "v703_survival_wk104_neut.csv")
nabFile704 <- file.path(dat2Dir, "VTN704_breakthrough_NAb_20200723.txt")
nabFile703 <- file.path(dat2Dir, "VTN703_breakthrough_NAb_20200722.txt")


# output directory and file names
pdfFileSave <- file.path(figDir, "amp_ic80ls_over_time_wk104.pdf")

# source input data
# surv <- read.csv(adataFile1, stringsAsFactors = FALSE)
survneut4 <- read.csv(adataFile2, stringsAsFactors = FALSE)
survneut3 <- read.csv(adataFile3, stringsAsFactors = FALSE)
nab4 <- read.csv(nabFile704, sep="\t")
nab3 <- read.csv(nabFile703, sep="\t")

# calculate time from enrollment to IC80 sample draw date
# and merge with IC80 value based on LS variant, primary endpoints only

nab <-
  rbind(nab3, nab4) %>%
  rename(pub_id=isolate_pubid, drawdy=isolate_drawdt) %>%
  arrange(pub_id, drawdy) %>%
  distinct(pub_id, .keep_all=TRUE) %>%
  select(pub_id, drawdy)

dat <-
  bind_rows(survneut3, survneut4) %>%
  filter(hiv1event==1 & nisolates>0) %>%
  mutate(gmt80lsn=as.numeric(gsub(">", "", gmt80ls))) %>%
  select(tx, pub_id, gmt80ls, gmt80lsn) %>%
  left_join(nab, by="pub_id") %>%
  mutate(#fudays = as.numeric(drawdt-enrdt),
         fudays = drawdy,
         fuwks = fudays/7,
         tx=factor(tx, levels=c("T2", "T1", "C3")))

dat2 <-
  dat %>%
  mutate(logtiterL = ifelse(grepl("^<", gmt80ls), NA, log(gmt80lsn)),
         logtiterR = ifelse(grepl("^>", gmt80ls), NA, log(gmt80lsn)),
         timeperiod = case_when(fuwks <= 16 ~ 16,
                                fuwks <= 32 ~ 32,
                                fuwks <= 48 ~ 48,
                                fuwks <= 64 ~ 64,
                                fuwks <= 80 ~ 80,
                                fuwks <= 96 ~ 96,
                                TRUE ~ 114))
  
# create figure footnote with GMT (95%) at 16 week intervals by arm

    # function to calculate GMT (95% CI)
    calc_gmt = function(dat, alpha=0.05){
      fit = survreg(Surv(dat$logtiterL, dat$logtiterR, type='interval2') ~ 1, dist = "gaussian")
      gmt = exp(coef(fit))
      up.gmt = exp(confint(fit, level=1-alpha)[2])
      lo.gmt = exp(confint(fit, level=1-alpha)[1])
      return(c(gmt, lo.gmt, up.gmt))
    }

    # loop through timepoints and treatment groups
    txgroups = c("T2", "T1", "C3")
    timeperiods = c(16, 32, 48, 64, 80, 96, 114)
    
    foot = as.data.frame(matrix(data=NA, nrow=length(txgroups)*length(timeperiods), ncol=5))
    colnames(foot) <- c("tx", "timeperiod", "gmt", "lo.gmt", "up.gmt")
    counter = 1
    
    for (i in txgroups){
      for (j in timeperiods){
        datsub=subset(dat2, tx == i & timeperiod == j)
        foot[counter, 1] <- i
        foot[counter, 2] <- j
        foot[counter, 3:5] <- calc_gmt(datsub)
        counter=counter+1
      }
    }
    
    foot <- 
      foot %>% 
      mutate(tx=factor(tx, levels=c("T2", "T1", "C3")),
             timeperiodX=timeperiod-8,
             gmtlab=format(round(gmt), nsmall=0),
             up.gmtlab=format(round(up.gmt, 0), nsmall=0),
             lo.gmtlab=format(round(lo.gmt, 0), nsmall=0))

# build plot with raw values

font_size=15

p1 <-
  ggplot(data=dat, aes(x=fuwks, y=gmt80lsn, color=tx)) +
  geom_point(size=3, position = position_jitter(w = 0, h = 0)) +
  scale_x_continuous(name="Weeks from Enrollment to IC80 Sample Collection Date",
                     limits=c(0, 120),
                     breaks = seq(0, 120, by=16),
                     labels = seq(0, 120, by=16)) +
  scale_y_continuous(name=expression("IC80 ("*mu*"g/mL)"),
                     limits=c(min(dat$gmt80lsn)-0.05, 100),
                     breaks=c(0.1, 0.3, 1, 3, 10, 30, 100),
                     labels=c("0.1", "0.3", "1", "3", "10", "30", expression("">="100")),
                     trans="log10") +
  geom_smooth(method="loess", se=FALSE) +
  scale_color_manual(values = c("T1" = "#1749FF",
                                "T2" = "#D92321",
                                "C3" = "#0AB7C9"),
                     name="",
                     labels = c("T1" = "VRC01 10 mg/kg",
                                "T2" = "VRC01 30 mg/kg",
                                "C3" = "Control")) +
  labs(title="Pooled AMP Trials - MITT Endpoints\nLeast Sensitive Variant") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=10)),
        axis.title.y = element_text(size=font_size, margin=margin(r=0)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,0,1), "cm"), #trbl
        legend.position = c(0.5, 0.04),
        legend.text=element_text(size=font_size*.7),
        legend.key.width = unit(2,"cm"),
        legend.background = element_rect(fill="transparent"),
        legend.direction = "horizontal") 

p2 <-
  ggplot(data=foot, aes(x=timeperiodX, y=gmt, ymin=lo.gmt, ymax=up.gmt, color=tx)) +
  geom_linerange(position=position_dodge(width=12)) + 
  geom_pointrange(position=position_dodge(width=12)) +
  geom_text(aes(y=gmt, label=gmtlab), position=position_dodge(width=12), hjust=1.35, vjust=-0.5, size=3, show.legend = FALSE) +
  geom_text(aes(y=lo.gmt, label=lo.gmtlab), position=position_dodge(width=12), hjust=1.25, size=3, show.legend = FALSE) +
  geom_text(aes(y=up.gmt, label=up.gmtlab), position=position_dodge(width=12), hjust=1.25, vjust=-0.5, size=3, show.legend = FALSE) +
  scale_x_continuous(name="Weeks from Enrollment to IC80 Sample Collection Date",
                     limits=c(0, 120),
                     breaks = seq(16, 112, by=16) - 8,
                     labels = c("(0,16]", "(16,32]", "(32,48]", "(48,64]", "(64,80]", "(80,96]", "(96,114]")) +
  scale_y_continuous(name=expression("Geometric Mean IC80 Titer, 95% CI ("*mu*"g/mL)"),
                     limits=c(0.1, 1300),
                     breaks=c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000),
                     labels=c("0.1", "0.3", "1", "3", "10", "30", "100", "300", "1000"),
                     trans="log10") +
  scale_color_manual(values = c("T1" = "#1749FF",
                                "T2" = "#D92321",
                                "C3" = "#0AB7C9"),
                     name="",
                     labels = c("T1" = "VRC01 10 mg/kg",
                                "T2" = "VRC01 30 mg/kg",
                                "C3" = "Control")) +
  # labs(title="Pooled AMP Trials - MITT Endpoints") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=10)),
        axis.title.y = element_text(size=font_size, margin=margin(r=0)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(0,1,1,1), "cm"), #trbl
        legend.position = c(0.5, 0.05),
        legend.text=element_text(size=font_size*0.8),
        legend.key.width = unit(2,"cm"),
        legend.background = element_rect(fill="transparent"),
        legend.direction = "horizontal") 
  
pdf(pdfFileSave, width=8.5, height=11)
grid.arrange(p1, p2, ncol = 1, nrow = 2)
dev.off()

q(save = "no")
