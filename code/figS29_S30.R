library(tidyverse)
library(lubridate)

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")


# dataFiles
dataFile703 = file.path(dat2Dir, 'v703_survival_wk104_neut.csv')
dataFile704 = file.path(dat2Dir, 'v704_survival_wk104_neut.csv')
dataVL703 = file.path(dat2Dir, 'v704_viral_loads.csv')
dataVL704 = file.path(dat2Dir, 'v703_viral_loads.csv')
dataSurv = file.path(dat2Dir, "amp_survival_postwk80.csv")
dataSieve = file.path(dat2Dir, 'd_wk80_wk104_survival_dataset_sieve.csv')

pdfFileSave = c(file.path(figDir, 'viral_load_plot_postwk80.pdf'),
                file.path(figDir, 'viral_load_vs_IC80_postwk80.pdf'))

dat.703 = read.csv(dataFile703)
dat.704 = read.csv(dataFile704)
dat.704$southAmerica = NULL
dat = rbind(dat.703, dat.704)
dat = subset(dat, select=c('protocol','pub_id','tx','hiv1survday','hiv1event','gmt80ls','gmt80ms'))
dat$tx = factor(dat$tx, levels = c('C3', 'T1', 'T2'), labels=c('Control', '10 mg/kg', '30 mg/kg'))
dat$gmt80ls = as.numeric(sub('>','',dat$gmt80ls))
dat$gmt80ms = as.numeric(sub('>','',dat$gmt80ms))
dat$least_sensitive = ifelse(dat$gmt80ls > 1, 'IC80ls > 1', 'IC80ls <= 1')

# vl data does not include non-MITT either
vl.704 = read.csv(dataVL703)
vl.703 = read.csv(dataVL704)
vl = rbind(vl.703, vl.704)

# get numeric Viral loads converting < values to half the LLOQ (20, 40 or 43)
# and 10,000,000 to 10,000,000
# note, vl is numeric or has < or > with unique values <20, <40, <43, or >10000000
# where <43 is for pub_ID 703-3000 after the first positive result
vl = vl %>% mutate(
  vln = case_when( vl=='<20' ~ '10',
                   vl=='<40' ~ '20',
                   vl=='<43' ~ '21.5',
                   vl=='>10000000' ~ '10000000',
                   TRUE ~ vl))
stopifnot(all(!is.na(as.numeric(vl$vln))))
vl$vln = as.numeric(vl$vln)


# remove negative results
vl = vl %>%
  filter(resultc != 'Not detected')

# only consider viral loads prior to art start
vl = vl %>% 
  filter(is.na(artstartdy) | drawdy < artstartdy)

# make sure data are ordered and take the first positive viral load
vl = vl %>%
  group_by(pub_id) %>%
  slice_head(n=1) %>%
  ungroup()

# restrict to 220 cases
df = dat %>%
  filter(hiv1event==1) %>%
  merge(vl)


pdat =
  df %>%
  filter(!is.na(gmt80ls)) %>%
  mutate(IC80le1 = ifelse(gmt80ls < 1, 1, 0),
         Trt = rx_code,
         VLfirstpos = vln,
         IC80 = gmt80ls,
         group = paste(Trt, IC80le1, statuswk80),
         groupplot = case_when(group=="C3 1 1" ~ "Placebo\nIC80 < 1\nprimary",
                               group=="C3 0 1" ~ "Placebo\nIC80 > 1\nprimary",
                               group=="T1 1 1" ~ "VRC01\n10mg/kg\nIC80 < 1\nprimary",
                               group=="T1 0 1" ~ "VRC01\n10mg/kg\nIC80 > 1\nprimary",
                               group=="T2 1 1" ~ "VRC01\n30mg/kg\nIC80 < 1\nprimary",
                               group=="T2 0 1" ~ "VRC01\n30mg/kg\nIC80 > 1\nprimary",
                               group=="C3 1 0" ~ "Placebo\nIC80 < 1\nnon-primary",
                               group=="C3 0 0" ~ "Placebo\nIC80 > 1\nnon-primary",
                               group=="T1 1 0" ~ "VRC01\n10mg/kg\nIC80 < 1\nnon-primary",
                               group=="T1 0 0" ~ "VRC01\n10mg/kg\nIC80 > 1\nnon-primary",
                               group=="T2 1 0" ~ "VRC01\n30mg/kg\nIC80 < 1\nnon-primary",
                               group=="T2 0 0" ~ "VRC01\n30mg/kg\nIC80 > 1\nnon-primary"),
         groupplot = factor(groupplot, levels=c("Placebo\nIC80 < 1\nprimary", "VRC01\n10mg/kg\nIC80 < 1\nprimary", "VRC01\n30mg/kg\nIC80 < 1\nprimary",
                                                "Placebo\nIC80 > 1\nprimary", "VRC01\n10mg/kg\nIC80 > 1\nprimary", "VRC01\n30mg/kg\nIC80 > 1\nprimary",
                                                "Placebo\nIC80 < 1\nnon-primary", "VRC01\n10mg/kg\nIC80 < 1\nnon-primary", "VRC01\n30mg/kg\nIC80 < 1\nnon-primary",
                                                "Placebo\nIC80 > 1\nnon-primary", "VRC01\n10mg/kg\nIC80 > 1\nnon-primary", "VRC01\n30mg/kg\nIC80 > 1\nnon-primary")),
         boxcol    = case_when(group=="C3 1 1" ~ "blue",
                               group=="C3 0 1" ~ "blue",
                               group=="T1 1 1" ~ "red3",
                               group=="T1 0 1" ~ "red3",
                               group=="T2 1 1" ~ "red3",
                               group=="T2 0 1" ~ "red3",
                               group=="C3 1 0" ~ "blue",
                               group=="C3 0 0" ~ "blue",
                               group=="T1 1 0" ~ "red3",
                               group=="T1 0 0" ~ "red3",
                               group=="T2 1 0" ~ "red3",
                               group=="T2 0 0" ~ "red3"),
         vlplot = ifelse(VLfirstpos < 40, 20, VLfirstpos))

pdat %>% group_by(group, groupplot) %>% tally()
pdat %>% group_by(groupplot) %>% tally()
pdat %>% group_by(group, groupplot) %>% summarise(minIC80=min(IC80),
                                                 maxIC80=max(IC80))


# get post week 80 cases
dat_postwk80 = read.csv(dataSurv)
pdat$status_postwk80 = ifelse(pdat$pub_id %in% subset(dat_postwk80, status_postwk80==1)$pub_id, 1, 0)
table(pdat$status_postwk80, pdat$statuswk80)

pdat = pdat %>% 
  filter(status_postwk80==1) %>%
  mutate(groupplot = 
           droplevels( 
             factor(groupplot,
             levels = levels(groupplot),
             labels= gsub('\\nnon-primary',  '', levels(groupplot)))))
  


counts = pdat %>%
  group_by(groupplot) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

# get values for manuscript for the 4 participants with IC80s < 1
tmp = subset(pdat, least_sensitive=="IC80ls <= 1", select=c('Trt','VLfirstpos'))
tmp$log10_vl = log10(tmp$VLfirstpos)
print(tmp)

# run linear regression adjusting for region 
# pull in sieve dataset for region information
sieve_dat = read.csv(dataSieve)
sieve_dat = sieve_dat %>% select(pub_id, protocol, southAmerica, southAfrica)
adat = pdat %>% 
  filter(least_sensitive == "IC80ls > 1" ) %>%
  select(pub_id, tx, VLfirstpos) %>%
  merge(sieve_dat) %>%
  mutate( region = case_when(
    protocol == 'HVTN 703' & southAfrica == 1 ~ 'RSA',
    protocol == 'HVTN 703' & southAfrica == 0 ~ 'notRSA',
    protocol == 'HVTN 704' & southAmerica == 1 ~ 'SA',
    protocol == 'HVTN 704' & southAmerica == 0 ~ 'notSA',
    TRUE ~ NA
  )) %>%
  mutate( region = factor(region),
          tx_pool = factor(ifelse(tx=='Control', "Placebo", "VRC01")) )

# sanity check for number of cases IC80ls > 1
stopifnot(nrow(adat)==31)

fit = lm(log10(VLfirstpos) ~ region + tx_pool, data=adat )
summary(fit)



set.seed(39573056)
p2 = ggplot(data=pdat, aes(x=groupplot, y=vlplot, color=boxcol)) +
  geom_violin(drop=FALSE) +
  geom_boxplot(aes(fill=boxcol), width=0.15, lwd=1, alpha = 0.3, outlier.shape=NA) + #remove outlier points
  geom_point(aes(shape=factor(status_postwk80)), size=3, position = position_jitter(w = 0.3, h = 0)) +
  geom_text(
    data = counts,
    aes(
      x = groupplot,
      y = 20000000,        # place at top of each box
      label = label
    ),
    vjust = -0.5,
    size = 5,
    col='black'
  ) +
  scale_y_continuous(limits=c(20,30000000), 
                     breaks=c(20, 100, 1000, 10000, 100000, 1000000, 10000000), 
                     labels=c("< 40", "100", expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression("">10^7)), 
                     trans="log10") +
  scale_color_manual(values = c("blue", "red3")) +
  scale_fill_manual(values = c("blue", "red3")) +
  labs(x=' ', y="Viral Load (copies/ml)", title=NULL) +
  theme_bw() +
  theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(margin = margin(t = 0, r=0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r=-10, b = 0, l = 0)))


pdf(file=pdfFileSave[1], width=17, height=7)
print(p2)
dev.off()

pdat$tx = factor(pdat$tx, 
                 levels=c("Control", "10 mg/kg", "30 mg/kg"),
                 labels=c("Placebo", "10 mg/kg", "30 mg/kg"))
pdat$protocol= factor(pdat$protocol,
                      levels=c('HVTN 703', 'HVTN 704'),
                      labels=c('HVTN 703/HPTN 081', 'HVTN 704/HPTN 085'))

p3 = ggplot(data=pdat, aes(x=gmt80ls, y=vlplot, shape=protocol, color=boxcol)) +
  geom_point() +
  scale_x_continuous(limits=c(0.1,100),
                     breaks=c(0.1, 1, 10, 100), 
                     labels=c('0.1','1','10','>100'), 
                     trans="log10") +
  scale_y_continuous(limits=c(20,10000000), 
                     breaks=c(20, 100, 1000, 10000, 100000, 1000000, 10000000), 
                     labels=c("< 40", "100", expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression("">10^7)), 
                     trans="log10") +
  scale_color_manual(values = c("blue", "red3")) +
  scale_fill_manual(values = c("blue", "red3")) +
  labs(shape=' ',
       x='IC80', 
       y="Viral Load (copies/ml)", 
       title="Week 80-104 cases: First Viral Load vs. IC80") +
  theme_bw() +
  guides(color = "none") +
  theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1, size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(margin = margin(t = 0, r=0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r=-10, b = 0, l = 0))) +
  facet_wrap(~tx, nrow=1)

pdf(file=pdfFileSave[2], width=11, height=5)
print(p3)
dev.off()


q(save='no')
