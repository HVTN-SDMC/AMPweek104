library(tidyverse)
library(lubridate)

library(here)
here::i_am("README.md")
repoDir <- here::here()
macroDir <- file.path(repoDir, "code/macro")
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")

# data files
dataFile704 = file.path(dat2Dir, 'v704_survival_wk80_tau_neut.csv')
dataDX = file.path(dat2Dir, 'amp_diagnostic_trajectories.csv')
dataIDT704 = file.path(dat2Dir, 'v704_infusion_dates.csv')
dataDBS1 = file.path(dat2Dir, 'VTN704_case_control_DBS01UC20200715E001.txt')
dataDBS2 = file.path(dat2Dir, 'VTN704_DBS01_20210521.txt')
dataDBS3 = file.path(dat2Dir, 'VTN704_DBS01_PILOT_UC_20181114.txt')

pdfFileSave = file.path(figDir, 'dbs_704primary_cases.pdf')
tabFileSave = file.path(tabDir, 'dbs_prep_use_estimates_704cases_week0to80.csv')

dat = read.csv(dataFile704)
dat$southAmerica = NULL
dat = subset(dat, select=c('pub_id','hiv1event','protocol','tx'))

case_ids = subset(dat, protocol=='HVTN 704' & hiv1event==1)$pub_id

dx = read.csv(dataDX)
dx = dx %>% select(pub_id, dxdy) %>% unique()

# infusion dates
idt = read.csv(dataIDT704)
idt$inf.number = as.numeric(factor(idt$visit))

# compute number of days between infusion and enrollment (ie, first infusion)
idt = idt %>% 
  group_by(pub_id) %>% 
  arrange(idt) %>%
  mutate(ndays=as.numeric(idt - first(idt)),
         ninf = 1:n())
idt = idt[order(idt$pub_id, idt$ndays),]

EFFECTIVE_CUTOFF = 700
LLOQ="<<"

# Data split across three datasets (combine)
dbs = read.csv(dataDBS1, sep='\t')
dbs_monitor = read.csv(dataDBS2, sep='\t')
dbs_pilot = read.csv(dataDBS3, sep='\t')
dbs = rbind(dbs, dbs_monitor, dbs_pilot)

dbs = dbs %>% 
  select("pub_id","visitno","drawdy", "concentration","concentration_units","analyte","concentration_oor_indicator") %>%
  filter(analyte=="TFV-DP")
dbs$detectable = as.numeric(!is.na(dbs$concentration) & dbs$concentration_oor_indicator!=LLOQ)
dbs$effective  = as.numeric(!is.na(dbs$concentration) & dbs$concentration >= EFFECTIVE_CUTOFF)
dbs$conc = ifelse(is.na(dbs$concentration), 12.5, dbs$concentration)

# check the number of primary 704 cases with dbs data 97 out of 98
sum(unique(dbs$pub_id) %in% case_ids)
dbs = dbs %>% filter(pub_id %in% case_ids)

# merge with dx
dbs = merge(dbs, dx, all.x=TRUE) %>%
  merge(dat %>% select(pub_id, tx, hiv1event)) %>%
  mutate(ndays = drawdy) %>%
  arrange(pub_id, ndays)

dbs = dbs %>%
  group_by(pub_id) %>%
  mutate(days2nearest_draw = min(abs(drawdy-dxdy))) %>%
  ungroup() %>%
  arrange(pub_id, ndays)

# conveniently all cases with evidence of dbs use have a diagnosis date and draw date that coincide
# which makes it simple to plot
stopifnot(all((dbs %>% filter(days2nearest_draw !=0))$conc == 12.5))
dbs$diagnosis_time = ifelse(dbs$drawdy==dbs$dxdy, 1, 0)

# most cases had no evidence of any use of PrEP
dbs_summ = dbs %>%
  group_by(pub_id, tx) %>%
  summarise(any_use=any(!is.na(concentration)),
            .groups = 'drop')
table(dbs_summ$any_use, dbs_summ$tx)

# conveniently all cases with evidence of dbs use have a diagnosis date and draw date that coincide
# which makes it simple to plot
stopifnot(all((dbs %>% filter(days2nearest_draw !=0))$conc == 12.5))
dbs$diagnosis_time = ifelse(dbs$drawdy==dbs$dxdy, 1, 0)


# for participants with no draw at diagnosis combine the last draw
# and the diagnosis draw data
# this data will be used to draw dashed lines for these ppts
dbs_nodraw = dbs %>%
  filter(days2nearest_draw != 0) %>%
  group_by(pub_id, tx) %>%
  summarise(
    i = which.min(abs(drawdy-dxdy)),
    ndays = dxdy[i],
    conc = conc[i],
    .groups = 'drop')
dbs_last_draw = lapply(1:nrow(dbs_nodraw), function(row) {
  id = dbs_nodraw$pub_id[row]
  i = dbs_nodraw$i[row]
  
  dbs %>%
    filter(pub_id==id) %>%
    mutate(row = row_number(),
           last_draw=TRUE) %>%
    filter(row == i) %>%
    select(pub_id, tx, ndays, conc, last_draw)
})
dbs_last_draw = do.call(rbind, dbs_last_draw)
dbs_nodraw = dbs_nodraw %>%
  mutate(last_draw=FALSE) %>%
  select(-i) %>%
  rbind(dbs_last_draw) %>%
  arrange(pub_id, ndays)

# Calculate "step" values for each participant
# where ppts with no detectible concentrations will have steps below 33 that
# can be used to assign concentrations to a STEP value in the plot
dbs_step = dbs %>%
  group_by(pub_id, tx, dxdy) %>%
  summarise(mconc = max(conc),
            max_days = max(ndays),
            any_use=any(!is.na(concentration)),
            .groups = 'drop') %>%
  mutate( days2dx = dxdy ) %>%
  group_by(pub_id, tx, mconc, any_use, days2dx) %>%
  summarise(total_max_days = max(max_days, days2dx),
            .groups = 'drop') %>%
  arrange(-mconc, total_max_days) %>%
  group_by(tx) %>%
  mutate(step = row_number(),
         gray = rank(-step)*any_use) %>% # gray could be used to gradate the lines with some use
  arrange(tx, step)

# merge the step data with the dbs data 
NROW = nrow(dbs)
dbs = dbs %>%
  mutate(order = row_number()) %>%
  merge(dbs_step) %>%
  arrange(order) %>%
  select(-order)
stopifnot(nrow(dbs)==NROW)
rm(NROW)

# table number with any use by treatment
table(dbs_step$any_use, dbs_summ$tx)
STEP = 10^(seq(log10(20), log10(2), length.out=max(table(dbs_step$tx))))
dbs$conc_step = ifelse( dbs$conc > 12.5, dbs$conc, STEP[dbs$step])

# merge the step data with the dbs data 
NROW = nrow(dbs_nodraw)
dbs_nodraw = dbs_nodraw %>%
  mutate(order = row_number()) %>%
  merge(dbs_step) %>%
  arrange(order) %>%
  select(-order)
stopifnot(nrow(dbs_nodraw)==NROW)
rm(NROW)

dbs_nodraw$conc_step = ifelse( dbs_nodraw$conc > 12.5, dbs_nodraw$conc, STEP[dbs_nodraw$step])

p = ggplot(dbs, aes(x = ndays/7, y = conc_step, group = pub_id, fill = tx, color = tx)) +
  geom_line(color = "grey70", alpha = 0.6) +
  geom_line(data=dbs_nodraw, color = "grey70", alpha = 0.6, linetype=3, alpha = 0.6) +
  geom_point(data = subset(dbs, diagnosis_time != 1), color = "gray70", size = 1, alpha = 0.6) +
  geom_point(data = subset(dbs, diagnosis_time == 1), size = 1) +
  geom_point(data = dbs_nodraw %>% filter( !last_draw ), size = 1, shape=17) +
  scale_y_log10(breaks=c(25,100,1000),
                labels=c('<25', 100, 1000)) +
  geom_hline(yintercept = 700, linetype = "dashed", color = "black") +
  annotate("text", x = max(dbs$ndays/7), y = 700, 
           label = "Effective (700)", color = "black", hjust = 1, vjust = -0.5, size = 4) +
  theme_bw() +
  labs(title="HVTN 704/HPTN 085 Cases Week 0-80", x = "Weeks since Enrollment", y = "TFV-DP (fmol/punch)", 
       fill = "Treatment",
       color = "Treatment") +
  scale_color_manual(
    values = c("C3" = "blue", "T1" = "red3", "T2" = "red3"),
    labels = c("C3" = "Placebo", "T1" = "10 mg/kg", "T2" = "30 mg/kg")
  ) +
  scale_fill_manual(
    values = c("C3" = "blue", "T1" = "red3", "T2" = "red3"),
    labels = c("C3" = "Placebo", "T1" = "10 mg/kg", "T2" = "30 mg/kg")
  ) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  facet_wrap(~factor(tx, levels=c("C3", "T1", "T2"), labels=c("Placebo", "10 mg/kg", "30 mg/kg")), nrow=3, scales='free_y')



pdf(pdfFileSave, width=11, height=8.5)
print(p)
dev.off()


#### Estimate the percent person years on effective PrEP among
# cases by treatment arm using the same method use to compute for NEJM paper

set.seed(23)

out = c()
trial = "HVTN 704"

for( treatment in c("C3", "T1", "T2") ) {
  
  dat = dbs %>% 
    filter(pub_id %in% case_ids, drawdy <= dxdy, tx==treatment)
  dat$days = dat$drawdy
  dat$visit = dat$visitno
  
  
  # generate df matrix (analogous to the 'df' matrix used in ../../code/dbs_prep_sampling.R) with
  # rows for ppts and columns for visits and entries are dates in number of days since the origin 
  # date 1 Jan 1970
  df = reshape(dat[,c("pub_id","days","visit")], v.names="days", idvar="pub_id", timevar="visit", direction="wide")
  rownames(df) = df$pub_id
  colnames(df) = sub("days.","", colnames(df))
  df$pub_id = NULL
  df = df[,order(as.numeric(colnames(df)))]
  df = as.matrix(df)
  
  
  # Define the sampling indicator Delta_{ij} for these data all 1s when entry of df is not NA
  
  Dij = df
  Dij[which(!is.na(df))] = 1
  stopifnot(all(rownames(df)==rownames(Dij)))
  stopifnot(all(colnames(df)==colnames(Dij)))
  
  # Generate t0 list similar to T0 in the prep sampling simulations.  T0 is the site specific activation date for version 2
  # of protocol.  Here just the first observed drawdat for each pub_id as a number
  t0 = lapply(rownames(df), function(id) {
    ss = dat %>% filter(pub_id==id)
    min(ss$drawdy)
  })
  
  
  # Compute \pi_{ij}
  piij = df
  piij[which(!is.na(df))] = 1
  
  # Person-years
  Pij = t(sapply(1:nrow(df), function(i) {
    t = df[i,]
    
    w = which(!is.na(t) & t>=t0[[i]])
    
    # observed sampling times within the DBS sampling timeframe [T0,T]
    tobs = t[w]
    t0i = t0[[i]] - 1 # subtract 1 simplifies the computation of person years
    tobs = c(t0i, tobs)
    
    # person-years between observed sample collection times
    pobs = tobs[-1] - tobs[-length(tobs)]
    stopifnot(all(pobs>=0))
    
    # vector p contains NA when sample wasn't collected for ppt i at visit j (ie, df[i,j] is NA)
    p = rep(NA, ncol(df))
    p[w] = pobs
    
    return(p)
  }))
  
  # Define 'onprep' 
  for( n in c("detectable","effective") ){
    dat$onprep = dat[[n]]
    
    xij = reshape(dat[,c("pub_id","onprep","visit")], v.names="onprep", idvar="pub_id", timevar="visit", direction="wide")
    rownames(xij) = xij$pub_id
    colnames(xij) = sub("onprep.","", colnames(xij))
    xij$pub_id = NULL
    xij = as.matrix(xij)
    xij = xij[,order(as.numeric(colnames(xij)))]
    stopifnot(all(rownames(df)==rownames(xij)))
    stopifnot(all(colnames(df)==colnames(xij)))
    
    # DBS assay results
    
    # numerator and denominator of \widehat{\Phi}
    w.enrolled = 1:nrow(df)
    N.enrolled = length(w.enrolled)
    
    A = (Dij * xij * (1/piij) * Pij)[w.enrolled,]
    B = Pij[w.enrolled,]
    
    w = 1:length(A)
    phihat = sum(A[w], na.rm=TRUE)/sum(B[w], na.rm=TRUE)
    BOOT = 1000
    N.sampled = sum(Dij[w.enrolled,][w], na.rm=TRUE)
    N.pos     = sum((xij * Dij)[w.enrolled,][w], na.rm=TRUE)
    N.tot     = length(w)
    
    phihat.boot = replicate(BOOT, {
      rows = sample(N.enrolled, size=N.enrolled, replace=TRUE)
      Asamp = A[rows,]
      Bsamp = B[rows,]
      
      w = 1:length(Asamp)
      sum(Asamp[w], na.rm=TRUE)/sum(Bsamp[w], na.rm=TRUE)
    })
    
    phihat.ci = quantile(phihat.boot, probs=c(0.025, 0.975))
    
    # compute bootstrap CIs for proportions (not adjusted for person-years) 
    # results should be similar
    Dpos = (Dij * xij)[w.enrolled,]
    D    = Dij[w.enrolled,]
    
    phat = N.pos/N.sampled
    
    phat.boot = replicate(BOOT, {
      rows = sample(N.enrolled, size=N.enrolled, replace=TRUE)
      Dsamp = D[rows,]
      Dpossamp = Dpos[rows,]
      
      w = 1:length(Dsamp)
      sum(Dpossamp[w], na.rm=TRUE)/sum(Dsamp[w], na.rm=TRUE)
    })
    
    phat.ci = quantile(phat.boot, probs=c(0.025, 0.975))
    
    out = rbind(out, data.frame(Trial=trial, tx=treatment, type=n, pos=N.pos, total=N.sampled,
                                phat=phat, phat.lci=phat.ci[1], phat.uci=phat.ci[2], 
                                phihat=phihat, phihat.lci = phihat.ci[1], phihat.uci=phihat.ci[2]))    
  }
}  

out 

write.csv(out, file=tabFileSave, row.names = FALSE)

q(save='no')













