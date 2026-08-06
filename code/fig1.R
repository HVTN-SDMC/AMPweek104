#Predicted PT80 over time for VRC01 and placebo arms in HVTN 703/HPTN 081 and HVTN 704/HPTN 085
rm(list=ls(all=TRUE))

library(here)
here::i_am("README.md")
repoDir <- here::here()
datDir <- file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")
outDir <- file.path(repoDir, "output")

library(tidyverse)
library(DescTools)
library(plyr)

#read in week 80 visit times for all noncases (82 pubids) in the MITT Case-Control Cohort
week80visits <- read.csv(file.path(datDir, "dat_ctl_wk80_time.csv"))[, c(1,2)]
#rename column `time` to week80VisitTime
names(week80visits)[names(week80visits) == "time"] <- "week80VisitTime"

#set up tibbles for estimated xx% of exposures with PT80 < 10, >100 over week 0-xx
exposure_perc <- tibble("Protocol" = character(), "Follow-up" = character(), "Dose" = character(), "PT80_category" = character(), "Perc" = numeric())
medianPT80 <- tibble("Protocol" = character(), "Follow-up" = character(), "Dose" = character(), "Median PT80" = numeric())
geometricMeanPT80 <- tibble("Protocol" = character(), "Follow-up" = character(), "Dose" = character(), "GM PT80" = numeric())

for(week in c(80, 104)){
  if(week == 104){
    #read in sieve mark file
    sievedata <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9_wk104_v2.csv")) %>%
      # two 703 ppts with missing sequences have also a missing time-to-event
      filter(!is.na(hiv1fpday)) %>%
      # stratification variable 
      mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                                protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                                protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                                protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
    #read in VRC01 serum concentration trajectories for all noncases (82 pubids) in the MITT Case-Control Cohort
    dat_est_daily_ctl <- read.csv(file.path(datDir, "dat_est_daily_82Ctl_wk112.csv"))
    dat_est_daily_ctl <- left_join(dat_est_daily_ctl, week80visits, by = "pub_id")
    
   
  }else{
    #read in sieve mark file
    sievedata <- read.csv(file.path(datDir, "amp_sieve_pooled_marks_final_v9c.csv")) %>%
      # two 703 ppts with missing sequences have also a missing time-to-event
      filter(!is.na(hiv1fpday)) %>%
      # stratification variable 
      mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                                protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                                protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                                protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
    #read in VRC01 serum concentration trajectories for all noncases (82 pubids) in the MITT Case-Control Cohort
    dat_est_daily_ctl <- read.csv(file.path(datDir, "dat_est_daily_82Ctl_wk112.csv"))
    dat_est_daily_ctl <- left_join(dat_est_daily_ctl, week80visits, by = "pub_id")
    dat_est_daily_ctl <- filter(dat_est_daily_ctl, time <= week80VisitTime)
    
  }
 
  #pre-process gmt80ls
  sievedata$gmt80ls[sievedata$gmt80ls == ">100"] <- "200" #impute >100 with 200
  sievedata$gmt80ls <- as.numeric(sievedata$gmt80ls)
  dat_ic80_raw <- sievedata[, c("protocol", "pub_id", "tx", "gmt80ls")]
  dat_ic80 <- dat_ic80_raw%>%
    dplyr::select(protocol,tx,pub_id, gmt80ls) %>%
    mutate(protocol=ifelse(protocol=="HVTN 703","HVTN 703/HPTN 081","HVTN 704/HPTN 085")
           ,rx_code=ifelse(tx=="T1","VRC01 10 mg/kg",ifelse(tx=="T2","VRC01 30 mg/kg","Control"))) 
  colnames(dat_ic80) <- c("study", "tx", "PUB_ID", "ic80", "dose")
  
  dat_ic80_placebo <- dat_ic80 %>%filter(dose=="Control", !is.na(ic80))
  dat_ic80_placebo <- dplyr::select(dat_ic80_placebo, all_of(c("study", "ic80")))
  
  #add dose information
  dat_est_daily_ctl <- left_join(dat_est_daily_ctl, sievedata[, c("pub_id", "tx")], by = "pub_id")
  dat_est_daily_ctl$tx <- ifelse(dat_est_daily_ctl$tx=="T1","VRC01 low dose", "VRC01 high dose")
  
  
  
  ## calculate daily-grid ID80 by conc/IC80
  pubid_704 <- filter(sievedata, protocol == "HVTN 704")$pub_id
  pubid_703 <- filter(sievedata, protocol == "HVTN 703")$pub_id
  
  dat_ID80_daily_ctl_704 <- filter(dat_est_daily_ctl, pub_id %in% pubid_704) %>%
    cross_join(filter(dat_ic80_placebo, study == "HVTN 704/HPTN 085")) %>%
    mutate(ID80=cc/ic80)
  
  dat_ID80_daily_ctl_703 <- filter(dat_est_daily_ctl, pub_id %in% pubid_703) %>%
    cross_join(filter(dat_ic80_placebo, study == "HVTN 703/HPTN 081")) %>%
    mutate(ID80=cc/ic80)
  
  dat_ID80_daily_ctl <- rbind(dat_ID80_daily_ctl_704, dat_ID80_daily_ctl_703)
  
  #An estimated xx% of exposures had PT80 > 41.3 that was connected to 75% prevention efficacy.
  print(sum(dat_ID80_daily_ctl$ID80> 41.3)/dim(dat_ID80_daily_ctl)[1])
  
  if(week == 104){
    for(f in c("Weeks 0-80", "Weeks 80-104")){
      for(protocol in c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081")){
        for(d in c("VRC01 low dose", "VRC01 high dose")){
          if(f == "Weeks 0-80" ){
            dat_sub <- filter(dat_ID80_daily_ctl, time <= week80VisitTime & study == protocol & tx == d)
          }else{
            dat_sub <- filter(dat_ID80_daily_ctl, time > week80VisitTime & study == protocol & tx == d)
          }
          exposure_perc <- add_row(.data = exposure_perc, "Protocol" = protocol, "Follow-up" = f, 
                                   "Dose" = d, "PT80_category" = "PT80 < 10", 
                                   "Perc" = sum( dat_sub$ID80 < 10)/nrow(dat_sub))
          exposure_perc <- add_row(.data = exposure_perc, "Protocol" = protocol, "Follow-up" = f,
                                   "Dose" = d, "PT80_category" = "PT80 >100", 
                                   "Perc" = sum( dat_sub$ID80 >100)/nrow(dat_sub))
          medianPT80 <- add_row(.data = medianPT80, "Protocol" = protocol, "Follow-up" = f, 
                            "Dose" = d, "Median PT80" = median(dat_sub$ID80))
          
          geometricMeanPT80 <- add_row(.data = geometricMeanPT80, "Protocol" = protocol, "Follow-up" = f, 
                                       "Dose" = d, "GM PT80" = 10^mean(log10(dat_sub$ID80)))
          
        }
      }}
  }
  
  
  
  #daily median and quantiles of ID80 across participants and viruses
  dat_ID80_daily_ctl_summary <- ddply(dat_ID80_daily_ctl, .(study, tx, time), function(df)quantile(df$ID80, probs = c(0.025, 0.5, 0.975)))
  colnames(dat_ID80_daily_ctl_summary) <- c("study", "tx", "time", "ll", "median", "ul")
  dat_ID80_daily_ctl_summary$group <- paste(dat_ID80_daily_ctl_summary$study, "\n",dat_ID80_daily_ctl_summary$tx, sep = "")
  
  dat_ID80_daily_ctl_summary$group <- factor(dat_ID80_daily_ctl_summary$group, 
                                             levels = c("HVTN 704/HPTN 085\nVRC01 low dose",
                                                        "HVTN 704/HPTN 085\nVRC01 high dose",
                                                        "HVTN 703/HPTN 081\nVRC01 low dose",
                                                        "HVTN 703/HPTN 081\nVRC01 high dose"))
  prettyNum0 <- function(x){sprintf("%.5g", x)}
  
  dat_ID80_daily_ctl_summary$median <- ifelse(dat_ID80_daily_ctl_summary$median <=0.01, 0.01, dat_ID80_daily_ctl_summary$median)
  dat_ID80_daily_ctl_summary$ll <- ifelse(dat_ID80_daily_ctl_summary$ll <=0.01, 0.01, dat_ID80_daily_ctl_summary$ll)
  dat_ID80_daily_ctl_summary$ul <- ifelse(dat_ID80_daily_ctl_summary$ul <=0.01, 0.01, dat_ID80_daily_ctl_summary$ul)
  
  dat_ID80_daily_ctl_summary$study <- factor(dat_ID80_daily_ctl_summary$study, levels = c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081"))
  dat_ID80_daily_ctl_summary$tx <- factor(dat_ID80_daily_ctl_summary$tx, levels = c("VRC01 low dose", "VRC01 high dose"))
  
  if(week == 104){
    breaks.x = seq(0, 112, 8)*7
    labels.x = seq(0, 112, 8)
    plotxlimit = 112*7
    annotatex = 116*7
    y0 = 10
    y50 = 32
    y90 = 194
  }else{
    breaks.x = seq(0, 88, 8)*7
    labels.x = seq(0, 88, 8)
    plotxlimit = 88*7 
    annotatex = 91*7+3
    y0 = 10
    y50 = 32
    y90 = 194
  }
  
  for(tx.x in c("VRC01 low dose", "VRC01 high dose")){
    for(study.x in c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081")){
      title <- case_when(tx.x == "VRC01 low dose" & study.x == "HVTN 704/HPTN 085" ~ "HVTN 704/HPTN 085 Americas trial: VRC01 low dose",
                         tx.x == "VRC01 high dose" & study.x == "HVTN 704/HPTN 085" ~ "HVTN 704/HPTN 085 Americas trial: VRC01 high dose",
                         tx.x == "VRC01 low dose" & study.x == "HVTN 703/HPTN 081" ~ "HVTN 703/HPTN 081 Africa trial: VRC01 low dose",
                         tx.x == "VRC01 high dose" & study.x == "HVTN 703/HPTN 081" ~ "HVTN 703/HPTN 081 Africa trial: VRC01 high dose"
      )
      xlimit = case_when(week == 104  ~ 112*7,
                         week == 80 & study.x == "HVTN 703/HPTN 081"  ~ 85.9*7,
                         week == 80 & study.x == "HVTN 704/HPTN 085"  ~ 87.0*7)
      
    
      p <- ggplot(data = filter(dat_ID80_daily_ctl_summary,
                                tx == tx.x & study == study.x & time <= xlimit)) +
        geom_line(aes(x = time, y = median, color = "median"), linewidth = 0.6) +
        geom_point(aes(x = time, y = median, color = "median"), size = 0.6) +
        geom_line(aes(x = time, y = ll, color = "ll"), linewidth = 0.6) +
        geom_point(aes(x = time, y = ll, color = "ll"), size = 0.6 ) +
        geom_line(aes(x = time, y = ul, color = "ul"), linewidth = 0.6) +
        geom_point(aes(x = time, y = ul, color = "ul"), size = 0.6) +
        geom_hline(yintercept = 10, linewidth = 0.7, color = "darkorange") + 
        geom_hline(yintercept = 32, linewidth = 0.7, color = "darkorange") + 
        geom_hline(yintercept = 194, linewidth = 0.7, color = "darkorange") + 
        annotate(geom = "text", x = annotatex, y = y0, label = "PE = 0%", hjust = 0, size = 3.5) + 
        annotate(geom = "text", x = annotatex, y = y50, label = "PE = 50%", hjust = 0, size = 3.5) + 
        annotate(geom = "text", x = annotatex, y = y90, label = "PE = 90%", hjust = 0, size = 3.5) + 
        scale_color_manual(name = "", breaks = c("median", "ll", "ul"), values = c("black", "blue", "blue"), 
                           labels = c("Median", "2.5th percentile", "97.5th percentile"))+
        guides(color = guide_legend(override.aes = list(size = 1.3, linewidth = 0.7, alpha = 1)))+
        ggtitle(title) +
        scale_y_log10(
          breaks = c(0.01,0.1,1,10,100,1000), minor_breaks = NULL,
          limits=c(min( dat_ID80_daily_ctl_summary$ll),max( dat_ID80_daily_ctl_summary$ul)),
          labels= c(expression(""<="0.01"), 0.1,1,10,100,1000))+
        scale_x_continuous(breaks = breaks.x, labels = labels.x, 
                           minor_breaks = NULL,
                           expand = expansion(mult = 0.03, add = 0.03))+
        coord_cartesian(xlim = c(0, plotxlimit), clip = "off")+
        theme_bw(base_size = 8)+
        ylab(expression(paste("Predicted VRC01 serum ", ID[80]," titers (", PT[80], ") ",  
                              "against placebo virus isolates")))+
        xlab("Weeks since first infusion")+
        theme(plot.margin = unit(c(0.25,0.7,0.25,0.25), "in"), # t, r, b, l
              legend.position = "bottom",
              legend.background = element_blank(),
              legend.key.width = unit(0.75, "cm"),
              plot.title = element_text(hjust = 0.5, size = 13),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              strip.text.x = element_text(size=12),
              strip.text.y = element_text(size=12),
              axis.title.y = element_text(margin = margin(r=6), size = 12),
              axis.title.x = element_text(margin = margin(t=8), size = 12))
      
      ggsave(file.path(figDir, paste0(ifelse(week==104, "Predicted_PT80_wk104_","Predicted_PT80_wk80_" ),
                                      gsub(" ", "_", tx.x), "_", ifelse(study.x=="HVTN 704/HPTN 085", "704", "703") , ".pdf")),
             plot=p, width=ifelse(week == 104, 7, 7), height=7)
      
    }
  }
  
  #####################################################################################
  #a single plot for pooled trials
  dat_ID80_daily_ctl <- dat_est_daily_ctl %>%cross_join(dat_ic80_placebo) %>%
    mutate(ID80=cc/ic80)
  #add median PT80 for trials pooled and dose pooled
  if(week == 104){
    for(f in c("Weeks 0-80", "Weeks 80-104")){
      if(f == "Weeks 0-80" ){
        dat_sub <- filter(dat_ID80_daily_ctl, time <= week80VisitTime)
      }else{
        dat_sub <- filter(dat_ID80_daily_ctl, time > week80VisitTime)
      }
      medianPT80 <- add_row(.data = medianPT80, "Protocol" = "Trials pooled", "Follow-up" = f, 
                            "Dose" = "VRC01 dose pooled", "Median PT80" = median(dat_sub$ID80))
      
    }
  }
  
  #daily median and quantiles of ID80 across participants and viruses
  dat_ID80_daily_ctl_summary <- ddply(dat_ID80_daily_ctl, .(time), function(df)quantile(df$ID80, probs = c(0.025, 0.5, 0.975)))
  colnames(dat_ID80_daily_ctl_summary) <- c("time", "ll", "median", "ul")
  
  prettyNum0 <- function(x){sprintf("%.5g", x)}
  
  dat_ID80_daily_ctl_summary$median <- ifelse(dat_ID80_daily_ctl_summary$median <=0.01, 0.01, dat_ID80_daily_ctl_summary$median)
  dat_ID80_daily_ctl_summary$ll <- ifelse(dat_ID80_daily_ctl_summary$ll <=0.01, 0.01, dat_ID80_daily_ctl_summary$ll)
  dat_ID80_daily_ctl_summary$ul <- ifelse(dat_ID80_daily_ctl_summary$ul <=0.01, 0.01, dat_ID80_daily_ctl_summary$ul)
  
  xlimit = case_when(week == 104 ~ 112*7,
                     week == 80  ~ 85.9*7)
  
  p <- ggplot(data = filter(dat_ID80_daily_ctl_summary, time <= xlimit)) +
    geom_line(aes(x = time, y = median, color = "median"), linewidth = 0.6) +
    geom_point(aes(x = time, y = median, color = "median"), size = 0.6) +
    geom_line(aes(x = time, y = ll, color = "ll"), linewidth = 0.6) +
    geom_point(aes(x = time, y = ll, color = "ll"), size = 0.6 ) +
    geom_line(aes(x = time, y = ul, color = "ul"), linewidth = 0.6) +
    geom_point(aes(x = time, y = ul, color = "ul"), size = 0.6) +
    geom_hline(yintercept = 10, linewidth = 0.7, color = "darkorange") + 
    geom_hline(yintercept = 32, linewidth = 0.7, color = "darkorange") + 
    geom_hline(yintercept = 194, linewidth = 0.7, color = "darkorange") + 
    annotate(geom = "text", x = annotatex, y = y0, label = "PE = 0%", hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = annotatex, y = y50, label = "PE = 50%", hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = annotatex, y = y90, label = "PE = 90%", hjust = 0, size = 3.5) + 
    scale_color_manual(name = "", breaks = c("median", "ll", "ul"), values = c("black", "blue", "blue"), 
                       labels = c("Median", "2.5th percentile", "97.5th percentile"))+
    guides(color = guide_legend(override.aes = list(size = 1.3, linewidth = 0.7, alpha = 1)))+
    scale_y_log10(
      breaks = c(0.01,0.1,1,10,100,1000), minor_breaks = NULL,
      limits=c(min( dat_ID80_daily_ctl_summary$ll),max( dat_ID80_daily_ctl_summary$ul)),
      labels= c(expression(""<="0.01"), 0.1,1,10,100,1000))+
    scale_x_continuous(breaks = breaks.x, labels = labels.x, 
                       minor_breaks = NULL,
                       expand = expansion(mult = 0.03, add = 0.03))+
    coord_cartesian(xlim = c(0, plotxlimit), clip = "off")+
    theme_bw(base_size = 8)+
    ylab(expression(paste("Predicted VRC01 serum ", ID[80]," titers (", PT[80], ") ",  
                          "against placebo virus isolates")))+
    xlab("Weeks since First Infusion")+
    theme(plot.margin = unit(c(0.25,0.75,0.25,0.25), "in"),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key.width = unit(0.75, "cm"),
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=12),
          axis.title.y = element_text(margin = margin(r=6), size = 12),
          axis.title.x = element_text(margin = margin(t=8), size = 12))
  
  ggsave(file.path(figDir, paste0(ifelse(week==104, "Predicted_PT80_wk104_","Predicted_PT80_wk80_" ), 
                                  "trials_tx_pooled" , ".pdf")),
         plot=p, width=ifelse(week == 104, 7, 7), height=7)
  
  
  
  
}

exposure_perc$Perc <- 100*round(exposure_perc$Perc, 3)
medianPT80$`Median PT80` <- round(medianPT80$`Median PT80`, 3)
geometricMeanPT80$`GM PT80` <- round(geometricMeanPT80$`GM PT80`, 3)

write.csv(exposure_perc, file.path(tabDir, "Predicted_PT80_exposure_perc.csv"), row.names = FALSE)
write.csv(medianPT80, file.path(tabDir, "Predicted_PT80_median.csv"), row.names = FALSE)

    