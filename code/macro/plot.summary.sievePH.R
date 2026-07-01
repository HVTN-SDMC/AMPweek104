library(ggpubr) # for combining two panels into one
library(ggpmisc) # for adding a table to the top panel
library(scales) # for pretty x-axis breaks
library(patchwork)

#' Sieve PH plot with two panels: a descriptive panel and a vaccine efficacy panel
#' @param fit.te a list of two data frames, each formatted as the output from summary.sievePH and used for
#' plotting a single VE curve
#' @param data1, a data frame specifying a continuous mark under column "markValue" and a categorical variable describing treatment (Vaccine/Placebo).
#' For subjects with a right-censored time-to-event, the value(s) in mark should be set to NA.
#' @param data2, a data frame specifying a continuous mark under column "markValue" and a categorical variable describing treatment (Vaccine/Placebo).
#' For subjects with a right-censored time-to-event, the value(s) in mark should be set to NA.
#' @param limits.x the range of the mark to be plotted
#' @param breaks.x the values where ticks will be placed for the mark
#' @param labels.x x-axis tick mark labels. If left unspecified, then values in \code{breaks.x} will be used.
#' @param xlab the x-axis label
#' @param sieveTestLabel a character string describing sieve test p-values
#' @param title title for the plot
ggplotSievePH <- function(fit.te, data1, data2, limits.x=NULL, breaks.x=NULL, labels.x=breaks.x, 
                          xlab, limits.y=NULL, sieveTestLabel=NULL, title=NULL){
  if (is.null(limits.x)){
    limits.x <- range(c(data1$mark, data2$mark), na.rm=TRUE)
  }
  
  if (is.null(limits.y)){
    limits.y <- c(-0.6, 1)
  }
  
  legend.x.just.p1 <- 0:1
 
  
  p1 <- ggplot() +
    geom_hline(yintercept=0, colour="gray50")
  
  for (i in 1:length(fit.te)){
    fit.te[[i]]$te$cohort <- factor(fit.te[[i]]$cohort, levels=c("Week 80", "Post-Week 80"))
    
    p1 <- p1 + 
      geom_ribbon(aes(x=mark, ymin=LB, ymax=UB, fill=cohort), data=fit.te[[i]]$te, alpha=0.2) +
      geom_line(aes(x=mark, y=TE, colour=cohort), data=fit.te[[i]]$te, na.rm=TRUE, size=1.9)
  }
  
  
  #browser()
  legendLabels <- c( "Weeks 0 to 80", "Weeks 80 to 104")
  vjustbottom <- ifelse(grepl("\n",xlab), -5, -0.5)
  legendpos <- ifelse(grepl("\n", xlab), -0.3, -0.25)
  breaks.y <- round(seq(-0.6, 1, 0.2), 1)
  p1 <- p1 +
    scale_colour_manual(values=c("blue4", "darkgreen"), name="", labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
    scale_fill_manual(values=c("blue4", "darkgreen"), name="", labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
    scale_y_continuous(name="Prevention Efficacy (%)", breaks=breaks.y, labels=breaks.y * 100) +
    coord_cartesian(ylim=limits.y) +
    xlab(xlab) +
    theme_bw() +
    theme(legend.key.size = unit(0.65, "cm"),
          legend.title=element_text(size = 15),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=14),
          legend.position = c(0.5,  legendpos),
          #legend.justification = c(0.5, 1),
          legend.key = element_blank(),
          legend.key.width = unit(1.4,"cm"),
          legend.background = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 2, size = 12),
          axis.title.x = element_text(size = 17, vjust= vjustbottom),
          axis.title.y = element_text(size = 17, hjust = 0.5),
          axis.text.x = element_text(size = 15, hjust = 0.5,vjust =  0.5, colour = "black"),
          axis.text.y = element_text(size = 15, colour = "black", hjust = 0.5),
          plot.margin=unit(c(0,0.2,1.5,0.2), "cm")) #t,r,b,l
  
  
  if (is.null(breaks.x)){
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = pretty_breaks(n = 5))
  } else {
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = breaks.x, labels=labels.x)
  }
 # browser()
 #to be edited
  
  data <- rbind(data1, data2) %>% 
    mutate(cohort=factor(c(rep("Week 80", length(data1$tx)), 
                           rep("Post-Week 80", length(data2$tx))), levels=c("Week 80", "Post-Week 80")))
           
  data$txchort <- factor(interaction(data$treatment, data$cohort), 
                                  levels=c(  "Placebo.Post-Week 80", "VRC01.Post-Week 80", "Placebo.Week 80","VRC01.Week 80"))
 
  p2 <- ggplot(data) +
    ggtitle(title) +
    geom_boxplot(aes(x=txchort, y=mark, colour=cohort), outlier.shape=NA, na.rm=TRUE) +
    geom_jitter(aes(x=txchort, y=mark, colour=cohort), na.rm=TRUE, width=0.2, height=0, shape=21, fill="white", size=2, stroke=1.3) +
    scale_colour_manual(values=c( "blue4", "darkgreen"), guide="none") +
    coord_flip() +
    scale_x_discrete(name=NULL, labels=c("VRC01.Week 80"="VRC01", "Placebo.Week 80"="Placebo", 
                                         "VRC01.Post-Week 80"="VRC01", "Placebo.Post-Week 80"="Placebo")) +
    
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=15),
          axis.text.y=element_text(size=12, colour="black"),
          axis.ticks.x=element_blank(),
          plot.margin=unit(c(0.2,0.2,0,0.2), "cm")) 
  
  if (is.null(breaks.x)){
    p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = pretty_breaks(n = 5), labels=NULL)
  } else {
    p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = breaks.x, labels=NULL)
  }

  # p <- ggarrange(p2, p1, heights=c(1, 2), ncol=1, nrow=2, align="v")
  # 
  # if (!is.null(title)){
  #   p <- annotate_figure(p, top=text_grob(title, color="black", face="bold", size=14, vjust=2.3))  
  # }
  
  p <- p2 + p1 + plot_layout(ncol=1, heights=c(1, 2))
  
  return(p)
}

