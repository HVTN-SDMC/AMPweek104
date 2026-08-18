library(tidyverse)
library(ggplot2)
library(gridExtra)

library(here)
here::i_am("README.md")
repoDir <- here::here()
figDir <- file.path(repoDir, "output/figures")

pdfFileSave = file.path(figDir, "amp_efficacy_models.pdf")

q = 0.7
p = 0.25
lambda = 0.03/365

# Get the 95% CI about an estimate of PE based on 1 - p1/p2 with sample sizes n1 and n2
pe_ci = function(p1, p2, n1, n2) {
  sd = sqrt(((1/p1) - 1)/n1 + ((1/p2) - 1)/n2)
  res = log(p1/p2) + c(-1,1)*1.96*sd
  res = exp(res)
  1 - res
}

# CDF for the sum S = T and D where T and D are exponentially distributed with lambdas VRC01 and D.
F_S = function(t, lambda_VRC01, lambda_D) {
  1 - ((lambda_VRC01 * exp(-lambda_D * t) - lambda_D * exp(-lambda_VRC01 * t)) / (lambda_VRC01 - lambda_D))
}

# Conditional distribution of S conditional on T <= 560 days (80 weeks)
F_S_cond_T = Vectorize(function(t, lambda_VRC01, lambda_D) {
  # Compute P(T <= 560)
  P_T_560 <- 1 - exp(-lambda_VRC01 * 560)
  
  if (t <= 560) {
    tmp = F_S(t, lambda_VRC01, lambda_D)
  } else {
    tmp = ((lambda_VRC01)/(lambda_D - lambda_VRC01)) * (exp(-lambda_D*(t-560)) * exp(-lambda_VRC01*560) - exp(-lambda_D*t))
    tmp =  P_T_560 - tmp
  }
  
  tmp/P_T_560
})

# Make a cartoon plot
make_cartoon_plot = function(df, title) {
  p = ggplot(df, aes(x = x)) +
    geom_line(aes(y = y1, color = "green"), linewidth = 1) +
    geom_line(aes(y = y2, color = "red"), linewidth = 1) +
    geom_line(aes(y = y3, color = "black"), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("green" = "green", "red" = "red", "black" = "black"), 
                       labels = c("Sensitive", "Resitant", "Overall"),
                       breaks = c("green", "red", "black") ) +
    labs(x = "Weeks", 
         y = "Dose-Pooled VRC01 vs. Placebo\nPrevention Efficacy (%)", 
         color='',
         title=title) +
    scale_y_continuous(
      limits = c(-1, 1),                        # Set y-axis limits
      breaks = seq(-1, 1, by = 0.5),            # Define tick positions
      labels = seq(-100, 100, by = 50)          # Label as -100 to 100
    ) +
    scale_x_continuous(
      limits = c(0, 728),                        # Set x-axis limits
      breaks = seq(0, 728, by = 56),            # Define tick positions
      labels = seq(0, 104, by = 8)         
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(angle = 0),     
      legend.position = c(0.3, 0.3),  # Move legend inside plot (adjust as needed)
      legend.background = element_rect(fill = alpha("white", 0)), # Semi-transparent background
      plot.title = element_text(size = 5),      # Reduce title size
      axis.title = element_text(size = 5),       # Reduce axis titles size
      axis.text = element_text(size = 5),        # Reduce axis text size
      legend.title = element_text(size = 5),     # Reduce legend title size
      legend.text = element_text(size = 5)       # Reduce legend text size
    )
  
  return(p)
}

# Make plot of cumulative PE with 95% CI
make_pe_ci_plot = function(df) {
  ggplot(df, aes(x = x, y = y, color = color)) +
    geom_vline(xintercept = 0, linetype = "dashed") +  # Reference line at 0
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, size = 0.5) +  # Confidence intervals
    geom_point(size = 2) +  # Point estimates
    geom_text(aes(label = label, x=-0.2), hjust = 1, vjust = -1.2, size=2, color = "black") +  # Add labels
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25),
      labels = seq(-100, 100, by = 25)
    ) +
    scale_y_continuous(limits = c(0.5, 3.5), breaks = NULL) + 
    labs(
      x = "Cumulative PE through week 104",
      y = "",
      color = "Legend"
    ) +
    scale_color_manual(values = c("green" = "green", "red" = "red", "black" = "black"), guide = "none") +  # Custom colors, no legend
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),  # Remove y-axis labels
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title = element_text(size = 8),       # Reduce axis titles size
      axis.text = element_text(size = 6),        # Reduce axis text size
      legend.title = element_text(size = 8),     # Reduce legend title size
      legend.text = element_text(size = 6)       # Reduce legend text size
    )
}




## panel A "No-delay + no-escape + no-emerge"
## lambda differs by treatment
## ratio of resistant to sensitive differs by treatment

# pvs during the infusion period (A) 
pvsA = ((p*(1-q))/(q + p - q*p)) * (1 - exp(-lambda*(q + p - q*p)*pmin(560, 1:728)))

# pvr during the infusion period (A) 
pvrA = (q/(q + p - q*p)) * (1-exp(-lambda*(q + p - q*p)*pmin(560, 1:728)))

# proportion at risk after week 80 in VRC01 arm
pvAtRisk = exp(-lambda*(q + p - q*p)*560)

# follow-up period
pvsB = pvAtRisk * (1-q) * (1-exp(-lambda*(pmax(0, (1:728)-560))))
pvrB = pvAtRisk * q * (1-exp(-lambda*(pmax(0, (1:728)-560))))

# overall
pvs = pvsA + pvsB
pvr = pvrA + pvrB

# placebos are simple
pps = (1-q) * (1 - exp(-lambda * (1:728)))
ppr = q * (1 - exp(-lambda * (1:728)))


# Model A plots
df = data.frame( x = 1:728, 
                 y1 = 1 - pvs / pps, 
                 y2 = 1 - pvr / ppr,
                 y3 = 1 - (pvs + pvr) / (pps + ppr))
pA = make_cartoon_plot(df, 'A. No-delay + No-escape + No-emerge')

cis = pe_ci(pvs[728], pps[728], 3000, 1500)
cir = pe_ci(pvr[728], ppr[728], 3000, 1500)
cio = pe_ci(pvs[728] + pvr[728], pps[728] + ppr[728], 3000, 1500)
df = data.frame(
  x = c(1 - pvs[728] / pps[728], 1 - pvr[728] / ppr[728], 1 - (pvs[728] + pvr[728]) / (pps[728] + ppr[728])),  # Point estimates
  y = c(3, 2, 1),  # Y positions
  lower = c(cis[1], cir[1], cio[1]),  # Lower bound of CI
  upper = c(cis[2], cir[2], cio[2]),  # Upper bound of CI
  color = c("green", "red", "black"),  # Colors
  label = c("Sensitive", "Resistant", "Overall")  # Labels for annotation
)
qA = make_pe_ci_plot(df)


## panel B "Delay + No-escape + No-emerge"
## lambda differs by treatment
## ratio of resistant to sensitive differs by treatment

# pvs during the infusion period (A) 
pvsA = ((p*(1-q))/(q + p - q*p)) * (1 - exp(-lambda*(q + p - q*p)*pmin(560, 1:728)))

# pvr during the infusion period (A) which will be F_T(560)*F_S_cond_T(1:728)
pvrA = (q/(q + p - q*p)) * (1-exp(-lambda*(q + p - q*p)*560)) * F_S_cond_T(1:728, lambda*(q + p - q*p), 1/56)

# proportion at risk after week 80 in VRC01 arm
pvAtRisk = exp(-lambda*(q + p - q*p)*560)

# follow-up period
pvsB = pvAtRisk * (1-q) * (1-exp(-lambda*(pmax(0, (1:728)-560))))
pvrB = pvAtRisk * q * (1-exp(-lambda*(pmax(0, (1:728)-560))))

# overall
pvs = pvsA + pvsB
pvr = pvrA + pvrB

# placebos are simple
pps = (1-q) * (1 - exp(-lambda * (1:728)))
ppr = q * (1 - exp(-lambda * (1:728)))

# Model B plots
df = data.frame( x = 1:728, 
                 y1 = 1 - pvs / pps, 
                 y2 = 1 - pvr / ppr,
                 y3 = 1 - (pvs + pvr) / (pps + ppr))
pB = make_cartoon_plot(df, 'B. Delay + No-escape + No-emerge')

cis = pe_ci(pvs[728], pps[728], 3000, 1500)
cir = pe_ci(pvr[728], ppr[728], 3000, 1500)
cio = pe_ci(pvs[728] + pvr[728], pps[728] + ppr[728], 3000, 1500)
df = data.frame(
  x = c(1 - pvs[728] / pps[728], 1 - pvr[728] / ppr[728], 1 - (pvs[728] + pvr[728]) / (pps[728] + ppr[728])),  # Point estimates
  y = c(3, 2, 1),  # Y positions
  lower = c(cis[1], cir[1], cio[1]),  # Lower bound of CI
  upper = c(cis[2], cir[2], cio[2]),  # Upper bound of CI
  color = c("green", "red", "black"),  # Colors
  label = c("Sensitive", "Resistant", "Overall")  # Labels for annotation
)
qB = make_pe_ci_plot(df)


## panel C "Delay + No-escape + Emerge"
## lambda equal by treatment
## ratio of resistant to sensitive equal by treatment

# pvs during the infusion period (A) 
pvsA = p * (1-q) * (1 - exp(-lambda*pmin(560, 1:728)))

# pvr during the infusion period (A) which will be F_T(560)*F_S_cond_T(1:728) 
pvrA = q * (1-exp(-lambda*560)) * F_S_cond_T(1:728, lambda, 1/56)

# proportion at risk after week 80 in VRC01 arm
pvAtRisk = exp(-lambda*560)

# follow-up period
pvsB = pvAtRisk * (1-q) * (1-exp(-lambda*(pmax(0, (1:728)-560))))
pvrB = pvAtRisk * q * (1-exp(-lambda*(pmax(0, (1:728)-560))))

# sensitive viruses that escape emerge during follow-up
pvsE = c(rep(0, 560), seq(0, ((1-p)/p)*pvsA[560], length.out=168))

# sensitive viruses
# overall
pvs = pvsA + pvsB + pvsE
pvr = pvrA + pvrB 

# placebos are simple
pps = (1-q) * (1 - exp(-lambda * (1:728)))
ppr = q * (1 - exp(-lambda * (1:728)))


# Model C plots
df = data.frame( x = 1:728, 
                 y1 = 1 - pvs / pps, 
                 y2 = 1 - pvr / ppr,
                 y3 = 1 - (pvs + pvr) / (pps + ppr))
pC = make_cartoon_plot(df, 'C. Delay + No-escape + Emerge')

cis = pe_ci(pvs[728], pps[728], 3000, 1500)
cir = pe_ci(pvr[728], ppr[728], 3000, 1500)
cio = pe_ci(pvs[728] + pvr[728], pps[728] + ppr[728], 3000, 1500)
df = data.frame(
  x = c(1 - pvs[728] / pps[728], 1 - pvr[728] / ppr[728], 1 - (pvs[728] + pvr[728]) / (pps[728] + ppr[728])),  # Point estimates
  y = c(3, 2, 1),  # Y positions
  lower = c(cis[1], cir[1], cio[1]),  # Lower bound of CI
  upper = c(cis[2], cir[2], cio[2]),  # Upper bound of CI
  color = c("green", "red", "black"),  # Colors
  label = c("Sensitive", "Resistant", "Overall")  # Labels for annotation
)
qC = make_pe_ci_plot(df)


## panel D "Delay + Escape + Emerge"
## lambda equal by treatment
## ratio of resistant to sensitive equal by treatment

# pvs during the infusion period (A) 
pvsA = p * (1-q) * (1 - exp(-lambda*pmin(560, 1:728)))

# pvr during the infusion period (A) which will be F_T(560)*F_S_cond_T(1:728) 
pvrA = q * (1-exp(-lambda*560)) * F_S_cond_T(1:728, lambda, 1/56)

# proportion at risk after week 80 in VRC01 arm
pvAtRisk = exp(-lambda*560)

# follow-up period
pvsB = pvAtRisk * (1-q) * (1-exp(-lambda*(pmax(0, (1:728)-560))))
pvrB = pvAtRisk * q * (1-exp(-lambda*(pmax(0, (1:728)-560))))

# sensitive viruses that escape during infusion period and emerge during follow-up
pvrE = c(rep(0, 560), seq(0, ((1-p)/p)*pvsA[560], length.out=168))

# sensitive viruses
# overall
pvs = pvsA + pvsB
pvr = pvrA + pvrB + pvrE

# placebos are simple
pps = (1-q) * (1 - exp(-lambda * (1:728)))
ppr = q * (1 - exp(-lambda * (1:728)))

# Model D plots
df = data.frame( x = 1:728, 
                 y1 = 1 - pvs / pps, 
                 y2 = 1 - pvr / ppr,
                 y3 = 1 - (pvs + pvr) / (pps + ppr))
pD = make_cartoon_plot(df, 'D. Delay + Escape + Emerge')

cis = pe_ci(pvs[728], pps[728], 3000, 1500)
cir = pe_ci(pvr[728], ppr[728], 3000, 1500)
cio = pe_ci(pvs[728] + pvr[728], pps[728] + ppr[728], 3000, 1500)
df = data.frame(
  x = c(1 - pvs[728] / pps[728], 1 - pvr[728] / ppr[728], 1 - (pvs[728] + pvr[728]) / (pps[728] + ppr[728])),  # Point estimates
  y = c(3, 2, 1),  # Y positions
  lower = c(cis[1], cir[1], cio[1]),  # Lower bound of CI
  upper = c(cis[2], cir[2], cio[2]),  # Upper bound of CI
  color = c("green", "red", "black"),  # Colors
  label = c("Sensitive", "Resistant", "Overall")  # Labels for annotation
)
qD = make_pe_ci_plot(df)

# Arrange figure
# grid.arrange(pA+ylab(NULL), pB+ylab(NULL), pC+ylab(NULL), pD+ylab(NULL), 
#             qA+xlab(NULL), qB+xlab(NULL), qC+xlab(NULL), qD+xlab(NULL), 
#             nrow=2,
#             heights = unit(c(2,1), "null"))

pdf(file=pdfFileSave, width=8, height=4)
grid.arrange(pA+ylab(NULL), pB+ylab(NULL), pC+ylab(NULL), pD+ylab(NULL),
             nrow=1,
             heights = unit(c(1), "null"))
dev.off()



