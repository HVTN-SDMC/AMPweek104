library(tidyverse)
library(patchwork)
library(grid)
library(gridExtra)
library(cowplot)

library(here)
here::i_am("README.md")
repoDir <- here::here()
macroDir <- file.path(repoDir, "code/macro")
datDir <- file.path(repoDir, "data")
dat2Dir <- "/Volumes/trials/vaccine/p704/analysis/public_use_data/postwk80/public_use_data_final" # file.path(repoDir, "data")
figDir <- file.path(repoDir, "output/figures")
tabDir <- file.path(repoDir, "output/tables")

# Data files
dataFile_postwk80 <- file.path(dat2Dir, "amp_survival_postwk80.csv")
dataFile <- file.path(dat2Dir, "amp_survival.csv")
data703cov <- file.path(dat2Dir, 'v703_subject_master.csv')
data704cov <- file.path(dat2Dir, 'v704_subject_master.csv')
dataWT <- file.path(dat2Dir, "amp_baseline_bodywt.csv")

pdfFileSave = file.path(figDir, 'HIV_prognostic_factors_categorized_by_postwk80_cohort.pdf')

# Post week 80 cohort
dat_postwk80 = read.csv(dataFile_postwk80)

# MITT cohort
dat = read.csv(dataFile)
dat = dat %>% filter(efficacy_flag==1)

# add cohort indicator to main data
df = dat %>%
  mutate(postwk80_cohort = ifelse(pub_id %in% dat_postwk80$pub_id, 1, 0))

# get covariate info
cov.703 = read.csv(data703cov)
cov.704 = read.csv(data704cov)
n = intersect(names(cov.703), names(cov.704))
cov = rbind(cov.703[,n], cov.704[,n])
cov = cov[,c("pub_id","country","age")]
cov$RSA = as.numeric(cov$country=='South Africa')
cov$SA  = as.numeric(cov$country %in% c('Brazil','Peru'))

# get weight
wt = read.csv(dataWT) %>%
  select(pub_id, visit, RIBwtkg)

df = df %>%
  merge(cov) %>%
  merge(wt) %>%
  mutate(postwk80_lab = factor(postwk80_cohort, levels=c(1,0), labels=c('Eligible', 'Non-eligible')),
         trial = ifelse(protocol=='HVTN 703', 'Africa', 'Americas'),
         trt = factor(rx_code, levels=c('C3', 'T1','T2'), labels=c('Placebo', 'Low-dose', 'High-dose')))


# Build a unified region label column
df_region <- df |>
  mutate(
    region = case_when(
      trial == "Africa"   & RSA == 1 ~ "RSA",
      trial == "Africa"   & RSA == 0 ~ "Non-RSA",
      trial == "Americas"  & SA  == 1 ~ "SA",
      trial == "Americas"  & SA  == 0 ~ "Non-SA"
    )
  )

# Summarise proportions by cohort × trial × treatment × region
df_prop <- df_region |>
  count(postwk80_lab, trial, trt, region) |>
  group_by(postwk80_lab, trial, trt) |>
  mutate(
    prop  = n / sum(n),
    denom = sum(n)
  ) |>
  ungroup()

# Max relative difference between eligible and non-eligible by treatment × region
cohort_diff <- df_prop |>
  select(postwk80_lab, trial, trt, region, prop) |>
  pivot_wider(names_from = postwk80_lab, values_from = prop) |>
  mutate(
    rel_diff = abs(Eligible - `Non-eligible`) / `Non-eligible`
  )

max_cohort <- cohort_diff |> slice_max(rel_diff, n = 1)
cat("Max relative difference (eligible vs non-eligible):\n")
print(max_cohort)

# Max percentage point difference between placebo and low- or high-dose
# among eligible only
trt_diff <- df_prop |>
  filter(postwk80_lab == "Eligible") |>
  select(trial, trt, region, prop) |>
  pivot_wider(names_from = trt, values_from = prop) |>
  mutate(
    diff_low  = abs(Placebo - `Low-dose`)  / Placebo,
    diff_high = abs(Placebo - `High-dose`) / Placebo
  ) |>
  pivot_longer(cols = c(diff_low, diff_high),
               names_to = "comparison", values_to = "rel_diff")

max_trt <- trt_diff |> slice_max(rel_diff, n = 1)

cat("\nMax relative difference (placebo vs low- or high-dose, eligible):\n")
print(max_trt)

# Shared aesthetics as a function to avoid repetition
region_plot <- function(data, trial_name) {
  ggplot(data, aes(x = postwk80_lab, y = prop, fill = postwk80_lab)) +
    geom_col(
      position = position_dodge(width = 0.7),
      width = 0.6, alpha = 0.85
    ) +
    geom_text(
      aes(label = paste0(n, "/", denom, "\n", scales::percent(prop, accuracy = 1)),
          group = postwk80_lab),
      position = position_dodge(width = 0.7),
      vjust = -0.4, size = 2.8, color = "grey30"
    ) +
    facet_grid(. ~ trt + region) +
    scale_fill_brewer(palette = "Set1", name = "Weeks 80-104", direction = -1) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x       = "",
      y       = "Proportion",
      title   = trial_name    ) +
    theme_bw() +
    theme(
      strip.background  = element_rect(fill = "grey92"),
      strip.text        = element_text(face = "bold"),
      axis.text.x       = element_text(angle = 30, hjust = 1),
      legend.position   = "top",
      plot.caption      = element_text(hjust = 0, size = 8, color = "grey40"),
      plot.caption.position = "plot"
    )
}

p_africa   <- region_plot(filter(df_prop, trial == "Africa"),  "Africa trial")
p_americas <- region_plot(filter(df_prop, trial == "Americas"), "Americas trial")


# Region proportions graph
# Bars show the proportion of participants in each sub-region by cohort and treatment arm.
# Labels show n/N (%) where N is the total within each treatment arm and region.
# RSA = South Africa; Non-RSA = other African sites.
  
p_region = p_africa / p_americas +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# Plot Weight
# Weight (kg) distribution by cohort (Weeks 80–104 eligibility status), 
# trial (Africa, Americas), and treatment arm (placebo, low dose, high dose). 
# Violins show the full distribution with kernel density estimation trimmed to 
# observed data range; internal box plots show median (center line), interquartile 
# range (box), and 1.5×IQR (whiskers). Outliers are suppressed.

p_wt = ggplot(df, aes(x = postwk80_lab, y = RIBwtkg, fill = postwk80_lab)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.9) +
  facet_grid(trial ~ trt) +
  scale_fill_brewer(palette = "Set1", name = "Weeks 80-104", direction = -1) +
  scale_y_continuous(limits = range(df$RIBwtkg), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "",
    y = "Weight (kg)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

# Plot Age
# Age (years) distribution by cohort (Weeks 80–104 eligibility status), 
# trial (Africa, Americas), and treatment arm (placebo, low dose, high dose). 
# Violins show the full distribution with kernel density estimation trimmed to 
# observed data range; internal box plots show median (center line), interquartile 
# range (box), and 1.5×IQR (whiskers). Outliers are suppressed.

p_age = ggplot(df, aes(x = postwk80_lab, y = age, fill = postwk80_lab)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.9) +
  facet_grid(trial ~ trt) +
  scale_fill_brewer(palette = "Set1", name = "Weeks 80-104", direction = -1) +
  scale_y_continuous(limits = range(df$age), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "",
    y = "Age (years)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

# Remove individual legends and titles from sub-plots
p_age_clean <- p_age +
  theme(legend.position = "none")

p_wt_clean <- p_wt +
  theme(legend.position = "none")

p_africa_clean <- p_africa +
  theme(legend.position = "none")

p_americas_clean <- p_americas +
  theme(legend.position = "none")


# Extract legend from any of the plots
shared_legend <- get_legend(
  p_wt + theme(legend.position = "bottom")
)

# Clean versions with no legends
p_region_clean <- p_africa_clean / p_americas_clean


# patchworkGrob for patchwork objects, ggplotGrob for plain ggplots
trim_margins <- theme(
  plot.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6)
)

grob_region <- patchworkGrob(p_region_clean & trim_margins)
grob_wt     <- ggplotGrob(p_wt_clean  + trim_margins)
grob_age    <- ggplotGrob(p_age_clean + trim_margins)


# Add tags manually
add_tag <- function(grob, label) {
  arrangeGrob(grob, top = textGrob(label, x = 0.02, hjust = 0,
                                   gp = gpar(fontface = "bold", fontsize = 14),
                                   vjust = 1.5))  # pull down toward panel top
}

grob_region <- add_tag(patchworkGrob(p_region_clean & trim_margins), "A")
grob_wt     <- add_tag(ggplotGrob(p_wt_clean  + trim_margins), "B")
grob_age    <- add_tag(ggplotGrob(p_age_clean + trim_margins), "C")

# Stack with explicit heights and legend at bottom
final <- plot_grid(
  grob_region,
  grob_wt,
  grob_age,
  shared_legend,
  ncol        = 1,
  rel_heights = c(2, 1, 1, 0.2)
)

pdf(pdfFileSave, width = 8.5, heigh=11)
print(final)
dev.off()

q(save='no')