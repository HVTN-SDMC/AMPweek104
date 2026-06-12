library(tidyverse)

# generate a forest plot for the PE at week 104
df = read.csv('../data/amp_cir_efficacy_postwk80_trunc.csv')

# output (note, TABLE data was added to plot in illustrator)
PLOT = '../output/figures/AMP_week80to104_cir_forest_plot_trunc.pdf'
TABLE = '../output/tables/AMP_week80to104_cir_forest_plot_trunc.csv'

# Match the order of set and label to Corey et al. Fig. 1C
df$label = factor(df$comparison,
                   levels = c("Pooled VRC01 vs. Control", "VRC01 10 mg/kg vs. Control", "VRC01 30 mg/kg vs. Control"),
                   labels = c('VRC01 pooled', 'VRC01, low dose', 'VRC01, high dose'))

df = df[order(df$label),]
df$ypos = nrow(df):1
df$pe = as.numeric(sub('%', '', df$pe, fixed=TRUE))
df$lci = unlist(lapply(df$bounds, function(v) {
  v = strsplit(v, ' to ')[[1]]
  v = v[1]
  v = sub('(', '', v, fixed=TRUE)
  v = sub('%', '', v, fixed=TRUE)
  as.numeric(v)
}))
df$uci = unlist(lapply(df$bounds, function(v) {
  v = strsplit(v, ' to ')[[1]]
  v = v[2]
  v = sub(')', '', v, fixed=TRUE)
  v = sub('%', '', v, fixed=TRUE)
  as.numeric(v)
}))


pdf(PLOT, width=8, height=6)
par(mar=c(4,12,1,1))
plot(0, type='n', xlim=c(-100,100), ylim=c(0,nrow(df)+1), axes=FALSE, xlab='', ylab='')
abline(v=0, lty=2)
axis(1, at=seq(-100, 100, by=20))
axis(2, at=df$ypos, labels=df$label, las=2)

for( i in 1:nrow(df) ) {
  if( is.na(df$pe[i]) ) next

  arrows(x0 = max(-98.0, df$lci[i]), y0 = df$ypos[i],
         x1 = df$uci[i], y1 = df$ypos[i],
         col = '#3A687E',
         code = ifelse(df$lci[i] < -98, 1, 0),      # 1 = arrow at the start (x0, y0), 2 = at end, 3 = both
         length = 0.1,  # Adjusts the size of the arrowhead
         lwd = 1)       # Adjusts line width
  points(df$pe[i], df$ypos[i], pch=15, col='#3A687E')
}
dev.off()

out = df %>%
  select(label, pe, lci, uci, pvalue) %>%
  mutate(PE = round(pe, 1),
         bounds = sprintf("(%0.0f, %0.1f)", lci, uci),
         pvalue = sprintf("%0.3f", pvalue)) %>%
  select(label, PE, bounds, pvalue)

write.csv(out, TABLE, row.names = FALSE)

q(save='no')





