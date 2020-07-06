library(tidyverse)

overlap_dat <- read_csv('internal_docs/overlap_dat.csv')
indiv_dat <- read_csv('internal_docs/indiv_dat.csv')


sites2use <- c('BART','KONZ','JORN')
site_name_labeller <- labeller(siteID = c(BART='Bartlett', KONZ='Konza', JORN='Jornada'))

#
colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = 24, replace = TRUE)

# Set the theme. Maybe don't include this in the final plotting function, but this is what we used for the MS.
theme_set(
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 18),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     legend.position = 'none',
                     strip.background = element_blank())
)

ggplot(filter(indiv_dat, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites2use))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
  facet_wrap(~ siteID, ncol = 1, labeller = site_name_labeller) +
  scale_fill_manual(values = colorvalues) +
  scale_x_continuous(name = 'Body Mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3)) +
  scale_y_continuous(name = 'Probability Density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = overlap_dat %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites2use)))
