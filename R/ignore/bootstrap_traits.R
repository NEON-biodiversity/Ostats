# Illustration of how to bootstrap O-statistics
# Using the example data from the vignette (small mammal body weights from HARV and JORN only)

library(Ostats)
library(dplyr)
library(tidyr)
library(ggplot2)

### Code copied from vignette to load a small example dataset and calculate O-stats on the observed dataset.
dat <- small_mammal_data[small_mammal_data$siteID %in% c('HARV', 'JORN'),
                         c('siteID', 'taxonID', 'weight')]
dat <- dat[!is.na(dat$weight), ]
dat$log_weight <- log10(dat$weight)

Ostats_example <- Ostats(traits = as.matrix(dat[,'log_weight', drop = FALSE]),
                         sp = factor(dat$taxonID),
                         plots = factor(dat$siteID),
                         run_null_model = FALSE)

# Convert observed normalized overlap to a data frame for plotting
Ostats_observed <- with(Ostats_example, data.frame(siteID = dimnames(overlaps_norm)[[1]], overlap_norm = overlaps_norm[,1]))

### Make a custom function that can be run for each bootstrap subsample.
### Following the suggestion of Maitner et al., we sample with replacement from each species' traits within each community
### **proportional to** that species' relative abundance in that community.
### We use the nonparametric bootstrap approach. In other words we are sampling directly from the raw data, instead of fitting a distribution
### and sampling from that (that would be parametric because we would be using the parameters of a fitted distribution to generate bootstrap samples).

boot_fn <- function(i) {
  # Bootstrap sample of same size as the data, within each community each row is sampled with replacement with equal probability
  # This means the proportion of individuals of each species in the bootstrap sample will be similar to the observed data
  # Currently uses dplyr. Can later translate this to base R code.
  bsample <- dat %>%
    group_by(siteID) %>%
    slice_sample(prop = 1, replace = TRUE)

  Ostats_bsample <- Ostats(traits = as.matrix(bsample[,'log_weight', drop = FALSE]),
                           sp = factor(bsample$taxonID),
                           plots = factor(bsample$siteID),
                           run_null_model = FALSE)

  # Extract only normalized overlap
  with(Ostats_bsample, data.frame(rep = i, data.frame(siteID = dimnames(overlaps_norm)[[1]], overlap_norm = overlaps_norm[,1])))
}

n_boot <- 100

Ostats_boot <- lapply(1:n_boot, boot_fn)
Ostats_boot <- do.call(rbind, Ostats_boot)

# Extract the results from the list and turn into a data frame. Plot and display summary statistics.
# In the plot, the solid line is the observed O-stat and the histogram is the bootstrap distribution of same.
ggplot(Ostats_boot, aes(x = overlap_norm, group = siteID, color = siteID, fill = siteID)) +
  geom_histogram(position = 'identity', alpha = 0.5, color = NA) +
  geom_vline(aes(xintercept = overlap_norm, color = siteID), data = Ostats_observed, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_color_manual(values = unname(palette.colors(3))[-1], aesthetics = c('colour', 'fill')) + # Colorblind accessible
  theme_bw()

# Calculate median 95% quantile interval of the bootstrap distribution of O-stats. Display alongside the observed.
Ostats_boot_summary_stats <- Ostats_boot %>%
  group_by(siteID) %>%
  reframe(data.frame(stat = c('median', '95% lower', '95% upper'), value = quantile(overlap_norm, probs = c(0.5, 0.025, 0.975)))) %>%
  pivot_wider(names_from = stat, values_from = value)

left_join(Ostats_observed, Ostats_boot_summary_stats) %>%
  setNames(c('siteID', 'Observed Ostat', 'Bootstrap median', 'Bootstrap 95% lower', 'Bootstrap 95% upper'))
