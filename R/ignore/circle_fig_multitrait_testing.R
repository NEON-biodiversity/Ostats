# Testing: create a figure with ant data

library(Ostats)
library(ggplot2)

data(ant_data)

## edit: add a couple more columns to ant_data with fake data
ant_data$time2 <- ant_data$time %% 12
ant_data$time3 <- ant_data$time %% 6

## Test all combinations

sp=ant_data$species
plots=ant_data$chamber
traits=ant_data[, c('time', 'time2', 'time3')]
color_values <- RColorBrewer::brewer.pal(3, 'Set1')

# Not circular
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = TRUE, means = TRUE)
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = FALSE, means = TRUE)
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = TRUE, limits_x = c(1, 1))
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = FALSE, limits_x = c(1, 1))

# Circular
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = TRUE, circular = TRUE)
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = FALSE, circular = TRUE)
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = TRUE, circular = TRUE, limits_x = c(1, 1))
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = FALSE, circular = TRUE, limits_x = c(1, 1))
