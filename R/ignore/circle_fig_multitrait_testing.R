# Testing: create a figure with ant data

library(Ostats)

data(ant_data)

## edit: add a couple more columns to ant_data with fake data
ant_data$time2 <- ant_data$time %% 12
ant_data$time3 <- ant_data$time %% 6

## Test all combinations

sp=ant_data$species
plots=ant_data$chamber
traits=ant_data[, c('time', 'time2', 'time3')]
color_values <- c("#E41A1C", "#377EB8", "#4DAF4A")

# Not circular
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = TRUE, means = TRUE) # OK
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = FALSE, means = TRUE) # OK
suppressWarnings(Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = TRUE)) # OK
Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = FALSE) # OK

# Circular
pc1 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = TRUE, circular = TRUE)
pc2 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = FALSE, circular = TRUE)
pc3 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = TRUE, circular = TRUE)
pc4 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = TRUE, normalize = FALSE, circular = TRUE)
pc5 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = TRUE, circular = TRUE, means = TRUE)
pc6 <- Ostats_plot(plots = plots, sp = sp, traits = traits, colorvalues = color_values, discrete = FALSE, normalize = FALSE, circular = TRUE, means = TRUE)
