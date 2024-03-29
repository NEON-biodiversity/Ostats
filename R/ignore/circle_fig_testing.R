# Testing: create a figure with ant data

library(Ostats)
library(ggplot2)

data(ant_data)

## edit: add a couple more columns to ant_data with fake data
ant_data$time2 <- ant_data$time %% 12
ant_data$time3 <- ant_data$time %% 6

# Discrete circular data

sp=ant_data$species
plots=ant_data$chamber
plot_dat=ant_data
traits=ant_data[, "time", drop=F]
adjust=2
scale="fixed"
n_col=2
x_limits=c(0,23)
name_x="Trait"
name_y="Number"
names(plot_dat)[2] <- "plots"
i=1
alpha <- 0.5
colorvalues <- sample(viridis::viridis(length(unique(sp)), alpha = alpha))

names(colorvalues) <- unique(sp)

ggplot2::theme_set(
  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                       axis.text = ggplot2::element_text(size = 12),
                                       axis.title = ggplot2::element_text(size = 12),
                                       axis.text.y = ggplot2::element_blank(),
                                       axis.ticks.y = ggplot2::element_blank(),
                                       strip.background = ggplot2::element_blank()))


# Discrete plot.
# Issue: If area is used, it keeps getting wider as you get further from the origin so it gives a biased depiction of the overlap.
(
  ggplot_dist <- ggplot2::ggplot(plot_dat) +
    ggplot2::geom_histogram(aes(x = time, group = sp, fill = sp), alpha = 1/3, position = 'identity', binwidth = 1) +
    #ggplot2::stat_density(adjust = adjust, ggplot2::aes_string(x = dimnames(traits)[[2]][i], group = 'sp', fill = 'sp'), alpha = alpha, geom='polygon', position = 'identity') +
    ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
    ggplot2::scale_fill_manual(values = colorvalues) +
    ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
    #ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
    ggplot2::coord_polar()
  #ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')
)

# Solution: manually bin the discrete data first, grouped by sp and plots,

## calculate manual jitter factor
jitter_width <- 10/diff(x_limits)
jitter_seq <- seq(from = -jitter_width, to = jitter_width, length.out = length(unique(sp)))

# If desired to normalize:
segment_heights <- by(plot_dat, list(sp, plots), function(x) cbind(sp = x$sp[1], plots = x$plots[1], as.data.frame.table(table(x$time)/sum(x$time), stringsAsFactors = FALSE)))

# If not desired to normalize:
segment_heights <- by(plot_dat, list(sp, plots), function(x) cbind(sp = x$sp[1], plots = x$plots[1], as.data.frame.table(table(x$time), stringsAsFactors = FALSE)))

plot_binned <- do.call(rbind, segment_heights)
plot_binned$Var1 <- as.numeric(plot_binned$Var1)

# Jitter manually
plot_binned$Var1 <- plot_binned$Var1 + jitter_seq[plot_binned$sp]

# This does a great job of resolving the issue except for the overplotting
# But if position_dodge is used it will cause them to curve
# Use manual jittering to deal with this

(
  ggplot_dist <- ggplot2::ggplot(plot_binned) +
    ggplot2::geom_segment(aes(x = Var1, xend = Var1, y = 0, yend = Freq, group = sp, color = sp), alpha = 1/2, size = 1.2) +
    ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
    ggplot2::scale_color_manual(values = colorvalues) +
    ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
    ggplot2::scale_y_continuous(name = name_y) +
    ggplot2::coord_polar() +
    ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')
)


### Use circular::density() to get kernel density estimate for continuous circular data and plot.
## test for one combination

library(circular)
x <- plot_dat$time[plot_dat$species == "Camponotus castaneus" & plot_dat$plots == 1]
xcirc <- circular(x, units = 'hours')

xcircdens <- density(xcirc, bw = 24)

xcircdens_data <- with(xcircdens, data.frame(x=x, y=y))

## use by() to get for all plots and species
calc_circ_dens <- function(dat) {
  xcirc <- circular(dat$time, units = 'hours')
  xcircdens <- density(xcirc, bw = max(x_limits))
  cbind(sp = dat$sp[1], plots = dat$plots[1], with(xcircdens, data.frame(x=x, y=y)))
}

circ_dens_data <- by(plot_dat, list(sp, plots), calc_circ_dens)
plot_binned <- do.call(rbind, circ_dens_data)

(
ggplot_dist <- ggplot2::ggplot(plot_binned, aes(x=x, y=y, fill=sp)) +
  ggplot2::geom_polygon(alpha = 1/3) +
    ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
  ggplot2::scale_fill_manual(values = colorvalues) +
  ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
  ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
  ggplot2::coord_polar()
)

### Means plotted as radii for continuous circular data ggplot_means.
calc_circ_mean <- function(x) {
  xcirc <- circular(x, units = 'hours')
  mean(xcirc, na.rm = TRUE)
}

taxon_mean <- stats::aggregate(traits, list(sp, plots), calc_circ_mean)
names(taxon_mean) <- c('sp', 'plots', dimnames(traits)[[2]])
taxon_mean <- taxon_mean[taxon_mean$plots %in% use_plots, ]

(
ggplot_means <- ggplot2::ggplot(taxon_mean) +
  ggplot2::geom_vline(ggplot2::aes_string(xintercept = dimnames(traits)[[2]][i], colour = 'sp', group='sp'), alpha = alpha, size=0.5, key_glyph = 'rect') +
  ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
  ggplot2::scale_colour_manual(values = colorvalues) +
  ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
  ggplot2::scale_y_continuous(expand = c(0,0)) +
  ggplot2::coord_polar() +
  ggplot2::theme(legend.position = if (!legend) 'none' else 'right')
)

### Discrete non-circular data (example)
bin_width = 1
legend = TRUE
means = FALSE

ggplot_dist <- ggplot2::ggplot(plot_dat) +
  ggplot2::geom_histogram(ggplot2::aes_string(x = dimnames(traits)[[2]][i], y = 'ggplot2::after_stat(density * width)', fill = 'sp'), alpha = alpha, position = 'identity', binwidth = bin_width)

ggplot_dist <- ggplot_dist +
  ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
  ggplot2::scale_fill_manual(values = colorvalues) +
  ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
  ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
  ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')
