#'@title  Plot community trait overlap
#'
#'@description This function plots the overlap of traits among
#'  species for each community. If there are multiple traits,
#'  each trait is plotted separately in one-dimensional space.
#'
#'@param plots Site identity: a vector of names of each community.
#'@param sp Taxon identity: a vector of species or taxa names.
#'@param traits A vector of trait measurements for each individual, or a matrix
#' or data frame with rows representing individuals and columns representing traits.
#'@param overlap_dat Optional: an object containing the output of \code{\link{Ostats}}.
#' If provided, overlap statistics will be displayed in the plot panels.
#'@param use_plots a vector of sites to plot. If NULL, the function will plot all the sites.
#'@param n_col Number of columns for layout of individual panels. Default is 1.
#'@param colorvalues Vector of color values for the density polygons.
#' Defaults to a viridis palette if none provided.
#'@param alpha defines the transparency level for the density polygons. Default is 0.5.
#'@param adjust the bandwidth adjustment of the density polygons. Default is 2.
#' See \code{\link[stats]{density}}. Only used if \code{discrete = FALSE}.
#'@param bin_width the width of each bin of the histograms. Default is 1.
#' Only used if \code{discrete = TRUE}.
#'@param limits_x Vector of length 2, with multiplicative factor to apply to the minimum
#' and maximum values of each trait to expand the limits of the x axis.
#' Default is \code{c(0.5, 1.5)}, or 0.5 times the minimum and 1.5 times the maximum
#' value of each trait, for continuous traits. For discrete traits the default is
#' \code{c(1, 1)} or no expansion of limits.
#'@param legend Whether to include a legend. Defaults to \code{FALSE}.
#'@param scale If you want the scale of x, y or both x and y axis to be independent,
#' set the argument to "free_x", "free_y" or "free" respectively.
#' Default = "fixed" which uses the same scale across all sites.
#' See \code{\link[ggplot2]{facet_grid}}.
#'@param name_x x-axis label. Default is 'trait value'
#'@param name_y y-axis label. Default is 'probability density'
#'@param normalize if \code{TRUE}, areas of density plots are normalized to be equal
#' across taxa; if \code{FALSE}, areas will be proportional to abundance.
#' Default is \code{TRUE}.
#'@param means if \code{TRUE}, trait means for each species are plotted in an additional plot
#' column next to the traits distribution plots for each site. Default is \code{FALSE}.
#'@param discrete if \code{TRUE}, plots histograms at discrete trait values instead
#' of smooth kernel density plots. Default is \code{FALSE}.
#'@param circular if \code{TRUE}, plots density plots or histograms using polar
#' coordinates, and estimates density using method for objects of class
#' \code{circular}. Default is \code{FALSE}.
#' @param circular_args optional list of additional arguments to pass to
#'  \code{\link[circular]{circular}}. Only used if \code{circular = TRUE} and
#'  \code{discrete = FALSE}. If no arguments are provided, default arguments to
#'  \code{\link[circular]{circular}} are used.
#'
#'@return Density plots of species trait distributions plotted together
#'  for each community to show how they overlap each other. Each community
#'  is plotted on a separate panel within a multipanel figure.
#'  The overlap value obtained as output from \code{\link{Ostats}}
#'  is labelled on each community graph, if provided by the user.
#'
#'  If trait values are discrete rather than continuous, histograms are
#'  plotted instead of kernel density plots.
#'
#'  If trait values are circular, a circular kernel density estimate for
#'  each species is plotted on a polar coordinate plot. If trait values are
#'  both circular and discrete, a "sunburst" plot is returned.
#'
#'  The class of the returned object is \code{Ostats_plot_object}. Calling
#'  \code{print} on this object will draw the plot using
#'  \code{\link[grid]{grid.draw}}.
#'
#'  If more than one trait is provided, a list of objects of class
#'  \code{Ostats_plot_object} will be returned.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'# set the arguments:
#'plots <- small_mammal_data$siteID
#'sp <- small_mammal_data$taxonID
#'traits <- log10(small_mammal_data$weight)
#'
#'# to plot only selected sites:
#'use_plots <- c('BART','KONZ','JORN')
#'
#'
#'Ostats_plot(plots = plots, sp = sp, traits = traits,
#'            overlap_dat = small_mammal_Ostats,
#'            use_plots = use_plots, means = TRUE)


#'@export
#'
Ostats_plot<-function(plots,
                      sp,
                      traits,
                      overlap_dat = NULL,
                      use_plots = NULL,
                      n_col = 1,
                      scale = "fixed",
                      colorvalues = NULL,
                      alpha = 0.5,
                      adjust = 2,
                      bin_width = 1,
                      limits_x = NULL,
                      legend = FALSE,
                      name_x = 'trait value',
                      name_y = 'probability density',
                      normalize = TRUE,
                      means = FALSE,
                      circular = FALSE,
                      discrete = FALSE,
                      circular_args = list()) {

  if (means & discrete) stop('Plotting trait means is not supported for discrete traits.')

  if (missing(limits_x) & !discrete) limits_x <- c(0.5, 1.5)
  if (missing(limits_x) & discrete) limits_x <- c(1, 1)

  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(use_plots)) {
    use_plots <- unique(plots)
  }

  # Check dimensions of traits. Create one-column matrix if it is a vector.
  if (is.vector(traits)) {
    traits <- as.matrix(traits)
  }

  # Give traits matrix names if it does not have any.
  if (is.null(dimnames(traits)[[2]])) {
    dimnames(traits)[[2]] <- paste('trait', 1:ncol(traits), sep = '_')
  }

  # Filter overlap statistics only for use_plots
  if (!is.null(overlap_dat)) {
    ostat_norm <- overlap_dat$overlaps_norm
    ostat_norm <- ostat_norm[rownames(ostat_norm) %in% use_plots, , drop = FALSE]
  }

  plot_dat <- cbind(as.data.frame(traits), sp = sp, plots = plots)
  plot_dat <- plot_dat[plots %in% use_plots, ]


  # Calculate mean value by taxon.
  if (circular) {
    mean_fn <- function(x) {
      xcirc <- do.call(circular::circular, c(list(x = x), circular_args))
      mean(xcirc, na.rm = TRUE)
    }
  } else {
    mean_fn <- function(x) mean(x, na.rm = TRUE)
  }
  taxon_mean <- stats::aggregate(traits, list(sp, plots), mean_fn)
  names(taxon_mean) <- c('sp', 'plots', dimnames(traits)[[2]])
  taxon_mean <- taxon_mean[taxon_mean$plots %in% use_plots, ]

  # If a color vector is not provided, create a default palette.
  if (is.null(colorvalues)) {
    colorvalues <- sample(viridis::viridis(length(unique(sp)), alpha = alpha))
  }

  names(colorvalues) <- unique(sp)

  plot_list <- list()

  ggplot2::theme_set(
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                         axis.text = ggplot2::element_text(size = 12),
                                         axis.title = ggplot2::element_text(size = 12),
                                         axis.text.y = ggplot2::element_blank(),
                                         axis.ticks.y = ggplot2::element_blank(),
                                         strip.background = ggplot2::element_blank()))

  for (i in 1:ncol(traits)) {

    if (!is.null(overlap_dat)) {
      overlap_labels <- data.frame(plots = row.names(ostat_norm),
                                   lab = paste('Overlap =', signif(ostat_norm[,i], 2)))
    }

    x_limits <- limits_x * range(traits[, i], na.rm = TRUE)

    if (!circular) {

      if (!discrete) {
        if (normalize) {
          ggplot_dist <- ggplot2::ggplot(plot_dat) +
            ggplot2::geom_density(adjust = adjust, ggplot2::aes(x = .data[[dimnames(traits)[[2]][i]]], group = sp, fill = sp), alpha = alpha, position = 'identity', color = NA)
        } else {
          ggplot_dist <- ggplot2::ggplot(plot_dat) +
            ggplot2::geom_density(adjust = adjust, ggplot2::aes(x = .data[[dimnames(traits)[[2]][i]]], y = ggplot2::after_stat(count), group = sp, fill = sp), alpha = alpha, position = 'identity', color = NA)
        }
      } else {
        if (normalize) {
          ggplot_dist <- ggplot2::ggplot(plot_dat) +
            ggplot2::geom_histogram(ggplot2::aes(x = .data[[dimnames(traits)[[2]][i]]], y = ggplot2::after_stat(density * width), fill = sp), alpha = alpha, position = 'identity', binwidth = bin_width)
        } else {
          ggplot_dist <- ggplot2::ggplot(plot_dat) +
            ggplot2::geom_histogram(ggplot2::aes(x = .data[[dimnames(traits)[[2]][i]]], fill = sp), alpha = alpha, position = 'identity', binwidth = bin_width)
        }
      }

      ggplot_dist <- ggplot_dist +
        ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
        ggplot2::scale_fill_manual(values = colorvalues) +
        ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
        ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
        ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')

      if (!is.null(overlap_dat)) {
        ggplot_dist <- ggplot_dist +
          ggplot2::geom_text(ggplot2::aes(label = lab), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1)
      }


      if (means) {
        ggplot_means <- ggplot2::ggplot(taxon_mean) +
          ggplot2::geom_vline(ggplot2::aes(xintercept = .data[[dimnames(traits)[[2]][i]]], colour = sp, group = sp), alpha = alpha, linewidth = 0.5, key_glyph = 'rect') +
          ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
          ggplot2::scale_colour_manual(values = colorvalues) +
          ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
          ggplot2::scale_y_continuous(expand = c(0,0)) +
          ggplot2::theme(legend.position = if (!legend) 'none' else 'right')

        ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ggplot_means, ncol = 2, widths = if (!legend) c(1, 1) else c(1, 1.3))
      } else {
        ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ncol = 1)
      }

    } else {
      if (!discrete) {

        calc_circ_dens <- function(dat) {
          xcirc <- do.call(circular::circular, c(list(x = dat[,i]), circular_args))
          xcircdens <- circular::density.circular(xcirc, bw = diff(x_limits))
          out <- cbind(sp = dat$sp[1], plots = dat$plots[1], with(xcircdens, data.frame(x=x, y=y)))
          if (!normalize) out$y <- out$y * nrow(dat)
          return(out)
        }

        circ_dens_data <- by(plot_dat, list(sp, plots), calc_circ_dens)
        plot_binned <- do.call(rbind, circ_dens_data)

        ggplot_dist <- ggplot2::ggplot(plot_binned, ggplot2::aes(x=x, y=y, fill=sp)) +
          ggplot2::geom_polygon(alpha = alpha) +
          ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
          ggplot2::scale_fill_manual(values = colorvalues) +
          ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
          ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
          ggplot2::coord_polar() +
          ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')

        if (means) {

          ggplot_means <- ggplot2::ggplot(taxon_mean) +
            ggplot2::geom_vline(ggplot2::aes(xintercept = .data[[dimnames(traits)[[2]][i]]], colour = sp, group = sp), alpha = alpha, linewidth = 0.5, key_glyph = 'rect') +
            ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
            ggplot2::scale_colour_manual(values = colorvalues) +
            ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
            ggplot2::scale_y_continuous(expand = c(0,0)) +
            ggplot2::coord_polar() +
            ggplot2::theme(legend.position = if (!legend) 'none' else 'right')

          ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ggplot_means, ncol = 2, widths = if (!legend) c(1, 1) else c(1, 1.3))
        } else {
          ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ncol = 1)
        }


      } else {

        # If data are discrete and circular, generate bin counts manually (not necessary for discrete non-circular)

        # calculate manual jitter factor (2% of data width)
        jitter_width <- diff(x_limits) * 0.02
        jitter_seq <- seq(from = -jitter_width, to = jitter_width, length.out = length(unique(sp)))

        if (normalize) {
          segment_heights <- by(plot_dat, list(sp, plots), function(x) cbind(sp = x$sp[1], plots = x$plots[1], as.data.frame.table(table(x[,i])/sum(x[,i]), stringsAsFactors = FALSE)))
        } else {
          segment_heights <- by(plot_dat, list(sp, plots), function(x) cbind(sp = x$sp[1], plots = x$plots[1], as.data.frame.table(table(x[,i]), stringsAsFactors = FALSE)))
        }

        plot_binned <- do.call(rbind, segment_heights)
        plot_binned$Var1 <- as.numeric(plot_binned$Var1)

        # Jitter manually
        plot_binned$Var1 <- plot_binned$Var1 + jitter_seq[plot_binned$sp]

        ggplot_dist <- ggplot2::ggplot(plot_binned) +
          ggplot2::geom_segment(ggplot2::aes(x = Var1, xend = Var1, y = 0, yend = Freq, group = sp, color = sp), alpha = 1/2, linewidth = 1.2) +
          ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
          ggplot2::scale_color_manual(values = colorvalues) +
          ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
          ggplot2::scale_y_continuous(name = name_y) +
          ggplot2::coord_polar() +
          ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')

        ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ncol = 1)
      }
    }

    # Assign class attribute Ostats_plot_object so that the plot has a default print method.
    attr(ggplot_dist, 'class') <- c('Ostats_plot_object', attr(ggplot_dist, 'class'))

    plot_list[[i]] <- ggplot_dist
  }

  names(plot_list) <- dimnames(traits)[[2]]

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }

}

#' @method print Ostats_plot_object
#' @export
print.Ostats_plot_object <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x, ...)
}
