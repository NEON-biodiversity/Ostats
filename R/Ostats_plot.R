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
#' See \code{\link[stats]{density}}.
#'@param limits_x Vector of length 2, with multiplicative factor to apply to the minimum
#' and maximum values of each trait to expand the limits of the x axis.
#' Default is 0.5 times the minimum and 1.5 times the maximum value of each trait.
#'@param legend Whether to include a legend. Defaults to \code{FALSE}.
#'@param scale If you want the scale of x, y or both x and y axis to be independent,
#' set the argument to "free_x", "free_y" or "free" respectively.
#' Default = "fixed" which uses the same scale across all sites.
#' See \code{\link[ggplot2]{facet_grid}}.
#'@param name_x x-axis label. Default is 'trait value'
#'@param name_y y-axis label. Default is 'probability density'
#'@param means if TRUE, trait means for each species are plotted in an additional plot
#' column next to the traits distribution plots for each site. Default is
#'
#'@return Density plots of species trait distributions plotted together
#'  for each community to show how they overlap each other.
#'  The overlap value obtained as output from \code{\link{Ostats}}
#'  is labelled on each community graph, if provided by the user.
#'
#'  The class of the returned object is \code{Ostats_plot_object}. Calling
#'  \code{print} on this object will invoke a method to draw the plot using
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
                      limits_x = c(0.5, 1.5),
                      legend = FALSE,
                      name_x = 'trait value',
                      name_y = 'probability density',
                      means = FALSE) {



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
  taxon_mean <- stats::aggregate(traits, list(sp, plots), mean, na.rm = TRUE)
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

    ggplot_dist <- ggplot2::ggplot(plot_dat) +
      ggplot2::stat_density(adjust = adjust, ggplot2::aes_string(x = dimnames(traits)[[2]][i], group = 'sp', fill = 'sp'), alpha = alpha, geom='polygon', position = 'identity') +
      ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
      ggplot2::scale_fill_manual(values = colorvalues) +
      ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
      ggplot2::scale_y_continuous(name = name_y, expand = c(0,0)) +
      ggplot2::theme(legend.position = if (!legend | means) 'none' else 'right')

    if (!is.null(overlap_dat)) {
      ggplot_dist <- ggplot_dist +
        ggplot2::geom_text(ggplot2::aes_string(label = 'lab'), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1)
    }


    if (means) {
      ggplot_means <- ggplot2::ggplot(taxon_mean) +
        ggplot2::geom_vline(ggplot2::aes_string(xintercept = dimnames(traits)[[2]][i], colour = 'sp', group='sp'), alpha = alpha, size=0.5, key_glyph = 'rect') +
        ggplot2::facet_wrap(~ plots, ncol = n_col, scales = scale) +
        ggplot2::scale_colour_manual(values = colorvalues) +
        ggplot2::scale_x_continuous(name = name_x, limits = x_limits) +
        ggplot2::scale_y_continuous(expand = c(0,0)) +
        ggplot2::theme(legend.position = if (!legend) 'none' else 'right')

      ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ggplot_means, ncol = 2, widths = if (!legend) c(1, 1) else c(1, 1.3))
    } else {
      ggplot_dist <- gridExtra::arrangeGrob(ggplot_dist, ncol = 1)
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
