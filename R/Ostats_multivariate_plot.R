#' Plot multivariate community overlap
#'
#' @description This function plots the overlap of traits among
#'  species for each community in multivariate space, showing
#'  projections of trait hypervolumes into two-dimensional space for
#'  all pairs of traits.
#'
#' @param plots Site identity: a vector of names of each community.
#' @param sp Taxon identity: a vector of species or taxa names.
#' @param traits A matrix or data frame with rows representing individuals
#'   and columns representing traits.
#' @param overlap_dat Optional: an object containing the output of
#'   \code{\link{Ostats_multivariate}}.
#'   If provided, overlap statistics will be displayed in the plot panels.
#' @param use_plots a vector of sites to plot. If NULL, the function will plot all the sites.
#' @param colorvalues Vector of color values for the density polygons.
#'   Defaults to a viridis palette if none provided.
#' @param plot_points whether to plot individual data points in addition to
#'   the hypervolume slices. Default is TRUE.
#' @param contour_level level at which to plot contour lines. If not provided
#'  by the user, a message is issued stating the default plotting level is 0.01.
#' @param axis_expansion multiplicative expansion factor by which to expand the x and y axes around
#'  the hypervolume contours before plotting. Default is 0.01.
#' @param contour_buffer_factor multiplicative expansion factor by which to expand the x and y axes in all directions,
#'   relative to the range of the axis, before calculating the hypervolume contours for plotting.
#'   If this is not set to a sufficiently large value, the contour lines of the hypervolumes will be cut off.
#'   Default value is 0.25 (25\% expansion of the axis limits in all directions).
#' @param panel_height height of the individual plot panels, in units given by \code{units}. Default is 3 cm.
#' @param panel_width height of the individual plot panels, in units given by \code{units}. Default is 3 cm.
#' @param units units for panel height and width. Default is centimeters.
#' @param hypervolume_args additional arguments to pass to \code{\link[hypervolume]{hypervolume}},
#'   such as \code{method}. If none are provided, default values are used.
#'
#' @details Some of the code for generating contour lines is modified from
#'   \code{\link[hypervolume:plot.Hypervolume]{plot.HypervolumeList}}.
#'
#' @return Two-dimensional projections of species trait hypervolumes for each pair of traits,
#'  plotted together for each community to show how they overlap each other.
#'  The overlap value obtained as output from \code{\link{Ostats_multivariate}}
#'  is labelled on each community graph, if provided by the user.
#'
#'  The class of the returned object is \code{Ostats_plot_object}. Calling
#'  \code{print} on this object will invoke a method to draw the plot using
#'  \code{\link[grid]{grid.draw}}.
#'
#'  If more than one community is provided, a list of objects of class
#'  \code{Ostats_plot_object} will be returned.
#'
#' @seealso \code{\link{Ostats_multivariate}} for generating multivariate O-statistics
#'
#' @export
Ostats_multivariate_plot <- function(plots,
                                     sp,
                                     traits,
                                     overlap_dat = NULL,
                                     use_plots = NULL,
                                     colorvalues = NULL,
                                     plot_points = TRUE,
                                     contour_level,
                                     axis_expansion = 0.01,
                                     contour_buffer_factor = 0.25,
                                     panel_height = 3,
                                     panel_width = 3,
                                     units = 'cm',
                                     hypervolume_args = list()) {

  if (missing(contour_level)) {
    contour_level <- 0.01
    message(paste('User did not specify density function value at which to plot contours. Using contour_level =', contour_level))
  }

  message('Estimating hypervolumes for plots. This may take a few minutes. . .')

  # Process input data
  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(use_plots)) {
    use_plots <- unique(plots)
  }

  # Filter overlap statistics only for use_plots
  if (!is.null(overlap_dat)) {
    ostat_norm <- overlap_dat$overlaps_norm
    ostat_norm <- ostat_norm[rownames(ostat_norm) %in% use_plots, , drop = FALSE]
  }

  # Use default method argument to hypervolume::hypervolume if none are provided
  if (!'method' %in% names(hypervolume_args)) {
    hypervolume_args[['method']] <- 'gaussian'
  }
  # Also set to verbose = FALSE if that argument is not provided.
  if (!'verbose' %in% names(hypervolume_args)) {
    hypervolume_args[['verbose']] <- FALSE
  }

  plots <- plots[plots %in% use_plots]
  sp <- sp[plots %in% use_plots]
  traits <- traits[plots %in% use_plots, ]

  # Panel widths and heights
  panel_height <- grid::unit(panel_height, units = units)
  panel_width <- grid::unit(panel_width, units = units)

  # Supply trait names if none exist
  if (is.null(dimnames(traits)[[2]])) {
    trait_names <- paste0('trait', 1:ncol(traits))
  } else {
    trait_names <- dimnames(traits)[[2]]
  }

  trait_combs <- utils::combn(trait_names, 2) # All combinations of traits

  # Create triangular layout
  layout_mat <- matrix(as.numeric(NA), ncol(traits) - 1, ncol(traits) - 1)
  layout_mat[lower.tri(layout_mat, diag = TRUE)] <- 1:ncol(trait_combs)
  layout_mat <- layout_mat[nrow(layout_mat):1, ncol(layout_mat):1]
  layout_mat[nrow(layout_mat), 1] <- ncol(trait_combs) + 1 # Location of legend

  # If no color values are provided, produce default colors.
  if (is.null(colorvalues)) {
    colorvalues <- sample(viridis::viridis(length(unique(sp)), alpha = 1))
  }

  sp_names <- rev(sort(unique(sp)))
  color_scale <- ggplot2::scale_color_manual(values = stats::setNames(colorvalues, sp_names))

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none')

  plot_list <- list()

  for (p in unique(plots)) {

    # Concatenate plot name and overlap statistic for the plot.
    if (is.null(overlap_dat)) {
      plot_label <- p
    } else {
      plot_label <- paste0(p, ': (overlap = ', round(ostat_norm[which(unique(plots) == p)], 2), ')')
    }

    # Generate legend for plot by drawing colored points.
    legend_panel <- ggplot2::ggplot(data.frame(x = 0.05, y = seq_along(sp_names), sp = sp_names),
                                    ggplot2::aes(x = x, y = y, label = sp)) +
      ggplot2::geom_text(hjust = 0) +
      ggplot2::geom_text(data = data.frame(x = 0.05, y = length(sp_names) + 1, label = plot_label), ggplot2::aes(label = label), hjust = 0, fontface = 'bold') +
      ggplot2::geom_point(x = 0, mapping = ggplot2::aes(color = sp), size = 3) +
      color_scale +
      ggplot2::scale_y_continuous(expand = c(0.5, 0.5)) +
      ggplot2::scale_x_continuous(limits = c(-0.1, 1), expand = ggplot2::expansion(mult = 0.01)) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = 'none')


    sp_plot <- sp[plots == p]
    traits_plot <- traits[plots == p, ]

    # Generate hypervolumes for each species at this plot.
    # Hypervolume for each species/site combination.
    sp_in_plot <- unique(sp_plot)
    hv_list <- replicate(length(sp_in_plot), NA, simplify = FALSE)

    for (i in 1:length(sp_in_plot)) {
      sp_plot_dat <- traits_plot[sp_plot == sp_in_plot[i], ]

      if (nrow(sp_plot_dat) > 2) {
        # Suppress all progress messages from hypervolume functions, including those from underlying C functions.
        invisible(utils::capture.output(suppressWarnings(suppressMessages({
          hv_list[[i]] <- do.call(hypervolume::hypervolume, args = c(list(data = sp_plot_dat), hypervolume_args))
        }))))
      }
    }

    # Generate contours for all hypervolumes
    contours_list <- lapply(hv_list, function(hv) if (class(hv) == 'Hypervolume') get_contours(hv, trait_combs, contour_buffer_factor, contour_level) else NA)
    # Join contours to data frame
    for (i in 1:length(contours_list)) contours_list[[i]][, 'sp'] = sp_in_plot[i]
    contours_df <- do.call(rbind, contours_list)

    # Generate plots, with contours
    trait_pairs_plot_list <- list()
    for (i in 1:ncol(trait_combs)) {

      dat_points <- data.frame(sp = sp_plot, x = traits_plot[, trait_combs[1, i]], y = traits_plot[, trait_combs[2, i]])
      dat_polygons <- contours_df[contours_df$trait_x == trait_combs[1, i] & contours_df$trait_y == trait_combs[2, i], ]

      plot_i <- ggplot2::ggplot() +
        ggplot2::geom_polygon(data = dat_polygons, ggplot2::aes(x = x, y = y, group = interaction(sp, polygon_id), color = sp), fill = 'transparent') +
        ggplot2::coord_cartesian(xlim = range(contours_df$x[contours_df$trait_x == trait_combs[1, i]]),
                                 ylim = range(contours_df$y[contours_df$trait_y == trait_combs[2, i]])) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = axis_expansion)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = axis_expansion)) +
        color_scale +
        plot_theme +
        ggplot2::labs(x = trait_combs[1, i], y = trait_combs[2, i])

      if (plot_points) plot_i <- plot_i + ggplot2::geom_point(data = dat_points, ggplot2::aes(x = x, y = y, group = sp, color = sp))

      trait_pairs_plot_list[[i]] <- plot_i

    }

    # Remove axis text and titles from plots not along the edge.
    for (i in layout_mat[upper.tri(layout_mat)]) {
      trait_pairs_plot_list[[i]] <- trait_pairs_plot_list[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                                axis.text.y = ggplot2::element_blank(),
                                                                                axis.title.x = ggplot2::element_blank(),
                                                                                axis.title.y = ggplot2::element_blank())
    }

    # Add legend panel to plot list
    trait_pairs_plot_list <- c(trait_pairs_plot_list, list(legend_panel))

    # Arrange plots
    # Sequentially bind the rows together, then the columns.
    trait_pairs_plots_rows <- list()

    # Create dummy grob to fill in the blank spaces, using a zeroGrob
    dummy_grob <- ggplot2::ggplotGrob(trait_pairs_plot_list[[1]])
    dummy_grob$grobs <- replicate(length(dummy_grob$grobs), ggplot2::zeroGrob(), simplify = FALSE)


    for (i in 1:nrow(layout_mat)) {
      trait_pairs_plots_rows[[i]] <- do.call(gridExtra::gtable_cbind, lapply(layout_mat[i, ], function(n) if (is.na(n)) dummy_grob else ggplot2::ggplotGrob(trait_pairs_plot_list[[n]])))
    }

    trait_pairs_plots_arranged <- do.call(gridExtra::gtable_rbind, trait_pairs_plots_rows)

    # Assign class attribute Ostats_plot_object so that the plot has a default print method.
    attr(trait_pairs_plots_arranged, 'class') <- c('Ostats_plot_object', attr(trait_pairs_plots_arranged, 'class'))

    plot_list[[length(plot_list) + 1]] <- trait_pairs_plots_arranged

  }

  names(plot_list) <- unique(plots)

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }

}

#' Unexported function to draw contours somewhat modified from hypervolume::plot.HypervolumeList
#' @noRd
get_contours <- function(hv, trait_combs, contour_buffer_factor, contour_level) {
  hv_density <- nrow(hv@RandomPoints)/hv@Volume
  hv_dimensionality <- hv@Dimensionality
  radius_critical <- hv_density^(-1/hv_dimensionality)
  # Calculate kernel density estimate for each combinations of two variables.
  contour_list <- list()
  for (i in 1:ncol(trait_combs)) {
    range_x <- range(hv@RandomPoints[, trait_combs[1, i]])
    range_x_buffered <- range_x + c(-1, 1) * contour_buffer_factor * diff(range_x)
    range_y <- range(hv@RandomPoints[, trait_combs[2, i]])
    range_y_buffered <- range_y + c(-1, 1) * contour_buffer_factor * diff(range_y)
    kde <- MASS::kde2d(hv@RandomPoints[, trait_combs[1, i]], hv@RandomPoints[, trait_combs[2, i]], n = 50, h = radius_critical, lims = c(range_x_buffered, range_y_buffered))
    contour_lines <- grDevices::contourLines(kde, levels = contour_level)
    contour_line_dfs <- list()
    for (j in 1:length(contour_lines)) {
      contour_line_dfs[[j]] <- with(contour_lines[[j]], data.frame(polygon_id = j, x = x, y = y))
    }
    contour_list[[i]] <- data.frame(trait_x = trait_combs[1, i], trait_y = trait_combs[2, i], do.call(rbind, contour_line_dfs))
  }
  do.call(rbind, contour_list)
}


