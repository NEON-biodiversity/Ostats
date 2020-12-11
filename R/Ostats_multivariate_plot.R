#' Plot multivariate community overlap
#'
#' @description Placeholder.
#'
#' @param plots Site identity: a vector of names of each community.
#' @param sp Taxon identity: a vector of species or taxa names.
#' @param traits A matrix or data frame with rows representing individuals and columns representing traits.
#' @param overlap_dat Optional: an object containing the output of \code{\link{Ostats_multivariate}}. If provided, overlap statistics will be displayed in the plot panels.
#' @param use_plots a vector of sites to plot. If NULL, the function will plot all the sites.
#' @param colorvalues Vector of color values for the density polygons. Defaults to a viridis palette if none provided.
#' @param alpha defines the transparency level for the density polygons. Default is 0.5.
#' @param panel_height height of the individual plot panels, in units given by \code{units}. Default is 3 cm.
#' @param panel_width height of the individual plot panels, in units given by \code{units}. Default is 3 cm.
#' @param units units for panel height and width. Default is centimeters.
#' @param hypervolume_args additional arguments to pass to \code{\link[hypervolume]{hypervolume}},
#' such as \code{method} If none are provided, default values
#' are used.
#'
#' MORE DOCUMENTATION GOES HERE
#'
#' @export
Ostats_multivariate_plot <- function() {
  # Plot 2x2 dimensions of the hypervolumes separately, a la hypervolume plot
  # Species will get colors.
  # Do this for each community on a separate page.
  # It will return a list if there are >1 communities.

  # Do not plot points, only the hypervolumes.

  # Panel widths and heights
  units <- 'cm'
  panel_height <- grid::unit(3, units = units)
  panel_width <- grid::unit(3, units = units)

  trait_combs <- combn(names(traits), 2) # All combinations of traits

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
  color_scale <- ggplot2::scale_color_manual(values = setNames(colorvalues, sp_names))

  # Generate common legend for all plots by writing species names in different colors.
  legend_panel <- ggplot2::ggplot(data.frame(x = 1, y = seq_along(sp_names), sp = sp_names),
                                  ggplot2::aes(x = x, y = y, label = sp, color = sp)) +
    ggplot2::geom_text(hjust = 0) +
    color_scale +
    ggplot2::scale_y_continuous(expand = c(0.5, 0.5)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = 'none')

  plot_list <- list()

  for (p in unique(plots)) {

    sp_plot <- sp[plots == p]
    traits_plot <- traits[plots == p, ]

    plot_theme <- ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'none')

    # Generate hypervolumes for each species at this plot.
    # Hypervolume for each species/site combination.
    sp_in_plot <- unique(sp_plot)
    hv_list <- replicate(length(sp_in_plot), NA, simplify = FALSE)

    for (i in 1:length(sp_in_plot)) {
      sp_plot_dat <- traits_plot[sp_plot == sp_in_plot[i], ]
      # Only generate hypervolume with at least 3 measurements
      # Generate UNSCALED hypervolume.
      if (nrow(sp_plot_dat) > 2) {
        hv_list[[i]] <- hypervolume::hypervolume(sp_plot_dat, method = 'gaussian')
      }
    }

    # Generate contours for all hypervolumes
    contours_list <- lapply(hv_list, function(hv) if (class(hv) == 'Hypervolume') get_contours(hv) else NA)
    # Join contours to data frame
    for (i in 1:length(contours_list)) contours_list[[i]][, 'sp'] = sp_in_plot[i]
    contours_df <- do.call(rbind, contours_list)

    # Generate plots, with contours
    trait_pairs_plot_list <- list()
    for (i in 1:ncol(trait_combs)) {

      dat_points <- data.frame(sp = sp_plot, x = traits_plot[, trait_combs[1, i]], y = traits_plot[, trait_combs[2, i]])
      dat_polygons <- contours_df[contours_df$trait_x == trait_combs[1, i] & contours_df$trait_y == trait_combs[2, i], ]

      trait_pairs_plot_list[[i]] <- ggplot2::ggplot(dat_points, ggplot2::aes(x = x, y = y, group = sp, color = sp)) +
        ggplot2::geom_polygon(data = dat_polygons, ggplot2::aes(group = interaction(sp, polygon_id)), fill = 'transparent') +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(limits = range(contours_df$x[contours_df$trait_x == trait_combs[1, i]])) +
        ggplot2::scale_y_continuous(limits = range(contours_df$y[contours_df$trait_y == trait_combs[2, i]])) +
        color_scale +
        plot_theme +
        ggplot2::labs(x = trait_combs[1, i], y = trait_combs[2, i])
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

    plot_list[[length(plot_list) + 1]] <- ggpubr::as_ggplot(trait_pairs_plots_arranged)

  }

  names(plot_list) <- unique(plots)

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }


}

#' Unexported function to draw contours somewhat modified from hypervolume::plot.HypervolumeList
get_contours <- function(hv) {
  hv_density <- nrow(hv@RandomPoints)/hv@Volume
  hv_dimensionality <- hv@Dimensionality
  radius_critical <- hv_density^(-1/hv_dimensionality)
  # Calculate kernel density estimate for each combinations of two variables.
  contour_list <- list()
  for (i in 1:ncol(trait_combs)) {
    kde <- MASS::kde2d(hv@RandomPoints[, trait_combs[1, i]], hv@RandomPoints[, trait_combs[2, i]], n = 50, h = radius_critical, lims = c(range(hv@RandomPoints[, trait_combs[1, i]]) * c(0.9, 1.1), range(hv@RandomPoints[, trait_combs[2, i]]) * c(0.9, 1.1)))
    contour_lines <- contourLines(kde, levels = 0.01)
    contour_line_dfs <- list()
    for (j in 1:length(contour_lines)) {
      contour_line_dfs[[j]] <- with(contour_lines[[j]], data.frame(polygon_id = j, x = x, y = y))
    }
    contour_list[[i]] <- data.frame(trait_x = trait_combs[1, i], trait_y = trait_combs[2, i], do.call(rbind, contour_line_dfs))
  }
  do.call(rbind, contour_list)
}


