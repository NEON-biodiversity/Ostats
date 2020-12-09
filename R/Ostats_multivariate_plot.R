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


}


#### TEST CODE BELOW THIS LINE

# Panel widths and heights
units <- 'cm'
panel_height <- grid::unit(3, units = units)
panel_width <- grid::unit(3, units = units)

# Plot only the points.

sp <- iris$Species
traits <- iris[, 1:4]
plots <- factor(rep(1:3, nrow(traits)/3))

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


color_scale <- ggplot2::scale_color_manual(values = setNames(colorvalues, sp_names))

plot_list <- list()

# Generate common legend for all plots by writing species names in different colors.
sp_names <- rev(sort(unique(sp)))
legend_panel <- ggplot2::ggplot(data.frame(x = 1, y = seq_along(sp_names), sp = sp_names),
                                ggplot2::aes(x = x, y = y, label = sp, color = sp)) +
  ggplot2::geom_text(hjust = 0) +
  color_scale +
  ggplot2::scale_y_continuous(expand = c(0.5, 0.5)) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = 'none')

for (p in unique(plots)) {

  sp_plot <- sp[plots == p]
  traits_plot <- traits[plots == p, ]

  # Find ranges of each trait so that panels will have a common range.
  traits_range <- apply(traits_plot, 2, range)

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none')

  # Generate plots
  trait_pairs_plot_list <- apply(trait_combs, 2, function(traits_to_plot) {
    dat_plot <- data.frame(sp = sp_plot, x = traits_plot[, traits_to_plot[1]], y = traits_plot[, traits_to_plot[2]])
    x_range <- traits_range[, traits_to_plot[1]]
    y_range <- traits_range[, traits_to_plot[2]]
    ggplot2::ggplot(dat_plot, ggplot2::aes(x = x, y = y, group = sp, color = sp)) +
      ggplot2::geom_point() +
      color_scale +
      plot_theme +
      ggplot2::labs(x = traits_to_plot[1], y = traits_to_plot[2])
  })

  # Remove axis text and titles from plots not along the edge.
  for (i in 1:ncol(trait_combs)) {
    if (!i %in% diag(layout_mat)) {
      trait_pairs_plot_list[[i]] <- trait_pairs_plot_list[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                                axis.text.y = ggplot2::element_blank(),
                                                                                axis.title.x = ggplot2::element_blank(),
                                                                                axis.title.y = ggplot2::element_blank())

    }
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
    trait_pairs_plots_rows[[i]] <- do.call(cbind, lapply(layout_mat[i, ], function(n) if (is.na(n)) dummy_grob else ggplot2::ggplotGrob(trait_pairs_plot_list[[n]])))
  }

  trait_pairs_plots_arranged <- do.call(rbind, trait_pairs_plots_rows)

  plot_list[[length(plot_list) + 1]] <- trait_pairs_plots_arranged

}

