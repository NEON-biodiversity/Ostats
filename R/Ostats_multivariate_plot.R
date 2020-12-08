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

# Plot only the points.

sp <- iris$Species
traits <- iris[, 1:4]
plots <- factor(rep(1:3, nrow(traits)/3))

trait_combs <- combn(names(traits), 2) # All combinations of traits

# Create triangular layout
layout_mat <- matrix(as.numeric(NA), ncol(traits) - 1, ncol(traits) - 1)
layout_mat[lower.tri(layout_mat, diag = TRUE)] <- 1:6
layout_mat <- layout_mat[nrow(layout_mat):1, ncol(layout_mat):1]

plot_list <- list()

for (p in unique(plots)) {

  sp_plot <- sp[plots == p]
  traits_plot <- traits[plots == p, ]

  # Find ranges of each trait so that panels will have a common range.
  traits_range <- apply(traits_plot, 2, range)

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none',
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())

  # Generate plots
  trait_pairs_plot_list <- apply(trait_combs, 2, function(traits_to_plot) {
    dat_plot <- data.frame(sp = sp_plot, x = traits_plot[, traits_to_plot[1]], y = traits_plot[, traits_to_plot[2]])
    x_range <- traits_range[, traits_to_plot[1]]
    y_range <- traits_range[, traits_to_plot[2]]
    ggplot2::ggplot(dat_plot, ggplot2::aes(x = x, y = y, group = sp, color = sp)) +
      ggplot2::geom_point() +
      plot_theme +
      ggplot2::labs(x = traits_to_plot[1], y = traits_to_plot[2])
  })



  # Arrange plots
  trait_pairs_plots_arranged <- gridExtra::grid.arrange(grobs = trait_pairs_plot_list, layout_matrix = layout_mat)

  plot_list[[length(plot_list) + 1]] <- trait_pairs_plots_arranged

}
