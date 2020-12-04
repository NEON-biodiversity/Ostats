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
