#' Community Overlap Calculation
#'
#' This function calculates the median or mean of pairwise overlaps between density
#' estimates of trait distributions of all species within a community, which can
#' be weighted by species abundances.
#'
#' @param traits a vector of trait measurements in the univariate case, or a
#'   matrix in the multivariate case where each column is a trait.
#' @param sp a vector with length equal to length(traits) that indicates the
#' taxon of each individual.
#' @param discrete whether trait data may take continuous or discrete values. Defaults to
#'   \code{FALSE} (all traits continuous). A single logical value or a logical
#'   vector with length equal to the number of columns in traits.
#' @param circular whether trait data are circular (e.g., hours or angles). Defaults to
#'   \code{FALSE} (all traits non-circular). A single logical value or a logical
#'   vector with length equal to the number of columns in traits.
#' @param normal if TRUE, the area under all density functions is normalized to 1,
#' if FALSE, the area under all density functions is proportional to the number of
#' observations in that group.
#' @param output specifies whether median or mean is calculated.
#' @param weight_type specifies weights to be used to calculate the median or mean.
#' @param randomize_weights If TRUE, randomize weights given to pairwise overlaps
#'   within a community. This can be used to generate null models.
#' @param unique_values Vector of all possible discrete values that \code{traits}
#'   can take. Only used if \code{discrete = TRUE} and \code{circular = TRUE}.
#' @param raw If \code{TRUE}, also return the raw individual pairwise overlaps
#'   used to calculate the community-level statistic. Default is \code{FALSE}.
#' @param circular_args optional list of additional arguments to pass to
#'  \code{\link[circular]{circular}}. Only used if \code{circular = TRUE} and
#'  \code{discrete = FALSE}.
#' @param density_args list of additional arguments to be passed to
#'   \code{\link[stats]{density}} if univariate, or
#'   \code{\link[hypervolume]{hypervolume}} if multivariate.
#' @param hypervolume_set_args list of additional arguments to be passed to
#'   \code{\link[hypervolume]{hypervolume_set}}. Used only in multivariate case.
#'
#' @details The function evaluates weighted mean or median of overlaps of density
#' estimates of all species in a community taking complete cases with species abundances
#' greater than 1 from the dataset. The default calculates the median of pairwise overlaps
#' for the whole community using the harmonic means of abundances of the species pairs as
#' weights, which minimizes the effect of outliners and rare species.If the argument
#' \code{weight_type = "none"}, no weights are used for the calculation of mean/median. If
#' \code{weight_type = "mean"}, arithmetic means of abundances are used as weights. To change the
#' output to mean, specify the argument \code{output = "mean"}.
#'
#' @return The function returns the O-statistic for the community as a numeric value. If
#' \code{raw = TRUE}, instead a list is returned, where the first element \code{value} is
#' the numeric value, and the second element \code{raw} is a data frame with all the raw
#' pairwise overlaps.
#'
#' @references Read, Q. D. et al. Among-species overlap in rodent body size
#'   distributions predicts species richness along a temperature gradient.
#'   Ecography 41, 1718-1727 (2018).
#'
#' @examples
#' library(Ostats)
#'
#' # Keep only the relevant part of small mammal data
#' dat <- small_mammal_data[small_mammal_data$siteID %in% c('HARV','JORN'), ]
#' dat <- dat[!is.na(dat$weight), ]
#' dat$log_weight <- log10(dat$weight)
#'
#' # Calculate median of pairwise overlaps for the community,weighted by harmonic means
#' # of abundances
#' community_overlap(traits = as.matrix(dat$log_weight),
#'    sp = factor(dat$taxonID))
#'
#' @export
community_overlap <- function(traits, sp, discrete = FALSE, circular = FALSE, normal = TRUE, output = "median", weight_type= "hmean", randomize_weights = FALSE, unique_values = NULL, raw = FALSE, circular_args = list(), density_args = list(), hypervolume_set_args = list()) {

  # Return error if circular or discrete are specified with multivariate data.
  if (circular & 'matrix' %in% class(traits)) {
    stop("circular data types are not supported with multivariate data.")
  }
  if (discrete & 'matrix' %in% class(traits)) {
    stop("discrete data types are not supported with multivariate data.")
  }

  # Clean input, removing missing values and species with <2 values.
  sp <- as.character(sp)
  dat <- cbind(as.data.frame(traits), sp = sp)
  dat <- dat[stats::complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat[, -ncol(dat)], dat$sp)
  uniquespp <- unique(dat$sp)
  nspp <- length(uniquespp)

  # Overlap cannot be calculated if there are less than 2 species with at least 2 individuals each.
  if (nspp < 2) return(NA)

  # Define common grid limits so that all density functions and hypervolumes are estimated across the same domain.
  grid_limits <- apply(dat[, -ncol(dat), drop = FALSE], 2, range)

  # Add a multiplicative factor to the upper and lower end of the range so that the tails aren't cut off.
  extend_grid <- c(-0.5, 0.5) %*% t(apply(grid_limits,2,diff))
  grid_limits <- grid_limits + extend_grid

  # In the discrete case for non-circular data, find the common unique values across which to build the histogram.
  if (discrete & !circular) {
    unique_values <- sort(unique(unlist(dat[, -ncol(dat)])))
  }

  # Create a list of univariate density functions (in univariate case) or hypervolumes (in multivariate case)
  density_list <- lapply(traitlist, function(x) trait_density(x, grid_limits, normal, discrete, circular, unique_values, density_args, circular_args))

  overlaps <- NULL
  abund_pairs <- NULL

  # All possible pairwise combinations.
  combs <- utils::combn(1:nspp, 2)

  for (idx in 1:ncol(combs)) {

      overlaps <- c(overlaps, pairwise_overlap(density_list[[combs[1, idx]]], density_list[[combs[2, idx]]], discrete, density_args, hypervolume_set_args))
      if (weight_type == "hmean")
        abund_pairs <- c(abund_pairs, 2/(1/abunds[combs[1, idx]] + 1/abunds[combs[2, idx]]))
      if (weight_type == "mean")
        abund_pairs <- c(abund_pairs, (abunds[combs[1, idx]] + abunds[combs[2, idx]]))

  }

  if (randomize_weights == TRUE){
    abund_pairs <- sample(abund_pairs)}
  if (output == "median" && weight_type == "none"){
    final_output <- stats::median(overlaps)}
  if (output == "median" && weight_type == "hmean"){
    final_output <- matrixStats::weightedMedian(x = as.vector(overlaps), w = abund_pairs)}
  if (output == "median" && weight_type == "mean"){
    final_output <- matrixStats::weightedMedian(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "hmean"){
    final_output <- stats::weighted.mean(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "mean"){
    final_output <- stats::weighted.mean(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "none"){
    final_output <- mean(overlaps)}

  if (raw == TRUE) {
    # Create a data frame of raw overlaps to return as well.
    overlaps_df <- as.data.frame(cbind(t(combs), overlaps))
    overlaps_df[,1] <- uniquespp[overlaps_df[,1]]
    overlaps_df[,2] <- uniquespp[overlaps_df[,2]]
    names(overlaps_df) <- c('species1', 'species2', 'overlap')

    return(list(value = final_output, raw = overlaps_df))
  } else {
    return(final_output)
  }
}


