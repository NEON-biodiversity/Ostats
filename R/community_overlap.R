#' Community Overlap Calculation
#'
#' This function calculates the median or mean of pairwise overlaps between density
#' estimates of trait distributions of all species within a community, which can
#' be weighted by species abundances.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#' taxon of each individual.
#' @param data_type data type can be "linear", "circular",or "circular_discrete".
#' Default to "linear".
#' @param normal if TRUE, the area under all density functions is normalized to 1,
#' if FALSE, the area under all density functions is proportional to the number of
#' observations in that group.
#' @param output specifies whether median or mean is calculated.
#' @param weight_type specifies weights to be used to calculate the median or mean.
#' @param randomize_weights If TRUE, randomize weights given to pairwise overlaps
#'   within a community. This can be used to generate null models.
#' @param circular_args list of additional arguments to be passed to
#'   \code{\link[circular]{circular}}. Only used if \code{data_type} is "circular".
#' @param density_args list of additional arguments to be passed to
#'   \code{\link[stats]{density}}.
#'
#' @details The function evaluates weighted mean or median of overlaps of density
#' estimates of all species in a community taking complete cases with species abundances
#' greater than 1 from the dataset. The default calculates the median of pairwise overlaps
#' for the whole community using the harmonic means of abundances of the species pairs as
#' weights, which minimizes the effect of outliners and rare species.If the argument
#' weight_type = "none", no weights are used for the calculation of mean/median. If
#' weight_type = "mean", arithmetic means of abundances are used as weights. To change the
#' output to mean, specify the argument output = "mean".
#'
#' @return The function returns overall species overlap for a community. At default, it returns the
#' median of pairwise overlaps weighted by harmonic means of abundances for the community.
#'
#' @references Read, Q. D. et al. Among-species overlap in rodent body size
#'   distributions predicts species richness along a temperature gradient.
#'   Ecography 41, 1718-1727 (2018).
#'
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical
#' density estimates.
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
community_overlap <- function(traits, sp, data_type = "linear", normal = TRUE, output = "median", weight_type= "hmean", randomize_weights = FALSE, circular_args = list(), density_args = list()) {

  # Return error if circular is specified with multivariate data.
  if (data_type %in% c('circular', 'circular_discrete') & 'matrix' %in% class(traits)) {
    stop("circular data types are not supported with multivariate data.")
  }

  # Clean input, removing missing values and species with <2 values.
  sp <- as.character(sp)
  dat <- cbind(as.data.frame(traits), sp = sp)
  dat <- dat[stats::complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat[, -ncol(dat)], dat$sp)
  nspp <- length(traitlist)

  # Overlap cannot be calculated if there are less than 2 species with at least 2 individuals each.
  if (nspp < 2) return(NA)

  # Define common grid limits so that all density functions and hypervolumes are estimated across the same domain.
  grid_limits <- apply(dat[, -ncol(dat), drop = FALSE], 2, range)

  # Add a multiplicative factor to the upper and lower end of the range so that the tails aren't cut off.
  extend_grid <- c(-0.5, 0.5) %*% t(apply(grid_limits,2,diff))
  grid_limits <- grid_limits + extend_grid

  # Create a list of univariate density functions (in univariate case) or hypervolumes (in multivariate case)
  density_list <- lapply(traitlist, function(x) trait_density(x, grid_limits, normal, data_type, density_args, circular_args))

  overlaps <- NULL
  abund_pairs <- NULL

  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {

      overlaps <- c(overlaps, pairwise_overlap(density_list[[sp_a]], density_list[[sp_b]], density_args))
      if (weight_type == "hmean")
        abund_pairs <- c(abund_pairs, 2/(1/abunds[sp_a] + 1/abunds[sp_b]))
      if (weight_type == "mean")
        abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))

    }
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

return(final_output)
}


