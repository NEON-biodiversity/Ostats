#' Calculate multivariate O-statistics
#'
#' Notes: (REPLACE THIS WITH ACTUAL DOCUMENTATION LATER)
#' I decided to make this a separate function for now, but it could potentially be merged into the main Ostats() function later
#' However, the data structures are somewhat different so it would take a bit of coding to merge them.
#' Currently we do not support circular multivariate data.
#' Also note that the argument hypervolume_args passes additional arguments to hypervolume::hypervolume.
#' That corresponds to the argument density_args in the function Ostats() --- which passes them to density().
#'
#' @param traits matrix of trait measurements. The number of rows in the matrix
#'   is the number of individuals,
#'   and the number of columns of the matrix is the number of traits.
#' @param plots a factor with length equal to nrow(traits) that indicates the
#'   community each individual belongs to.
#' @param sp a factor with length equal to nrow(traits) that indicates the taxon
#'   of each individual.
#' @param output specifies whether median or mean is calculated.
#' @param weight_type specifies weights to be used to calculate the median or mean.
#' @param run_null_model whether to run a null model (if \code{TRUE}) and evaluate the
#'  O-statistics against it, or simply return the raw O-statistics (if \code{FALSE}).
#'  Defaults to \code{TRUE}.
#' @param nperm the number of null model permutations to generate. Defaults to 99.
#' @param nullqs numeric vector of probabilities with values in [0,1] to set
#'   effect size quantiles. Defaults to \code{c(0.025, 0.975)}.
#' @param shuffle_weights If TRUE, shuffle weights given to pairwise overlaps
#'   within a community when generating null models.
#' @param swap_means If TRUE, swap means of body sizes within a community when
#'   generating null models.
#' @param random_seed User may supply a random seed to enable reproducibility
#'   of null model output. A warning is issued, and a random seed is generated
#'   based on the local time, if the user does not supply a seed.
#' @param hypervolume_args additional arguments to pass to \code{\link[hypervolume]{hypervolume}},
#' such as \code{method} If none are provided, default values
#' are used.
#'
#' @details TBD (can copy from Ostats)
#'
#' @return The function returns a list containing 4 objects:
#' \item{overlaps_norm}{a matrix showing the median of weighted pairwise overlaps (at default) of trait
#'   distributions of all species in each community (at default output
#'   and weight_type), with the area under all density functions normalized to 1.}
#' \item{overlaps_unnorm}{a matrix showing O-stats calculated with the area under all density
#'   functions proportional to the number of observations in that group.}
#' \item{overlaps_norm_ses}{5 matrices of effect size statistics against a null model
#'   with the area under all density functions normalized to 1.}
#' \item{overlaps_unnorm_ses}{5 matrices of effect size statistics against a null model
#'   with the area under all density functions proportional to the number
#'   of observations in that group.}
#'
#' @seealso \code{\link{Ostats}} for univariate data.
#'
#' @export
Ostats_multivariate <- function(traits, plots, sp, output = "median", weight_type = "hmean", run_null_model = TRUE, nperm = 99, nullqs = c(0.025, 0.975), shuffle_weights = FALSE, swap_means = FALSE, random_seed = NULL, hypervolume_args = list()) {
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits),
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits),

  # warning and error messages to check the inputs
  if(is.numeric(traits) == FALSE) stop("traits must be a numeric matrix.")
  if(length(unique(sp)) == 1) warning("only one taxon is present; overlap cannot be calculated.")

  # If user did not supply a random seed, generate one and print a message.
  if (run_null_model && missing(random_seed)) {
    random_seed <- round(as.numeric(Sys.time()) %% 12345)
    message(paste("Note: argument random_seed was not supplied; setting seed to", random_seed))
  }

  # If the abundances of species are different, print a message.
  if (length(unique(table(sp))) > 1) message("Note: species abundances differ. Consider sampling equivalent numbers of individuals per species.")

  # Set random seed.
  set.seed(random_seed)

  # Declaration of data structures to hold the results
  # Data structures for observed O-Stats
  overlaps_norm <- matrix(nrow = nlevels(plots), ncol = 1)
  overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = 1)

  # Data structures for null O-Stats
  overlaps_norm_null <- array(dim = c(nlevels(plots), 1, nperm))
  overlaps_unnorm_null <- array(dim = c(nlevels(plots), 1, nperm))

  # Name rows of the outputs (in multivariate case we do not name columns after traits because there is only one column in output).
  dimnames(overlaps_norm) <- list(levels(plots))
  dimnames(overlaps_unnorm) <- list(levels(plots))

  # Calculation of observed O-Stats

  print('Calculating observed local O-stats for each community . . .')
  pb <- utils::txtProgressBar(min = 0, max = nlevels(plots), style = 3)

  for (s in 1:nlevels(plots)) {
    overlap_norm_s <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], output = output, weight_type = weight_type, normal=TRUE, density_args = hypervolume_args), TRUE)
    overlaps_norm[s, 1] <- if (inherits(overlap_norm_s, 'try-error')) NA else overlap_norm_s
    overlap_unnorm_s <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], output = output, weight_type = weight_type, normal = FALSE, density_args = hypervolume_args), TRUE)
    overlaps_unnorm[s, 1] <- if (inherits(overlap_unnorm_s, 'try-error')) NA else overlap_unnorm_s

    utils::setTxtProgressBar(pb, s)
  }

  close(pb)

  if (run_null_model) {

    print('Calculating null distributions of O-stats . . . ')
    pb <- utils::txtProgressBar(min = 0, max = nperm, style = 3)

    # Null model generation and calculation of null O-Stats

    # Local null model: generation and calculation done in the same loop
    for (i in 1:nperm) {
      utils::setTxtProgressBar(pb, i)
      for (s in 1:nlevels(plots)) {

        if (shuffle_weights == FALSE & swap_means == FALSE) overlap_norm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sample(sp[plots == levels(plots)[s]]), output = output, weight_type = weight_type, normal=TRUE, density_args = hypervolume_args), TRUE)
        if (shuffle_weights == TRUE) overlap_norm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], output = output, weight_type = weight_type, normal=TRUE, randomize_weights = TRUE, density_args = hypervolume_args), TRUE)
        if (swap_means == TRUE) {
          traits_s <- traits[plots==levels(plots)[s], ]
          sp_s <- sp[plots==levels(plots)[s]]

          traitmeans <- tapply(traits_s,sp_s,mean)
          traitdeviations <- traits_s-traitmeans[sp_s]

          # Sort the trait means out randomly.
          traitmeans_null <- sample(traitmeans)
          sp_null <- rep(names(traitmeans_null), table(sp_s))
          traits_null <- traitdeviations + traitmeans_null[sp_null]
          overlap_norm_si <- try(community_overlap_merged(traits = traits_null, sp = sp_null, output = output, weight_type = weight_type, normal=TRUE, randomize_weights = FALSE, density_args = hypervolume_args), TRUE)
        }

        overlaps_norm_null[s, 1, i] <- if (inherits(overlap_norm_si, 'try-error')) NA else overlap_norm_si
        overlap_unnorm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sample(sp[plots == levels(plots)[s]]), output = output, weight_type = weight_type, normal=FALSE, density_args = hypervolume_args), TRUE)
        overlaps_unnorm_null[s, 1, i] <- if (inherits(overlap_unnorm_si, 'try-error')) NA else overlap_unnorm_si

      }
    }

    close(pb)
    print('Extracting null quantiles to get standardized effect sizes (almost done!) . . .')

    # Extract quantiles to get standardized effect sizes for the overlap stats
    overlaps_norm_ses <- get_ses(overlaps_norm, overlaps_norm_null, nullqs)
    overlaps_unnorm_ses <- get_ses(overlaps_unnorm, overlaps_unnorm_null, nullqs)
    list(overlaps_norm=overlaps_norm, overlaps_unnorm=overlaps_unnorm,
         overlaps_norm_ses=overlaps_norm_ses, overlaps_unnorm_ses=overlaps_unnorm_ses)

  } else {
    list(overlaps_norm=overlaps_norm, overlaps_unnorm=overlaps_unnorm)
  }

}
