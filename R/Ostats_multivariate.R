# Test version of multivariate O-stats
# QDR 19 Oct 2020
# I decided to make this a separate function for now, but it could potentially be merged into the main Ostats() function later
# However, the data structures are somewhat different so it would take a bit of coding to merge them.
# This does not support circular data.
# Also note at the moment that the argument hypervolume_args is passed to density_args for consistency with argument names.
# It is all internal so not really an issue.

Ostats_multivariate <- function(traits, plots, sp, output = "median", weight_type= "hmean", nperm = 99, nullqs = c(0.025, 0.975), shuffle_weights = FALSE, swap_means = FALSE, hypervolume_args = list()) {
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits),
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits),

  # warning and error messages to check the inputs
  if(is.numeric(traits) == FALSE) stop("the function only evaluates numerical data")
  if(length(unique(plots)) == 1) warning("only one community is evaluated")
  if(length(unique(sp)) == 1) warning("only one taxon is present")

  # Declaration of data structures to hold the results
  # Data structures for observed O-Stats
  overlaps_norm <- matrix(nrow = nlevels(plots), ncol = 1)
  overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = 1)

  # Data structures for null O-Stats
  overlaps_norm_null <- array(dim = c(nlevels(plots), 1, nperm))
  overlaps_unnorm_null <- array(dim = c(nlevels(plots), 1, nperm))

  # Name rows of the outputs (in multivariate case we do not name columns after traits).
  dimnames(overlaps_norm) <- list(levels(plots))
  dimnames(overlaps_unnorm) <- list(levels(plots))

  # Calculation of observed O-Stats

  print('Calculating observed local O-stats for each community . . .')
  pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)

  for (s in 1:nlevels(plots)) {
      overlap_norm_s <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], data_type = data_type, output = output, weight_type = weight_type, normal=TRUE, density_args = hypervolume_args), TRUE)
      overlaps_norm[s, 1] <- if (inherits(overlap_norm_s, 'try-error')) NA else overlap_norm_s
      overlap_unnorm_s <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], data_type = data_type, output = output, weight_type = weight_type, normal=FALSE, density_args = hypervolume_args), TRUE)
      overlaps_unnorm[s, 1] <- if (inherits(overlap_unnorm_s, 'try-error')) NA else overlap_unnorm_s

    setTxtProgressBar(pb, s)
  }

  close(pb)

  print('Calculating null distributions of O-stats . . . ')
  pb <- txtProgressBar(min = 0, max = nperm, style = 3)

  # Null model generation and calculation of null O-Stats

  # Local null model: generation and calculation done in the same loop
  for (i in 1:nperm) {
    setTxtProgressBar(pb, i)
    for (s in 1:nlevels(plots)) {

        if (shuffle_weights == FALSE & swap_means == FALSE) overlap_norm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sample(sp[plots == levels(plots)[s]]), data_type=data_type, output = output, weight_type = weight_type, normal=TRUE, circular_args = circular_args, density_args = hypervolume_args), TRUE)
        if (shuffle_weights == TRUE) overlap_norm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sp[plots == levels(plots)[s]], data_type=data_type, output = output, weight_type = weight_type, normal=TRUE, randomize_weights = TRUE, circular_args = circular_args, density_args = hypervolume_args), TRUE)
        if (swap_means == TRUE) {
          traits_s <- traits[plots==levels(plots)[s], ]
          sp_s <- sp[plots==levels(plots)[s]]

          traitmeans <- tapply(traits_s,sp_s,mean)
          traitdeviations <- traits_s-traitmeans[sp_s]

          # Sort the trait means out randomly.
          traitmeans_null <- sample(traitmeans)
          sp_null <- rep(names(traitmeans_null), table(sp_s))
          traits_null <- traitdeviations + traitmeans_null[sp_null]
          overlap_norm_si <- try(community_overlap_merged(traits = traits_null, sp = sp_null, data_type=data_type, output = output, weight_type = weight_type,normal=TRUE, randomize_weights = FALSE, circular_args = circular_args, density_args = hypervolume_args), TRUE)
        }

        overlaps_norm_null[s, 1, i] <- if (inherits(overlap_norm_si, 'try-error')) NA else overlap_norm_si
        overlap_unnorm_si <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], ], sp = sample(sp[plots == levels(plots)[s]]),data_type=data_type, output = output, weight_type = weight_type, normal=FALSE, circular_args = circular_args, density_args = hypervolume_args), TRUE)
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

}
