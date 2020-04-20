#' Calculate O-statistics (community-level pairwise niche overlap statistics)
#'
#' This is the primary function in the Ostats package. It evaluates the
#' O-statistics against a local null model.
#'
#' This test version created by QDR on 20 April 2020 is more flexible in two ways:
#' First, you can choose any measure of central tendency for calculating the
#' community-level overlap statistic from the vector of pairwise overlaps.
#' Second, you can choose any function for calculating the overlap. This is relevant
#' if you are doing the overlaps in circular distributions.
#'
#' #### STILL IN DEVELOPMENT -- NOT DONE!!! -Q 4/20
#'
#'
#' @export
Ostats2 <- function(traits, plots, sp, # The data
                   overlap_function = pairwise_overlap, # Defaults to use the basic pairwise overlap fn
                   central_tendency = c('median', 'mean'), # Can be median or mean
                   abundance_weighted = TRUE, # whether the median/mean is abundance-weighted
                   weighting_type = c('harmonic', 'arithmetic'), # Whether to weight by harmonic or arithmetic mean of pair's abund.
                   nperm = 99, nullqs = c(0.025, 0.975), # Null model options
                   shuffle_weights = FALSE, swap_means = FALSE) { # Other null model options
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits),
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits).

  # Set the central tendency function
  if (!abundance_weighted) weighting_type <- 'none'

  central_tendency_fn_name <- paste0('community_overlap_',
                                     ifelse(weighting_type == 'harmonic', 'harmonic', ''),
                                     central_tendency)
  community_overlap_function <- get(central_tendency_fn_name)

  # Declaration of data structures to hold the results

  # Data structures for observed O-Stats
  overlaps_norm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
  overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))

  # Data structures for null O-Stats
  overlaps_norm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
  overlaps_unnorm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))

  # Name rows and columns of the outputs.
  dimnames(overlaps_norm) <- list(levels(plots), dimnames(traits)[[2]])
  dimnames(overlaps_unnorm) <- list(levels(plots), dimnames(traits)[[2]])

  # Calculation of observed O-Stats

  print('Calculating observed local O-stats for each community . . .')
  pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)

  for (s in 1:nlevels(plots)) {
    for (t in 1:ncol(traits)) {
      overlap_norm_st <- try(community_overlap_function(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
      overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
      overlap_unnorm_st <- try(community_overlap_function(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
      overlaps_unnorm[s, t] <- if (inherits(overlap_unnorm_st, 'try-error')) NA else overlap_unnorm_st
    }
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
      for (t in 1:ncol(traits)) {
        if (!shuffle_weights & !swap_means) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
        if (shuffle_weights) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE, randomize_weights = TRUE), TRUE)
        if (swap_means) {
          traits_st <- traits[plots==levels(plots)[s], t]
          sp_st <- sp[plots==levels(plots)[s]]

          traitmeans <- tapply(traits_st,sp_st,mean)
          traitdeviations <- traits_st-traitmeans[sp_st]

          # Sort the trait means out randomly.
          traitmeans_null <- sample(traitmeans)
          sp_null <- rep(names(traitmeans_null), table(sp_st))
          traits_null <- traitdeviations + traitmeans_null[sp_null]
          overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits_null, sp = sp_null, norm=TRUE, randomize_weights = FALSE), TRUE)
        }

        overlaps_norm_null[s, t, i] <- if (inherits(overlap_norm_sti, 'try-error')) NA else overlap_norm_sti
        overlap_unnorm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=FALSE), TRUE)
        overlaps_unnorm_null[s, t, i] <- if (inherits(overlap_unnorm_sti, 'try-error')) NA else overlap_unnorm_sti
      }
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
