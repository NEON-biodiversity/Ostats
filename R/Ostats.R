#' Calculate O-statistics (community-level pairwise niche overlap statistics)
#'
#' This is the primary function in the Ostats package. It evaluates the
#' O-statistics against a local null model.
#'
#' @seealso \code{\link{Ostats_regional}} to evaluate against a regional null model.
#'
#' @export
Ostats <- function(traits, plots, sp, nperm = 99, nullqs = c(0.025, 0.975), shuffle_weights = FALSE, swap_means = FALSE) {
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits),
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits),

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
      overlap_norm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
      overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
      overlap_unnorm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
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

#' Calculate O-statistics (community-level pairwise niche overlap statistics) evaluated against regional null model
#'
#' This function evaluates the O-statistics against a regional null model.
#' In contrast to \code{Ostats}, it requires trait and species information
#' about the regional species pool.
#'
#' @seealso \code{\link{Ostats}} to evaluate against a local null model.
#'
#' @export
Ostats_regional <-function(traits, plots, sp, reg_pool_traits, reg_pool_sp, nperm = 99, nullqs = c(0.025, 0.975)) {
	# Declaration of data structures to hold the results

	# Data structures for observed O-Stats
	overlaps_reg <- matrix(nrow = nlevels(plots), ncol = ncol(traits))

	# Data structures for null O-Stats
	allpool_overlaps_reg_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
	bysp_overlaps_reg_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))

	# Name rows and columns of the outputs.
	dimnames(overlaps_reg) <- list(levels(plots), dimnames(traits)[[2]])

	# Calculation of observed O-Stats

	print('Calculating observed regional O-stats for each community . . .')
	pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)

	for (s in 1:nlevels(plots)) {
		for (t in 1:ncol(traits)) {
			overlap_reg_st <- try(pairwise_overlap(a = traits[plots == levels(plots)[s], t], b = `[[`(reg_pool_traits, levels(plots)[s])[, t], norm=TRUE), TRUE)
			overlaps_reg[s, t] <- if(inherits(overlap_reg_st, 'try-error')) NA else overlap_reg_st[1]
		}
		setTxtProgressBar(pb, s)
	}

	close(pb)

	print('Calculating total-pool and species-pool null distributions of regional O-stats . . . ')
	pb <- txtProgressBar(min = 0, max = nperm, style = 3)

	# Null model generation and calculation of null O-Stats

	# Total-pool null model: generation and calculation done in the same loop
		for (i in 1:nperm) {
			for (s in 1:nlevels(plots)) {
				for (t in 1:ncol(traits)) {
					trnull_sti <- sample(x=`[[`(reg_pool_traits, levels(plots)[s])[, t], size=sum(plots == levels(plots)[s]), replace=FALSE)
					allpool_overlap_reg_sti <- try(pairwise_overlap(a = trnull_sti, b = `[[`(reg_pool_traits, levels(plots)[s])[, t], norm = TRUE), TRUE)
					allpool_overlaps_reg_null[s, t, i] <- if (inherits(allpool_overlap_reg_sti, 'try-error')) NA else allpool_overlap_reg_sti[1]
				}
			}

	# By-species null model: generation of null communities and calculation of O-stats (in same loop)


		for (s in 1:nlevels(plots)) {
			for (t in 1:ncol(traits)) {
				sp_sti <- sp[plots == levels(plots)[s]]
				trnull_sti <- c()
				spnames_sti <- unique(sp_sti)
				sptable_sti <- table(sp_sti)

				for (j in 1:length(spnames_sti)) {
					z <- `[[`(reg_pool_traits, levels(plots)[s])[`[[`(reg_pool_sp, levels(plots)[s]) == as.character(spnames_sti[j]), t]
					if (any(!is.na(z))) trnull_sti <- c(trnull_sti, sample(x=z, size=sum(sp_sti == spnames_sti[j]), replace=FALSE))
				}

				bysp_overlap_reg_sti <- try(pairwise_overlap(a = trnull_sti, b = `[[`(reg_pool_traits, levels(plots)[s])[, t], norm = TRUE), TRUE)
				bysp_overlaps_reg_null[s, t, i] <- if (inherits(bysp_overlap_reg_sti, 'try-error')) NA else bysp_overlap_reg_sti[1]

			}
		}
		setTxtProgressBar(pb, i)
	}

	close(pb)
	print('Extracting null quantiles to get standardized effect sizes (almost done!) . . .')

	overlaps_reg_allpool_ses <- get_ses(overlaps_reg, allpool_overlaps_reg_null, nullqs)
	overlaps_reg_bysp_ses <- get_ses(overlaps_reg, bysp_overlaps_reg_null, nullqs)

	list(overlaps_reg=overlaps_reg, overlaps_reg_allpool_ses=overlaps_reg_allpool_ses, overlaps_reg_bysp_ses=overlaps_reg_bysp_ses)


}
