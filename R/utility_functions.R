# Utility functions, not exported
# -------------------------------

#' @title  Calculation of standardized effect sizes of overlap statistics against null distributions
#'
#' @description Extract quantiles to get standardized effect sizes for the overlap stats
#'
#' @details Implementation of null model that swaps means of distributions to test whether distributions are more evenly spaced than expected by chance - It is part of the Ostats function
#'
#' @seealso \code{\link{Ostats}} to Calculate O-statistics (community-level pairwise niche overlap statistics)
#'
#' @param o an matrix data containing the overlapping values of normal distribution of the observed data (the observed local O-stats for each community)
#' @param o_n an array object containing the overlapping values of normal distribution of the randomize data - (the null distributions of O-stats function)
#' @param qs quantiles of the null distribution to extract, usually 0.025 and 0.975
#'
#' @return
#' \item{ses}{a matrix containing values of standardized effect sizes for observed data of each of the locals (rows) and traits(columns) avaliated. To see if the overlap values can be expected by chance see returns ses_upper and ses_lower. If the value corresponding to the local in ses is between ses_upper and ses_lower, the overlapping between traits in community can be expected by chance.}
#' \item{ses_lower}{a matrix containing null values extracted for the quantil of qs[1] (usually 0.025)}
#' \item{ses_upper}{a matrix containing null values extracted for the quantil of qs[2] (usually 0.975)}
#' \item{raw_lower}{a matrix containing null values extracted for the quantil of qs[1] (usually 0.025) without subtracting the mean and divided by standard deviation (see equation for calculating the standardize effect size on the function code)}
#' \item{raw_upper}{a matrix containing null values extracted for the quantil of qs[2] (usually 0.975) without subtracting the mean and divided by standard deviation (see equation for calculating the standardize effect size on the function code)}
#'
#'@example
#' #import data
#' final_mammal_data <- read.csv('final_NEON_mammal_data.csv', stringsAsFactors = FALSE)
#' mammal_logweight <- final_mammal_data[,c('logweight'), drop = FALSE]
#' traits<-mammal_logweight
#' plots <- factor(final_mammal_data$siteID)
#' sp <- factor(final_mammal_data$taxonID)
#'
#' downloading what is in the github
#' devtools::install_github('NEON-biodiversity/Ostats')
#'
#' nperm<-99 #number of permutations
#' #Data structures for observed O-Stats
#' overlaps_norm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
#' overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
#'
#' #Data structures for null O-Stats
#' overlaps_norm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
#' overlaps_unnorm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
#'
#' #Name rows and columns of the outputs.
#' dimnames(overlaps_norm) <- list(levels(plots), dimnames(traits)[[2]])
#' dimnames(overlaps_unnorm) <- list(levels(plots), dimnames(traits)[[2]])
#'
#' ###set the "o" argument - the overlap between observed data
#' for (s in 1:nlevels(plots)) {
#'  for (t in 1:ncol(traits)) {
#'    overlap_norm_st <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], t], #' sp = sp[plots == levels(plots)[s]], data_type = "linear"), TRUE)
#'    overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
#'    overlap_unnorm_st <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s], #' t], sp = sp[plots == levels(plots)[s]], data_type = "linear"), TRUE)
#'    overlaps_unnorm[s, t] <- if (inherits(overlap_unnorm_st, 'try-error')) NA else overlap_unnorm_#' st
#'  }
#' }
#'
#'
#' ###set the "o_n" argument - the overlap between randomized data
#'
#' Local null model: generation and calculation done in the same loop
#' shuffle_weights=FALSE
#' swap_means=FALSE
#' for (i in 1:nperm) {
#'  for (s in 1:nlevels(plots)) {
#'    for (t in 1:ncol(traits)) {
#'      if (!shuffle_weights & !swap_means) overlap_norm_sti <- try(community_overlap_merged(traits #' = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), data_type = #' "linear"), TRUE)
#'      if (shuffle_weights) overlap_norm_sti <- try(community_overlap_merged(traits = traits[plots #' == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE, randomize_weights = TRUE, #' data_type = "linear"), TRUE)
#'      if (swap_means) {
#'        traits_st <- traits[plots==levels(plots)[s], t]
#'        sp_st <- sp[plots==levels(plots)[s]]
#'
#'        traitmeans <- tapply(traits_st,sp_st,mean)
#'        traitdeviations <- traits_st-traitmeans[sp_st]
#'
#'        # Sort the trait means out randomly.
#'        traitmeans_null <- sample(traitmeans)
#'        sp_null <- rep(names(traitmeans_null), table(sp_st))
#'        traits_null <- traitdeviations + traitmeans_null[sp_null]
#'        overlap_norm_sti <- try(community_overlap_merged(traits = traits_null, sp = sp_null, norm #' =TRUE, randomize_weights = FALSE, data_type = "linear"), TRUE)
#'      }
#'
#'      overlaps_norm_null[s, t, i] <- if (inherits(overlap_norm_sti, 'try-error')) NA else overlap_norm_sti
#'      overlap_unnorm_sti <- try(community_overlap_merged(traits = traits[plots == levels(plots)[s] #' , t], sp = sample(sp[plots == levels(plots)[s]]), data_type = "linear"), TRUE)
#'      overlaps_unnorm_null[s, t, i] <- if (inherits(overlap_unnorm_sti, 'try-error')) NA else #' overlap_unnorm_sti
#'    }
#'  }
#' }
#'
#'
#' SES<-get_ses(o = overlaps_norm,o_n =overlaps_norm_null ,qs = c(0.025, 0.975))
#'
#'
#'@export
#'
get_ses <- function(o, o_n, qs) {
  ses <- ses_lower <- ses_upper <- raw_lower <- raw_upper <- matrix(NA, nrow=nrow(o), ncol=ncol(o))

  for (i in 1:nrow(o)) {
    for (j in 1:ncol(o)) {
      if(!is.na(o[i,j])) {
        obs <- o[i,j]
        nullvals <- na.omit(o_n[i, j, ])
        ses[i,j] <- (obs - mean(nullvals))/sd(nullvals)
        ses_lower[i,j] <- (quantile(nullvals, probs=qs[1]) - mean(nullvals))/sd(nullvals)
        ses_upper[i,j] <- (quantile(nullvals, probs=qs[2]) - mean(nullvals))/sd(nullvals)
        raw_lower[i,j] <- quantile(nullvals, probs=qs[1])
        raw_upper[i,j] <- quantile(nullvals, probs=qs[2])
      }
    }
  }
  dimnames(ses) <- dimnames(ses_lower) <- dimnames(ses_upper) <- dimnames(raw_lower) <- dimnames(raw_upper) <- dimnames(o)
  return(list(ses=ses, ses_lower=ses_lower, ses_upper=ses_upper, raw_lower=raw_lower, raw_upper=raw_upper))
}
