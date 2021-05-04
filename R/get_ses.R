#'@title  Calculation of standardized effect sizes of overlap statistics against
#'  null distributions
#'
#'@description Extract quantiles to get standardized effect sizes for the
#'  overlap stats
#'
#'@details This is an internal function not intended to be called directly.
#'
#'@seealso \code{\link{Ostats}} to calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@param o a matrix containing the observed local O-stats for each
#'  community
#'@param o_n an array object containing the the null distributions of O-stats
#'@param qs quantiles of the null distribution to extract
#'
#'@return \item{ses}{a matrix containing values of standardized effect sizes for
#'each site (rows) and trait (columns) evaluated.}
#'\item{ses_lower}{a matrix containing z-transformed null values
#' extracted for the quantile of \code{qs[1]} (default 0.025)}
#'\item{ses_upper}{a matrix containing z-transformed null values
#' extracted for the quantile of \code{qs[2]} (default 0.975)}
#'\item{raw_lower}{a matrix containing raw null values
#' extracted for the quantile of \code{qs[1]} (default 0.025)}
#'\item{raw_upper}{a matrix containing null values
#' extracted for the quantile of \code{qs[2]} (default 0.975)}
#'
#'@noRd
get_ses <- function(o, o_n, qs) {
  ses <- ses_lower <- ses_upper <- raw_lower <- raw_upper <- matrix(NA, nrow=nrow(o), ncol=ncol(o))

  for (i in 1:nrow(o)) {
    for (j in 1:ncol(o)) {
      if(!is.na(o[i,j])) {
        obs <- o[i,j]
        nullvals <- stats::na.omit(o_n[i, j, ])
        ses[i,j] <- (obs - mean(nullvals))/stats::sd(nullvals)
        ses_lower[i,j] <- (stats::quantile(nullvals, probs=qs[1]) - mean(nullvals))/stats::sd(nullvals)
        ses_upper[i,j] <- (stats::quantile(nullvals, probs=qs[2]) - mean(nullvals))/stats::sd(nullvals)
        raw_lower[i,j] <- stats::quantile(nullvals, probs=qs[1])
        raw_upper[i,j] <- stats::quantile(nullvals, probs=qs[2])
      }
    }
  }
  dimnames(ses) <- dimnames(ses_lower) <- dimnames(ses_upper) <- dimnames(raw_lower) <- dimnames(raw_upper) <- dimnames(o)
  return(list(ses=ses, ses_lower=ses_lower, ses_upper=ses_upper, raw_lower=raw_lower, raw_upper=raw_upper))
}
