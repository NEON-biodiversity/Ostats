#'@title  Calculation of standardized effect sizes of overlap statistics against
#'  null distributions
#'
#'@description Extract quantiles to get standardized effect sizes for the
#'  overlap stats
#'
#'@details Implementation of null model that swaps means of distributions to
#'  test whether distributions are more evenly spaced than expected by chance -
#'  it is part of the Ostats function
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@param o an matrix data containing the overlapping values of normal
#'  distribution of the observed data (the observed local O-stats for each
#'  community)
#'@param o_n an array object containing the overlapping values of normal
#'  distribution of the randomize data - (the null distributions of O-stats
#'  function)
#'@param qs quantiles of the null distribution to extract, usually 0.025 and
#'  0.975
#'
#'@return \item{ses}{a matrix containing values of standardized effect sizes for
#'observed data of each of #' the locals (rows) and traits(columns) avaliated.
#'To see if the overlap values can be expected by #' chance see returns
#'ses_upper and ses_lower. If the value corresponding to the local in ses is
#'between ses_upper and ses_lower, the overlapping between traits in community
#'can be expected by chance.} \item{ses_lower}{a matrix containing null values
#'extracted for the quantil of qs[1] (usually 0.025)} \item{ses_upper}{a matrix
#'containing null values extracted for the quantil of qs[2] (usually 0.975)}
#'\item{raw_lower}{a matrix containing null values extracted for the quantil of
#'qs[1] (usually 0.025) without subtracting the mean and divided by standard
#'deviation (see equation for calculating the standardize effect size on the
#'function code)} \item{raw_upper}{a matrix containing null values extracted for
#'the quantil of qs[2] (usually 0.975) without subtracting the mean and divided
#'by standard deviation (see equation for calculating the standardize effect
#'size on the function code)}
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
