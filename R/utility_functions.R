# Utility functions, not exported
# -------------------------------

# Calculate standardized effect sizes of overlap statistics against null distributions
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
