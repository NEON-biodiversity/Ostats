### TEST VERSION: MULTIVARIATE
pairwise_overlap_mv <- function(a, b, normal = TRUE, density_args = list()) {

  # clean input
  a <- a[complete.cases(a), ]
  b <- b[complete.cases(b), ]

  # Scale input if normalized
  if (normal) {
    a <- scale(a)
    b <- scale(b)
  }

  # FIXME allow user to pass arguments here
  # Convert each of the input matrices to hypervolume.
  hv_a <- hypervolume::hypervolume(a, method = 'gaussian', verbose = FALSE)
  hv_b <- hypervolume::hypervolume(b, method = 'gaussian', verbose = FALSE)

  # FIXME allow user to pass arguments here
  # Calculate hypervolume set operations
  hv_set_ab <- hypervolume::hypervolume_set(hv_a, hv_b, num.points.max = NULL, verbose = FALSE, check.memory = FALSE, distance.factor = 1)

  # Calculate hypervolume overlap statistic
  hv_overlap_ab <- hypervolume::hypervolume_overlap_statistics(hv_set_ab)

  # Return Sorensen similarity (equivalent to univariate overlap statistic)
  # Return other values to correspond with our univariate version
  out <- c(hv_overlap_ab['sorensen'], 1 - hv_overlap_ab['frac_unique_1'], 1 - hv_overlap_ab['frac_unique_2'])
  names(out) <- c('overlap_average', 'overlap_a', 'overlap_b')
  return(out)
}
