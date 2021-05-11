#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' For the univariate case,
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#' @param a a numeric vector or matrix. Overlap is calculated between a and b.
#'   If a and b are vectors, \code{\link[stats]{density}} is used to calculate unidimensional overlap.
#'   If a and b are matrices, with each column representing a trait or dimension,
#'   \code{\link[hypervolume]{hypervolume}} is used to calculate multidimensional overlap.
#' @param b a numeric vector or matrix. The number of columns of a and b must be equal.
#' @param density_args Additional arguments to pass to \code{\link[stats]{density}}
#'   in the univariate case, or \code{\link[hypervolume]{hypervolume}} in
#'   the multivariate case.
#'   See \code{\link{Ostats}} and \code{\link{Ostats_multivariate}}.
#' @param hypervolume_set_args Additional arguments to pass to
#'   \code{\link[hypervolume]{hypervolume_set}}. See
#'   \code{\link{Ostats_multivariate}}.
#'
#' @details This is an internal function not intended to be called directly.
#'
#' @return A single numeric value that may range between 0 and 1.
#'
#' @noRd
pairwise_overlap <- function(a, b, discrete, density_args = list(), hypervolume_set_args = list()) {

  # Check structure of inputs a and b.
  # If they are unidimensional use density(), if >1 dimension use hypervolume()
  # If the dimensions don't match (number of columns in a != number of columns in b), return error.
  if (!inherits(a, "Hypervolume")) {
    # Univariate case
    # calculate intersection density
    w <- pmin(a$y, b$y)

    if (!discrete) {
      # If continuous, integrate the areas under curves
      total <- sfsmisc::integrate.xy(a$x, a$y) + sfsmisc::integrate.xy(b$x, b$y)
      intersection <- sfsmisc::integrate.xy(a$x, w)
    } else {
      # If discrete, take the overlaps at each discrete point
      total <- sum(a$y + b$y)
      intersection <- sum(w)
    }

      # compute overlap coefficient (Sorensen)
      overlap_average <- 2 * intersection / total

      return(overlap_average)

  } else {
    # Multivariate case

    # Set to verbose = FALSE if that argument is not provided.
    if (!'verbose' %in% names(density_args)) {
      density_args[['verbose']] <- FALSE
    }

    # If num.points.max and distance.factor are not provided to pass to hypervolume_set, use default.
    if (!'num.points.max' %in% names(hypervolume_set_args)) {
      hypervolume_set_args[['num.points.max']] <- NULL
    }

    if (!'distance.factor' %in% names(hypervolume_set_args)) {
      hypervolume_set_args[['distance.factor']] <- 1
    }

    # Suppress all progress messages from hypervolume functions, including those from underlying C functions.
    invisible(utils::capture.output(suppressWarnings(suppressMessages({

      # Calculate hypervolume set operations
      hv_set_ab <- do.call(hypervolume::hypervolume_set,
                           c(list(hv1 = a,
                                  hv2 = b,
                                  verbose = density_args[['verbose']],
                                  check.memory = FALSE
                           ),
                           hypervolume_set_args)
      )

      # Calculate hypervolume overlap statistic
      hv_overlap_ab <- hypervolume::hypervolume_overlap_statistics(hv_set_ab)

    }))))

    # Return Sorensen similarity (equivalent to univariate overlap statistic)
    return(hv_overlap_ab['sorensen'])

  }

}

