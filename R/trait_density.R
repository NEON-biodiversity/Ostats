#' Univariate density functions (linear and circular) or multivariate hypervolumes
#'
#' This function accepts a numeric vector in the univariate case or a numeric
#' matrix in the multivariate case and calculates the kernel density function
#' or the hypervolume, respectively. It can also accept circular data in the
#' univariate case, whether continuous or discrete.
#'
#' @param x Numeric vector or matrix.
#' @param grid_limits 2 x n numeric matrix. Each column contains the minimum and
#'   maximum value over which each trait's density function will be estimated.
#' @param normal Passed from \code{\link{Ostats}}.
#' @param discrete Passed from \code{\link{Ostats}}.
#' @param circular Passed from \code{\link{Ostats}}.
#' @param unique_values Vector of unique possible values for \code{x}.
#'   Only used for discrete data types.
#' @param density_args Passed from \code{\link{Ostats}}.
#' @param circular_args Passed from \code{\link{Ostats}}.
#'
#' @details This is an internal function not intended to be called directly.
#'
#' @return In the univariate case, returns a data frame with columns \code{x}
#' and \code{y}, where \code{x} is an evenly spaced sequence of values between
#' the minimum and maximum limit, and \code{y} is the kernel density function
#' evaluated at \code{x}.
#'
#' In the multivariate case, returns an object of class \code{"hypervolume"}.
#'
#' @noRd
trait_density <- function(x, grid_limits, normal, discrete, circular, unique_values, density_args, circular_args) {
  if (is.vector(x)) {
    # Univariate case
    if (!circular & !discrete) {
      # clean input
      x <- as.numeric(stats::na.omit(x))

      # generate kernel densities
      # Bandwidth method defaults to nrd0 if not given, n defaults to 512 if not given
      if ('bw' %in% names(density_args)) {
        bw <- density_args[['bw']]
      } else {
        bw <- 'nrd0'
      }
      if ('n' %in% names(density_args)) {
        n <- density_args[['n']]
      } else {
        n <- 512
      }

      dens <- stats::density(x, from=grid_limits[1, 1], to=grid_limits[2, 1], bw=bw, n=n)

      d <- data.frame(x=dens$x, y=dens$y)

      # If not normalized, multiply each density entry by the length of each vector
      if (!normal) {
        d$y <- d$y * length(x)
      }

      return(d)

    }

    if (circular & !discrete) {
      # clean input
      x <- as.numeric(stats::na.omit(x))

      # convert input to circular
      xcirc <- do.call(circular::circular, c(list(x = x, units = 'radians'), circular_args))

      # generate kernel densities
      # Bandwidth defaults to 1 if not given, n defaults to 512 if not given
      if ('bw' %in% names(density_args)) {
        bw <- density_args[['bw']]
      } else {
        bw <- 1
      }
      if ('n' %in% names(density_args)) {
        n <- density_args[['n']]
      } else {
        n <- 512
      }

      dens <- circular::density.circular(xcirc, bw=bw, n=n)

      d <- data.frame(x=dens$x, y=dens$y)

      # If not normalized, multiply each density entry by the length of each vector
      if (!normal) {
        d$y <- d$y * length(x)
      }

      return(d)

    }

    if (discrete) {
      x_weights <- calc_weight(x, normal, unique_values)

      d <- data.frame(x = x_weights[,'points'], y = x_weights[,'weights'])

      return(d)
    }

  } else {
    # Multivariate case
    # Use default method argument to hypervolume::hypervolume if none are provided
    if (!'method' %in% names(density_args)) {
      density_args[['method']] <- 'gaussian'
    }
    # Also set to verbose = FALSE if that argument is not provided.
    if (!'verbose' %in% names(density_args)) {
      density_args[['verbose']] <- FALSE
    }

    # clean input
    x <- x[stats::complete.cases(x), ]

    # Scale input if normalized
    if (normal) {
      x <- scale(x)
    }

    # Find hypervolume.
    invisible(utils::capture.output(suppressWarnings(suppressMessages({
      hv <- do.call(hypervolume::hypervolume, args = c(list(data = x), density_args))
    }))))

    return(hv)
  }
}

#' Function to calculate weights for discrete data.
#' @noRd
calc_weight <- function(x, normal, x_values) {
  tab <- table(factor(x, levels=as.character(x_values)),
               useNA="ifany")

  dimnames(tab) <- NULL
  if (normal) {
    weights <- tab / sum(tab)
  } else {
    weights <- tab
  }
  mat <- cbind(weights = weights, points = x_values)
  return(mat)
}
