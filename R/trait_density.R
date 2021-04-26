#' Univariate density functions (linear and circular) or multivariate hypervolumes
#' @noRd
trait_density <- function(x, grid_limits, normal, data_type, unique_values, density_args, circular_args) {
  if (is.vector(x)) {
    # Univariate case
    if (data_type %in% 'linear') {
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

    if (data_type %in% 'circular') {
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

    if (data_type %in% 'circular_discrete') {
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

#' Function to calculate hourly weights
#' @noRd
calc_weight <- function(x, normal, x_values) {
  tab <- table(factor(x,  levels=as.character(x_values)),
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
