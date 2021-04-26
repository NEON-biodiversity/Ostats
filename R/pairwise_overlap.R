#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#' FIXME new documentation could be added here, if we want to still export this function.
#'
#' @noRd
pairwise_overlap <- function(a, b, density_args = list(), hypervolume_set_args = list()) {

  # Check structure of inputs a and b.
  # If they are unidimensional use density(), if >1 dimension use hypervolume()
  # If the dimensions don't match (number of columns in a != number of columns in b), return error.
  if (!inherits(a, "Hypervolume")) {
    # Univariate case
    # calculate intersection density
    w <- pmin(a$y, b$y)

    # integrate areas under curves
    total <- sfsmisc::integrate.xy(a$x, a$y) + sfsmisc::integrate.xy(b$x, b$y)
    intersection <- sfsmisc::integrate.xy(a$x, w)

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

