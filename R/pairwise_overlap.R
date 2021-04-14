#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#' FIXME new documentation could be added here, if we want to still export this function.
#'
#' @export
pairwise_overlap <- function(a, b, density_args = list()) {

  # Check structure of inputs a and b.
  # If they are unidimensional use density(), if >1 dimension use hypervolume()
  # If the dimensions don't match (number of columns in a != number of columns in b), return error.
  if (!inherits(a, "Hypervolume")) {
    # Univariate case
    # calculate intersection densities
    w <- pmin(a$y, b$y)

    # integrate areas under curves
    total <- sfsmisc::integrate.xy(a$x, a$y) + sfsmisc::integrate.xy(b$x, b$y)
    intersection <- sfsmisc::integrate.xy(a$x, w)

    # compute overlap coefficient
    overlap_average <- 2 * intersection / total

    return(overlap_average)

  } else {

    # Set to verbose = FALSE if that argument is not provided.
    if (!'verbose' %in% names(density_args)) {
      density_args[['verbose']] <- FALSE
    }

    # Suppress all progress messages from hypervolume functions, including those from underlying C functions.
    invisible(utils::capture.output(suppressWarnings(suppressMessages({

      # Calculate hypervolume set operations
      # This uses default arguments except for verbose.
      # Later we may implement the ability to modify these arguments.
      hv_set_ab <- hypervolume::hypervolume_set(a, b, num.points.max = NULL, verbose = density_args[['verbose']], check.memory = FALSE, distance.factor = 1)

      # Calculate hypervolume overlap statistic
      hv_overlap_ab <- hypervolume::hypervolume_overlap_statistics(hv_set_ab)

    }))))

    # Return Sorensen similarity (equivalent to univariate overlap statistic)
    return(hv_overlap_ab['sorensen'])

  }

}

