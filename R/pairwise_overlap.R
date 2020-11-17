#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#'
#' @param a a numeric vector or matrix. Overlap is calculated between a and b.
#'   If a and b are vectors, \code{\link[stats]{density}} is used to calculate unidimensional overlap.
#'   If a and b are matrices, with each column representing a trait or dimension,
#' \code{\link[hypervolume]{hypervolume}} is used to calculate multidimensional overlap.
#' @param b a numeric vector or matrix. The number of columns of a and b must be equal.
#' @param normal if TRUE, the area under all density functions is normalized to 1,
#'  if FALSE, the area under all density functions is proportional to the number of
#'  observations in that group.
#' @param density_args list of additional arguments to be passed to
#'   \code{\link[stats]{density}} in the univariate case, or
#'   \code{\link[hypervolume]{hypervolume}} in the multivariate case.
#'
#' @details This function generates kernel density estimates using the
#' \code{\link[stats]{density}} function for two datasets on a common grid.
#' Default values for \code{bw} and \code{n} are used if not provided in \code{density_args}.
#' Intersection density function
#' is then calculated by taking the integral of the minimum of the two functions,
#' from which the overlap outputs are calculated.
#'
#' @return The funtion returns a vector of three values:
#' \item{overlap_average}{the average overlap of a and b, calculated by the overlap
#' area *2 divided by the sum of areas under the two functions.}
#' \item{overlap_a}{the proportion of a that overlaps with b, calculated by the overlap
#' area divided by area under the function generated from a.}
#' \item{overlap_b}{the proportion of b that overlaps with a, calculated by the overlap
#' area divided by area under the function generated from b.}
#'
#' @references http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities
#'
#' @seealso \code{\link{circular_overlap}} to calculate circular overlap.
#'
#' @examples
#' #overlap of miles per gallon between 4-cylinder and 6-cylinder cars
#' a <- mtcars[mtcars$cyl == 6, ]
#' a <- a[,1]
#' b <- mtcars[mtcars$cyl == 4, ]
#' b <- b[,1]
#' pairwise_overlap(a,b)
#'
#' @export
pairwise_overlap <- function(a, b, normal = TRUE, density_args = list()) {

  # Check structure of inputs a and b.
  # If they are unidimensional use density(), if >1 dimension use hypervolume()
  # If the dimensions don't match (number of columns in a != number of columns in b), return error.
  if (is.vector(a) & is.vector(b)) {
    # Univariate case

    # clean input
    a <- as.numeric(stats::na.omit(a))
    b <- as.numeric(stats::na.omit(b))

    # define limits of a common grid, adding a buffer so that tails aren't cut off
    lower <- min(c(a, b)) - 1
    upper <- max(c(a, b)) + 1

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

    da <- stats::density(a, from=lower, to=upper, bw=bw, n=n)
    db <- stats::density(b, from=lower, to=upper, bw=bw, n=n)
    d <- data.frame(x=da$x, a=da$y, b=db$y)

    # If not normalized, multiply each density entry by the length of each vector
    if (!normal) {
      d$a <- d$a * length(a)
      d$b <- d$b * length(b)
    }

    # calculate intersection densities
    d$w <- pmin(d$a, d$b)

    # integrate areas under curves
    total <- sfsmisc::integrate.xy(d$x, d$a) + sfsmisc::integrate.xy(d$x, d$b)
    intersection <- sfsmisc::integrate.xy(d$x, d$w)

    # compute overlap coefficient
    overlap_average <- 2 * intersection / total
    overlap_a <- intersection / sfsmisc::integrate.xy(d$x, d$a)
    overlap_b <- intersection / sfsmisc::integrate.xy(d$x, d$b)

    return(c(overlap_average = overlap_average, overlap_a = overlap_a, overlap_b = overlap_b))

  } else {
    if (dim(a)[2] != dim(b)[2]) {
      stop('Number of columns across species do not match.')
    }

    # Use default method argument to hypervolume::hypervolume if none are provided
    if (!'method' %in% names(density_args)) {
      density_args[['method']] <- 'gaussian'
    }
    # Also set to verbose = FALSE if that argument is not provided.
    if (!'verbose' %in% names(density_args)) {
      density_args[['verbose']] <- FALSE
    }

    # clean input
    a <- a[stats::complete.cases(a), ]
    b <- b[stats::complete.cases(b), ]

    # Scale input if normalized
    if (normal) {
      a <- scale(a)
      b <- scale(b)
    }

    # Suppress all progress messages from hypervolume functions, including those from underlying C functions.
    invisible(utils::capture.output(suppressWarnings(suppressMessages({

      # Convert each of the input matrices to hypervolume.
      # User-input arguments are passed using do.call. This may be error prone so could be improved later.
      hv_a <- do.call(hypervolume::hypervolume, args = c(list(data = a), density_args))
      hv_b <- do.call(hypervolume::hypervolume, args = c(list(data = b), density_args))

      # Calculate hypervolume set operations
      # This uses default arguments except for verbose. I think it is not a good idea to allow these arguments to be changed.
      hv_set_ab <- hypervolume::hypervolume_set(hv_a, hv_b, num.points.max = NULL, verbose = density_args[['verbose']], check.memory = FALSE, distance.factor = 1)

      # Calculate hypervolume overlap statistic
      hv_overlap_ab <- hypervolume::hypervolume_overlap_statistics(hv_set_ab)

    }))))

    # Return Sorensen similarity (equivalent to univariate overlap statistic)
    # Return other values to correspond with our univariate version
    out <- c(hv_overlap_ab['sorensen'], 1 - hv_overlap_ab['frac_unique_1'], 1 - hv_overlap_ab['frac_unique_2'])
    names(out) <- c('overlap_average', 'overlap_a', 'overlap_b')
    return(out)

  }

}

