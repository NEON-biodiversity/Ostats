#' Overlap of Two Circular Distributions
#'
#' This function converts input vectors to circular objects, generates kernel
#' density estimates to calculate the area of overlap between the two
#' distributions.
#'
#' @param a a vector dataset with nrows= n individuals, with a column containing
#'   one measurement of a certain trait.
#' @param b another vector dataset with the same trait measurements to be compared
#'   against a.
#' @param circular_units units of the measures. Passed to \code{\link[circular]{circular}}.
#' Defaults to \code{"radians"}.
#' @param normal if TRUE, the area under all density functions is normalized to 1,
#' if FALSE, the area under all density functions is proportional to the number of
#' observations in that group.
#' @param circular_args list of additional arguments to be passed to
#' \code{\link[circular]{circular}}.
#' @param density_args list of additional arguments to be passed to
#' \code{\link[circular]{density.circular}}.
#'
#'
#' @details Circular conversion is carried out by using the function \code{\link[circular]{circular}}.
#' Kernel density estimates are then generated. Additional options for kernel density
#' estimation should be passed to \code{density_args}; otherwise, default values are used.
#' Intersection density function is then calculated by taking the integral of
#' the minimum of the two functions, from which the overlap outputs are calculated. The
#' user must specify the bandwidth for the kernel density estimates as well as options
#' for the circular conversion.
#'
#' @return The function returns a vector of three values:
#' \item{overlap_average}{the average overlap of a and b, calculated by the overlap
#' area *2 divided by the sum of areas under the two functions.}
#' \item{overlap_a}{the proportion of a that overlaps with b, calculated by the overlap
#' area divided by area under the function generated from a.}
#' \item{overlap_b}{the proportion of b that overlaps with a, calculated by the overlap
#' area divided by area under the function generated from b.}
#'
#' @seealso \code{\link{pairwise_overlap}} to calculate linear overlap between two empirical
#' density estimates.
#' @seealso \code{\link{circular_overlap_24hour}} to calculate overlap for discrete hourly
#' data.
#'
#' @examples
#' # circular overlap of two random uniform vectors of radian angles
#'
#' x <- runif(n = 100, min = 0, max = pi*3/2)
#' y <- runif(n = 100, min = pi, max = 2*pi)
#' circular_overlap(x, y, circular_units = "radians", density_args = list(bw = 1))
#'
#' @export
circular_overlap <- function(a, b, circular_units = 'radians', normal = TRUE, circular_args = list(), density_args = list()) {

  # clean input
  a <- as.numeric(na.omit(a))
  b <- as.numeric(na.omit(b))

  # convert input to circular
  acirc <- do.call(circular::circular, c(list(x = a, units = circular_units), circular_args))
  bcirc <- do.call(circular::circular, c(list(x = b, units = circular_units), circular_args))

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

  if (is.null(n)) n <- 512 # Default value if not given
  da <- circular::density.circular(acirc, bw=bw, n=n)
  db <- circular::density.circular(bcirc, bw=bw, n=n)
  d <- data.frame(x=da$x, a=da$y, b=db$y)

  # If not normalized, multiply each density entry by the length of each vector
  if (!normal) {
    d$a <- d$a * length(a)
    d$b <- d$b * length(b)
  }

  # calculate intersection densities
  d$w <- pmin(d$a, d$b)

  # integrate areas under curves
  integral_a <- sfsmisc::integrate.xy(d$x, d$a)
  integral_b <- sfsmisc::integrate.xy(d$x, d$b)
  total <- integral_a + integral_b
  intersection <- sfsmisc::integrate.xy(d$x, d$w)

  # compute overlap coefficient
  overlap_average <- 2 * intersection / total
  overlap_a <- intersection / integral_a
  overlap_b <- intersection / integral_b

  return(c(overlap_average = overlap_average, overlap_a = overlap_a, overlap_b = overlap_b))

}

#' Overlap of Two Hourly Circular Distributions
#'
#' This function calculates circular overlap based on discrete hourly observations from
#' 0 to 23.
#'
#' @param a a vector dataset with nrows= n individuals, with a column containing
#'   one measurement of a certain trait.
#' @param b another matrix dataset with the same trait measurements to be compared
#'   against a.
#' @param normal if TRUE, the area under all density functions is normalized to 1, if FALSE,
#'  the area under all density functions is proportional to the number of observations in
#'  that group.
#'
#' @details This function works for discrete data collected on the hour.Data is converted
#' to a factor with levels(0:23). It then by manually calculates the density by taking
#' the proportion of each hour. Intersection density function is then calculated by
#' taking the integral of the minimum of the two functions, from which the overlap
#' outputs are calculated.
#'
#' @return The funtion returns a vector of three values:
#' \item{overlap_average}{the average overlap of a and b, calculated by the overlap
#' area *2 divided by the sum of areas under the two functions.}
#' \item{overlap_a}{the proportion of a that overlaps with b, calculated by the overlap
#' area divided by area under the function generated from a.}
#' \item{overlap_b}{the proportion of b that overlaps with a, calculated by the overlap
#' area diveided by area under the function generated from b.}
#'
#' @seealso \code{\link{pairwise_overlap}} to calculate linear overlap between two empirical
#' density estimates.
#' @seealso \code{\link{circular_overlap}} to calculate continous circular overlap between
#' two empirical density estimates.
#'
#' @examples
#' library(dplyr)
#' a <- filter(ant_data, chamber == 1)
#' b <- filter(ant_data, chamber == 2)
#' circular_overlap_24hour(a$time,b$time)
#'
#' @export
circular_overlap_24hour <- function(a, b, normal = TRUE) {
  calc_weight <- function(x) { # a vector of hours
    tab <- table(factor(x,  levels=as.character(0:23)),
                 useNA="ifany")

    dimnames(tab) <- NULL
    if (normal) {
      weights <- tab / sum(tab)
    } else {
      weights <- tab
    }
    mat <- cbind( weights=weights, points=0:23 )
    mat
  }

  A <- calc_weight(a)
  B <- calc_weight(b)

  d <- data.frame(x = A[,'points'], a = A[,'weights'], b = B[,'weights'])

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

}
