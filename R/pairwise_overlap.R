#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#'
#' @param a a vector dataset with nrows= n individuals, with a column containing
#'   one measurement of a certain trait.
#' @param b another matrix dataset with the same trait measurements to be compared
#'   against a.
#' @param normal If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density
#'   entry by the length of each vector.
#' @param density_args list of additional arguments to be passed to
#'   \code{\link[stats]{density}}.
#'
#' @details This function generates kernel density estimates using the
#' \code{\link[stats]{density}}  function for two datasets on a common grid.
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
#' area diveided by area under the function generated from b.}
#'
#' @references http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities
#'
#' @seealso \code{\link{circular_overlap}} to calculate circular overlap.
#'
#' @examples
#' #overlap of miles per gallon between 4-cylinder and 6-cylinder cars
#' library(dplyr)
#' a <- mtcars %>%
#'   filter(cyl==6)
#' a <- a[,1]
#' b <- mtcars %>%
#'   filter(cyl==4)
#' b <- b[,1]
#' pairwise_overlap(a,b)
#'
#' @export
pairwise_overlap <- function(a, b, normal = TRUE, density_args = list()) {

  # clean input
  a <- as.numeric(na.omit(a))
  b <- as.numeric(na.omit(b))

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

  da <- density(a, from=lower, to=upper, bw=bw, n=n)
  db <- density(b, from=lower, to=upper, bw=bw, n=n)
  d <- data.frame(x=da$x, a=da$y, b=db$y)

  # If not normalized, multiply each density entry by the length of each vector
  if (normal!=TRUE) {
    d$a <- d$a * length(a)
    d$b <- d$b * length(b)
  }

  # calculate intersection densities
  d$w <- pmin(d$a, d$b)

  # integrate areas under curves
  total <- sfsmisc::integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- sfsmisc::integrate.xy(d$x, d$w)

  # compute overlap coefficient
  overlap_average <- 2 * intersection / total
  overlap_a <- intersection / sfsmisc::integrate.xy(d$x, d$a)
  overlap_b <- intersection / sfsmisc::integrate.xy(d$x, d$b)

  return(c(overlap_average = overlap_average, overlap_a = overlap_a, overlap_b = overlap_b))

}
