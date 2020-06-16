#' Overlap between Two Empirical Density Estimates
#'
#' This function generates kernel density estimates for two datasets on a common
#' grid to calculate the area of overlap between the two estimates.
#' See \url{http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities},
#' this function is derived from code posted in the answer by user mmk.
#'
#'
#' @param a a matrix dataset with nrows= n individuals, with a column containing 
#'   one measurement of a certain trait.
#' @param b another matrix dataset with the same trait measurements to be compared 
#'   against a.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#'
#' @details This function generates kernel density estimates using the density() 
#' function for two datasets on a common grid.If bw = NULL, the default 'nrd0'is 
#' used. If n = NULL, the default value of 512 is used.Intersection density function 
#' is then calculated by taking the integral of the minimum of the two functions, 
#' from which the coefficients are calculated.
#'
#' @return The funtion returns a vector of three values: 
#' \item{overlap}{the overlap area *2 divided by the sum of areas under the two 
#' functions.} 
#' \item{overlap_a}{the overlap area divided by area under the function generated from 
#' a.}
#' \item{overlap_b}{the overlap area diveided by area under the function generated 
#' from b.}
#'
#' @references http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities
#'   
#' @examples
#' #overlap of miles per gallon between 4-cylinder and 6-cylinder cars
#' library(ggplot2)
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
#'
#' 
# Overlap between two empirical density estimates

pairwise_overlap <- function(a, b, normal = TRUE, bw = NULL, N = NULL) {
  
  # clean input
  a <- as.numeric(na.omit(a))
  b <- as.numeric(na.omit(b))
  
  # define limits of a common grid, adding a buffer so that tails aren't cut off
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  # add option to use user-defined bandwidth and n
  if (is.null(bw)) bw <- 'nrd0' # Defaults to traditional method if not given
  if (is.null(N)) N <- 512 # Default value if not given
  da <- density(a, from=lower, to=upper, bw=bw, n=N)
  db <- density(b, from=lower, to=upper, bw=bw, n=N)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # If not normalized, multiply each density entry by the length of each vector
  if (normal!=TRUE) {
    d$a <- d$a * length(a)
    d$b <- d$b * length(b)
  }
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  suppressMessages(require(sfsmisc))
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  overlap <- 2 * intersection / total
  overlap_a <- intersection / integrate.xy(d$x, d$a)
  overlap_b <- intersection / integrate.xy(d$x, d$b)
  
  return(c(overlap = overlap, overlap_a = overlap_a, overlap_b = overlap_b))
  
}
