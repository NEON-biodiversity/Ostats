#' Overlap of Two Circular Distributions
#'
#' This function converts input vectors to circular objects, generates kernel
#' density estimates to calculate the area of overlap between the two
#' distributions.
#' 
#' @param a a vector dataset with nrows= n individuals, with a column containing 
#'   one measurement of a certain trait.
#' @param b another matrix dataset with the same trait measurements to be compared 
#'   against a.
#' @param circular_units units of the measures.
#' @param circular_template how the data should be plotted. This set modulo, zero 
#'   and rotation of the function \code{circular} to some suitable values. 
#' @param normal If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param N the number of equally spaced points at which the density is to be
#'   estimated.
#'
#'
#' @details Circular conversion is carried out by using the function \code{circular}.
#' Kernel density estimates are then generated.If n = NULL, the default value of 512 
#' is used.Intersection density function is then calculated by taking the integral of 
#' the minimum of the two functions, from which the overlap outputs are calculated.The 
#' user must specify the bandwidth for the kernel density estimates as well as options 
#' for the circular conversion.
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
#' @seealso \code{\link{circular_overlap_24hour}} to calculate overlap for discrete hourly 
#' data.
#' 
#' @examples 
#' # circular overlap of randomly generated circular data
#' 
#' x <- rvonmises(n=100, mu=circular(pi), kappa=2)
#' y <- rvonmises(n=100, mu=circular(pi/2), kappa=2)
#' circular_overlap(x,y,circular_units = "radians", circular_template = "none",bw = 1)
#'
#' @export
circular_overlap <- function(a, b, circular_units, circular_template, normal = TRUE, bw, n = NULL) {

  # clean input
  a <- as.numeric(na.omit(a))
  b <- as.numeric(na.omit(b))

  # convert input to circular
  acirc <- circular::circular(a, units = circular_units, template = circular_template)
  bcirc <- circular::circular(a, units = circular_units, template = circular_template)

  # generate kernel densities
  # add option to use user-defined n
  # Must specify bandwidth
  if (is.null(n)) n <- 512 # Default value if not given
  da <- circular::density.circular(a, bw=bw, n=n)
  db <- circular::density.circular(b, bw=bw, n=n)
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
#' @param normal If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
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
#' # waiting for datasets
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

#' Community-level weighted median overlap using circular distributions
#'
#' Defaults to use 24-hour clock as units.
#' Uses manual calculation of density.
#'
#' @export
community_overlap_circular <- function(traits, sp, norm = TRUE, randomize_weights = FALSE) {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)

  if (nspp < 2) return(NA)

  overlaps <- numeric(0)
  abund_pairs <- numeric(0)

  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      o <- circular_overlap_24hour(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm)
      overlaps <- c(overlaps, o[1])
      harmonic_mean <- 2/(1/abunds[sp_a] + 1/abunds[sp_b])
      abund_pairs <- c(abund_pairs, harmonic_mean)
    }
  }

  if (randomize_weights) abund_pairs <- sample(abund_pairs)

  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)

}
