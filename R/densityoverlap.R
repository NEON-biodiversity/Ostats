#' Median Pairwise Overlap across All Species in a Community
#' 
#' This function calculates the median of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'   taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' 
#' @details The funtion evaluates pairwise overlaps of density estimates of all 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset to calculate the median for the community.
#' 
#' @return The function returns a median of pairwise overlaps of all species in 
#' the community.
#' 
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#' @seealso \code{\link{community_overlap_harmonicwmedian}} to see a function that
#' uses harmonic means of abundances of species pairs as weights to calculate median.
#' 
#' @export
community_overlap <- function(traits, sp, normal = TRUE, bw = NULL, N = NULL) {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  dat <- dat[dat$sp %in% names(abunds)[abunds>1], ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)
  
  if (nspp < 2) return(NA)
  
  overlaps <- numeric(0)
  
  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], normal = normal, bw = bw, N = N)
      overlaps <- c(overlaps, o[2:3])
    }
  }
  
  median(overlaps)
  
}
#' Abundance-weighted Mean of Pairwise Overlaps across All Species in a Community
#' 
#' This function calculates the mean of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community, weighted by 
#' the mean of abundances of the species pairs.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'   taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' 
#' @details The funtion evaluates pairwise overlaps of density estimates of all 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset, using the mean of abundances of species 
#' pairs as weights to calculate the weighted mean of overlaps for the community.
#' 
#' @return The function returns a mean of pairwise overlaps of all species in 
#' the community, weighted by means of abundances of species pairs.
#' 
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#' @seealso \code{\link{community_overlap_harmonicwmedian}} to see a function that
#' uses harmonic means of abundances of species pairs as weights to calculate median 
#' of overlaps.
#' 
#' @export
community_overlap_wm <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
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
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
	  abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))
    }
  }

  weighted.mean(x = overlaps, w = abund_pairs)

}

#' Abundance-weighted Median of Pairwise Overlaps across All species in a Community
#' 
#' This function calculates the median of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community, weighted by 
#' the mean of abundances of the species pairs.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'   taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' 
#' @details The funtion evaluates pairwise overlaps of density estimates of all 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset, using the mean of abundances of species 
#' pairs as weights to calculate the weighted median of overlaps for the community.
#' 
#' @return The function returns a median of pairwise overlaps of all species in 
#' the community, weighted by means of abundances of species pairs.
#' 
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#' @seealso \code{\link{community_overlap_harmonicwmedian}} to see a function that
#' uses harmonic means of abundances of species pairs as weights to calculate median 
#' of overlaps.
#' 
#' @export
community_overlap_wmedian <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
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
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], normal = norm, bw = bw, N = n)
      overlaps <- c(overlaps, o[1])
      abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))
    }
  }

  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)

}

#' Abundance-weighted Median of Pairwise Overlap of Species Trait Distributions 
#' in a Community
#'
#' This function calculates the median of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community, weighted 
#' by the harmonic mean of the abundances of each species pair.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' @param randomize_weights If TRUE, randomize weights given to pairwise overlaps
#'   within a community. This can be used to generate null models.
#'
#' @details The funtion evaluates weighted pairwise overlaps of density estimates of all 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset. The median of pairwise overlaps is calculated 
#' for the whole community using the harmonic means of abundances of the species 
#' pairs as weights.
#'
#' @return The function returns a median of pairwise overlaps weighted by harmonic 
#' means of abundances for the community.
#'
#' @references Read, Q. D. et al. Among-species overlap in rodent body size
#'   distributions predicts species richness along a temperature gradient.
#'   Ecography 41, 1718-1727 (2018).
#'   
#' @note The function so far only supports overlap statistics for one trait.
#'
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#'
#' @examples
#' 
#' library(tidyverse)
#' library(Ostats)
#' 
#' # Load data from web archive and Keep only the relevant part of data
#' dat <- read_csv('https://ndownloader.figshare.com/files/9167548') %>%
#'    filter(siteID %in% 'HARV') %>%
#'    select(siteID, taxonID, weight) %>%
#'    filter(!is.na(weight)) %>%
#'    mutate(log_weight = log10(weight))
#' 
#' # Calculate median of pairwise overlaps for the community
#' community_overlap_harmonicwmedian(traits = as.matrix(dat$log_weight), 
#' sp = factor(dat$taxonID))
#' @export
#'
community_overlap_harmonicwmedian <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL, randomize_weights = FALSE) {
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
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
      harmonic_mean <- 2/(1/abunds[sp_a] + 1/abunds[sp_b])
      abund_pairs <- c(abund_pairs, harmonic_mean)
    }
  }

  if (randomize_weights) abund_pairs <- sample(abund_pairs)

  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)

}

#' Abundance-weighted Median of Pairwise Distance between Species Trait Means 
#' in a Community
#' 
#' This function calculates the median of pairwise distances between species 
#' trait means within a community, weighted by the mean of abundances of the 
#' species pairs. This reduces each species to a trait mean and the intraspecific 
#' trait variation (ITV) is ignored.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'   taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' 
#' @details The funtion evaluates means of trait measurements of each 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset, using the mean of abundances of species 
#' pairs as weights to calculate the weighted median of pairwise distances between
#' species for the community.
#' 
#' @return The function returns a median of pairwise distances of all species in 
#' the community, weighted by means of abundances of species pairs.
#' 
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#' @seealso \code{\link{community_overlap_harmonicwmedian}} to see a function that
#' uses harmonic means of abundances of species pairs as weights to calculate median 
#' of overlaps.
#' 
#' @export
community_overlap_noitv <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
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
      o <- abs(mean(traitlist[[sp_a]], na.rm=T) - mean(traitlist[[sp_b]], na.rm=T))
      overlaps <- c(overlaps, o[1])
      abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))
    }
  }

  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)

}

#' Display All Pairwise Overlap Values
#'
#' This function displays all pairwise overlaps between species pairs in a 
#' community.
#'
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'   taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' 
#' @details The funtion evaluates pairwise overlaps of density estimates of all 
#' species in a community taking complete cases with species abundances greater 
#' than 1 from the original dataset, and displays all pairwise overlap values
#' within the community.
#' 
#' @return The function returns all pairwise overlap values of trait distributions
#' of all species in the community.
#' 
#' @seealso \code{\link{pairwise_overlap}} to calculate overlap between two empirical 
#' density estimates.
#' @seealso \code{\link{community_overlap_harmonicwmedian}} to see a function that
#' uses harmonic means of abundances of species pairs as weights to calculate median 
#' of overlaps.
#' 
#' @export
#' 
community_overlap_distributions <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  dat <- dat[dat$sp %in% names(abunds)[abunds>1], ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)

  if (nspp < 2) return(NA)

  overlaps <- numeric(0)

  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[2:3])
    }
  }

  return(overlaps)

}

