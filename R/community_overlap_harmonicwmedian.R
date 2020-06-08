#' Median Pairwise Overlap of Species Trait Distributions in a Community
#'
#' This function calculates the median of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community, weighted 
#' by the harmonic mean of the abundances of each species pair.
#'
#' @param traits a matrix dataset with nrows = n individuals, ncols = n traits.
#' @param sp a factor with length equal to nrow(traits) that indicates the taxon
#'   of each individual.
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
#' @details The funtion evaluates pairwise overlaps of density estimates of all 
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
#' 

# Median pairwise overlap of trait distributions of all species in a community
# The median is weighted by the harmonic mean of the abundances of each species pair for which overlap is calculated.
community_overlap_harmonicwmedian <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL, randomize_weights = FALSE) {
  require(matrixStats)
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
  
  weightedMedian(x = overlaps, w = abund_pairs)
  
}
