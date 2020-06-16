#' Community Overlap Calculation
#'
#' This function calculates the median or mean of pairwise overlaps between density 
#' estimates of trait distributions of all species within a community, which can 
#' be weighted by species abundances. It can also evaluate pairwise distances 
#' calculated from trait means.
#' 
#' @param traits a vector of trait measurement.
#' @param sp a vector with length equal to length(traits) that indicates the
#'taxon of each individual.
#' @param norm If TRUE, assume data are normally distributed; if FALSE,
#'   additional normalization step is carried out by multiplying each density 
#'   entry by the length of each vector.
#' @param output specifies whether median or mean is calculated.
#' @param weights specifies weights to be used to calculate the median or mean.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#'   that this is the standard deviation of the smoothing kernel.
#' @param n the number of equally spaced points at which the density is to be
#'   estimated.
#' @param itv If TRUE, pairwise overlaps are evaluated;if FALSE, pairwise distances
#' are calculated, which means each species is reduced to a trait mean and the 
#' intraspecific trait variation (ITV) is ignored.
#' @param randomize_weights If TRUE, randomize weights given to pairwise overlaps
#'   within a community. This can be used to generate null models.
#' @param display If TRUE, display all pairwise overlaps.
#'
#' @details The funtion evaluates weighted pairwise distances or overlaps of density 
#' estimates of all species in a community taking complete cases with species abundances 
#' greater than 1 from the dataset. The default calculates the median of pairwise overlaps
#' for the whole community using the harmonic means of abundances of the species pairs as 
#' weights, which minimizes the effect of outliners.If the argument weights == NULL, no 
#' weights are used for the calculation of mean/median. If weights == "mean", means of 
#' abundances are used as weights. To change the output to mean, specify the argument 
#' output == "mean".
#' 
#' @return At default, the function returns a median of pairwise overlaps weighted by 
#' harmonic means of abundances for the community. 
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
#' # Calculate median of pairwise overlaps for the community,weighted by harmonic median
#' of abundances
#' community_overlap(traits = as.matrix(dat$log_weight), 
#' sp = factor(dat$taxonID))
#' @export
#'
community_overlap_merged <- function(traits, sp, norm = TRUE, output = "median", weights = "hmean", bw = NULL, n = NULL, itv = TRUE, randomize_weights = FALSE, display = FALSE) {
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
      if (itv==TRUE) {
        o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      }
      else(o <- abs(mean(traitlist[[sp_a]], na.rm=T) - mean(traitlist[[sp_b]], na.rm=T)))
        
      overlaps <- c(overlaps, o[1])
        if (weights == "hmean")
          harmonic_mean <- 2/(1/abunds[sp_a] + 1/abunds[sp_b])
          abund_pairs <- c(abund_pairs, harmonic_mean)
        if (weights == "mean")
          abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))
      
    }
  }
  if (display==TRUE) return(overlaps)
  
  if (randomize_weights) abund_pairs <- sample(abund_pairs)
  
  if (output == "median" & weights != NULL)
    matrixStats::weightedMedian(x = overlaps, w = abund_pairs)
  if (output == "median" & weights == NULL) median(overlaps)
  if (output == "mean" & weights != NULL)
    weighted.mean(x = overlaps, w = abund_pairs)
  if (output == "mean" & weights == NULL) mean(overlaps)
  
}

