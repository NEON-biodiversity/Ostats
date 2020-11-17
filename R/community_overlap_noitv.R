

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
#' @param output specifies whether median or mean is calculated.
#' @param weight_type specifies weights to be used to calculate the median or mean.
#'
#'
#' @details The function evaluates the means of trait measurements of each
#' species in a community taking complete cases with species abundances greater
#' than 1 from the original dataset, and calculates weighted mean or median of pairwise
#' distances between all species in a community. The default calculates the median
#' using the harmonic means of abundances of the species pairs as weights, which
#' minimizes the effect of outliners and rare species.If the argument
#' weight_type == "none", no weights are used for the calculation of mean/median. If
#' weight_type == "mean", arithmetic means of abundances are used as weights. To change the
#' output to mean, specify the argument output == "mean".
#'
#' @return The function returns a mean or median of pairwise distances of all species in
#' the community, which can be weighted by harmonic means or arithmetic means of abundances
#' of species pairs.
#'
#' @note This function so far only works for linear datasets.
#'
#' @seealso \code{\link{community_overlap_merged}} to see a function that takes
#' intraspecific trait variation (ITV) into account to calculate overlap statistics
#' for a community.
#'
#' @examples
#' library(Ostats)
#'
#' # Load data from web archive
#' dat <- read.csv('https://ndownloader.figshare.com/files/9167548')
#' # Keep only the relevant part of data
#' dat <- dat[dat$siteID %in% c('HARV','JORN') & !is.na(dat$weight), c('siteID', 'taxonID', 'weight')]
#' dat <- dat[!is.na(dat$weight), ]
#' dat$log_weight <- log10(dat$weight)
#'
#' # Calculate median of pairwise distances for the community, weighted by harmonic means
#' # of abundances
#' community_overlap_noitv(traits = as.matrix(dat$log_weight),
#'    sp = factor(dat$taxonID))
#'
#' @export

community_overlap_noitv <- function(traits, sp, output = "median", weight_type= "hmean") {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[stats::complete.cases(dat), ]
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
      if (weight_type == "hmean")
        abund_pairs <- c(abund_pairs, 2/(1/abunds[sp_a] + 1/abunds[sp_b]))
      if (weight_type == "mean")
        abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))

    }
  }

  if (output == "median" && weight_type == "none")
    return(stats::median(overlaps))
  else if (output == "median" && weight_type == "hmean")
    return(matrixStats::weightedMedian(x = overlaps, w = abund_pairs))
  else if (output == "median" && weight_type == "mean")
    return(matrixStats::weightedMedian(x = overlaps, w = abund_pairs))
  else if (output == "mean" && weight_type == "hmean")
    return(stats::weighted.mean(x = overlaps, w = abund_pairs))
  else if (output == "mean" && weight_type == "mean")
    return(stats::weighted.mean(x = overlaps, w = abund_pairs))
  else if (output == "mean" && weight_type == "none")
    return(mean(overlaps))
  else return("invalid arguments for weight_type or output")

}
