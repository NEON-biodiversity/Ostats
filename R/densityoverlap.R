#' Median pairwise overlap across all species in a community
#'
#' Input: vector of traits, vector of species identities
#'
#' @export
community_overlap <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
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

  median(overlaps)

}

#' Abundance-weighted mean pairwise overlap across all species in a community
#'
#' Input: vector of traits, vector of species identities
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

#' Abundance-weighted median pairwise overlap across all species in a community
#'
#' Input: vector of traits, vector of species identities
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
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
      abund_pairs <- c(abund_pairs, (abunds[sp_a] + abunds[sp_b]))
    }
  }

  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)

}

#' Abundance-weighted (harmonic) median pairwise overlap across all species in a community
#'
#' Input: vector of traits, vector of species identities
#'
#' This version uses the harmonic mean of abundances in each species pair as weights.
#'
#' @export
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

#' Abundance-weighted median pairwise distance between species trait means in a community
#'
#' Input: vector of traits, vector of species identities
#'
#' This reduces each species to a trait mean so that you can compare the inference from overlap to a
#' situation where intraspecific trait variation (ITV) is ignored.
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

#' All pairwise overlap values between pairs of species in a community
#'
#' Input: vector of traits, vector of species identities
#'
#' @export
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

