#'@title  Plotting Community Overlap
#'
#'@description This function plots the overlap of traits among
#'  species for each community.
#'
#'@param indiv_dat Individual data: a data frame containing individual measurments for a certain
#'traits,its species identity, and the community identity it belongs to.
#'@param siteID Site identity: a column in indiv_dat data frame that indicates the names of each community.
#'@param taxonID Taxon identity: a column in indiv_dat data frame that indicates species or taxa names.
#'@param trait The trait you want to overlap among species: a column in indiv_dat data frame containing trait measurements for each individual.
#'@param overlap_dat an object containing the results from the Ostats function - see \code{\link{Ostats}}
#'@param sites2use a vector that select the sites you want to plot. If NULL, the function will plot all the sites.
#'@param n_col Number of columns for layout of individual panels. Default is 3.
#'@param colorvalues Vector of color values for the density polygons. Defaults to a viridis palette if none provided.
#'
#'
#'@return Density plots of species trait distribution plotted on the same graph
#'  for each community to show how they overlap each other. The overlap value obtained as output from the function \code{\link{Ostats}}, is labelled on each community graph.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'library(tidyverse)
#'library(ggplot2)
#'library(dplyr)
#'
#'overlap_dat <- Ostats_example #your results from the Ostat function - see \code{\link{Ostats}}
#'indiv_dat <- read_csv('https://ndownloader.figshare.com/files/9167548', col_names = T)
#'siteID <- indiv_dat$siteID
#'taxonID <- indiv_dat$taxonID
#'trait <- indiv_dat$logweight
#'#to plot only selected sites:
#'sites2use<- c('BART','KONZ','JORN') #used in the example
#'
#'#to plot all sites:
#'sites2use<- NULL

#'Ostats_plot(indiv_dat = indiv_dat, siteID = siteID, taxonID = taxonID, trait = trait, overlap_dat = overlap_dat, sites2use = sites2use)
#'@export
#'
Ostats_plot<-function(indiv_dat, siteID, taxonID, trait, overlap_dat, sites2use = NULL, n_col = 3, colorvalues = NULL) {

  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(sites2use)) {
    sites2use <- unique(indiv_dat$siteID)
  }


  ostat_norm<-overlap_dat$overlaps_norm
  ostat_norm <- subset(ostat_norm, rownames(ostat_norm) %in% sites2use)
  indiv_dat <- subset(indiv_dat, siteID %in% sites2use)
  siteID <- indiv_dat$siteID #FIXME this will not work if indiv_dat does not contain siteID. Also, this will overwrite the siteID argument the user supplies.
  taxonID <- indiv_dat$taxonID #FIXME this will not work if indiv_dat does not contain taxonID. Also, this will overwrite the taxonID argument the user supplies.
  trait <- indiv_dat$logweight #FIXME this will not work if the trait is not called logweight. Also, this will overwrite the trait argument the user supplies.

  # If a color vector is not provided, create a default palette.
  if (is.null(colorvalues)) {
    colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = length(unique(taxonID)), replace = TRUE)
  }

  ggplot2::theme_set(
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                         axis.text = ggplot2::element_text(size = 12),
                                         axis.title = ggplot2::element_text(size = 18),
                                         axis.text.y = ggplot2::element_blank(),
                                         axis.ticks.y = ggplot2::element_blank(),
                                         legend.position = 'none',
                                         strip.background = ggplot2::element_blank()))

  overlap_labels <- data.frame(siteID = row.names(ostat_norm),
                               lab = paste('Overlap =', round(ostat_norm[,1], 2)))
  indiv_dat$siteID <- factor(indiv_dat$siteID)


  ggplot2::ggplot(indiv_dat) +
    ggplot2::stat_density(adjust = 2, size = 1, ggplot2::aes(x = trait, group = taxonID, fill = taxonID), alpha = 0.5, geom='polygon', position = 'identity') + #FIXME add options to make alpha, adjust, and size variable.
    ggplot2::facet_wrap(~ siteID, ncol = n_col) +
    ggplot2::scale_fill_manual(values = colorvalues) +
    ggplot2::geom_text(ggplot2::aes(label = lab), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1) +
    ggplot2::scale_x_continuous(name = 'Trait value', limits = c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))) + #FIXME add options for the user to manually supply limits, use these as defaults otherwise.
    ggplot2::scale_y_continuous(name = 'Probability Density', expand = c(0,0)) #FIXME add options to change default x and y axis names, adding sensible defaults if none provided



}
