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
#'@param overlap_dat a data frame where the first column contains the names of each community and the second column contains the values of total overlap among species in each of these communities obtained from the \code{\link{Ostats}} function  NEED TO FIGURE OUT HOW TO USE OVERLAP OUTPUT
#'@param trait_name a character that specifies the name of the trait, to be
#'  labelled on the x axis.
#'
#'
#'@return Density plots of species trait distribution plotted on the same graph
#'  for each community to show how they overlap each other. The overlap value obtained as output from the function \code{\link{Ostats}}, is labelled on each community graph.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'
#' #FUNCTION NOT COMPLETE!!!
#' #example use same dataset as ostats (need to change later)
#'
#' library(tidyverse)
#' overlap_dat <- read.csv('overlap_dat.csv')
#' indiv_dat <- read.csv('indiv_dat.csv')
#' indiv_dat <- indiv_dat %>%filter(siteID %in% c('HARV','JORN'))
#' overlap_dat <- overlap_dat %>%filter(siteID %in% c('HARV','JORN'))
#'
#' indiv_dat <- mutate(indiv_dat, log_weight = log10(weight))
#'
#' #change column name to test whether the function works
#' names(indiv_dat) <- c("siteID", "species","weight", "log_weight")
#'
#' Ostats_plot(indiv_dat = indiv_dat,siteID=indiv_dat$siteID, taxonID = indiv_dat$species, trait = indiv_dat$log_weight,overlap_dat = overlap_dat, trait_name = "log_weight")
#'@export
#'
Ostats_plot<-function(indiv_dat, siteID, taxonID, trait, overlap_dat, trait_name){
  colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = 24, replace = TRUE)
  theme_set(
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 18),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       legend.position = 'none',
                       strip.background = element_blank())
  )
  ggplot(indiv_dat %>% mutate(siteID = factor(siteID))) +
    stat_density(adjust = 2, size = 1,aes(x = trait, group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
    facet_wrap(~ siteID, ncol = 1) +
    scale_fill_manual(values = colorvalues) +
    scale_x_continuous(name = trait_name) +
    scale_y_continuous(name = 'Probability Density',expand = c(0,0)) +
    geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 1.5, y = 8.5), color = 'black', data = overlap_dat %>% mutate(siteID = factor(siteID)))
}



