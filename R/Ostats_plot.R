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
#'@param alpha_o defines the colors trasparency level for the density polygons. Default is 0.5
#'@param adjust_o multiplicate the bandwidth adjustment of the density polygons. The less, the tiny your density polygons will be. Default is 2.
#'@param limits_xo the limits (min and max values) of the x axis. Default is \code{c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))}
#'@param scale_o If you want the scale of x, y or both x and y axis to be adjusted according to each site density probability set the argument to "free_x", "free_y" or "free" respectively. Default= "fixed" which uses the same scale across all sites.
#'@param name_x a character indicating the name of your x axis (i.e. the name of your trait). Default is 'Trait value'
#'@param name_y a character indicating the name of your y axis. Default is 'Probability Density'
#'@param means_o if TRUE it plot traits means for each species in an additionlan plot column next to the traits distribution plots for each site. Default is FALSE, which make the function plot only the traits distribution for each site.
#'@return Density plots of species trait distribution plotted on the same graph
#'  for each community to show how they overlap each other. The overlap value obtained as output from the function \code{\link{Ostats}}, is labelled on each community graph.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'library(tidyverse)
#'library(Ostats)
#'
#'indiv_dat <- read_csv('https://ndownloader.figshare.com/files/9167548', col_names = T)
#'
#'Ostats_plot(indiv_dat = indiv_dat, siteID = indiv_dat$siteID, taxonID = indiv_dat$taxonID, trait = indiv_dat$logweight, overlap_dat = Ostats_example, sites2use = c('HARV','JORN'))

#'@export
#'
Ostats_plot<-function(indiv_dat, siteID, taxonID, trait, overlap_dat, sites2use = NULL, n_col=3, scale_o = "fixed", colorvalues = NULL, alpha_o = 0.5, adjust_o = 2, limits_xo =c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE)), name_x = 'Trait value', name_y = 'Probability Density', means_o=FALSE) {



  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(sites2use)) {
    sites2use <- unique(indiv_dat$siteID)
  }

  #filter ostats results
  ostat_norm<-overlap_dat$overlaps_norm
  ostat_norm <- subset(ostat_norm, rownames(ostat_norm) %in% sites2use)

  #filter only for sites2use
  trait <- subset(trait, siteID %in% sites2use)
  indiv_dat <- subset(indiv_dat, siteID %in% sites2use)
  taxonID<-subset(taxonID, siteID %in% sites2use) #filter the taxons in the sites2use
  siteID<-subset(siteID, siteID %in% sites2use)

  #organize data in a table
  teble_trait_taxon<-data.frame(trait, taxonID, siteID)


  # If the user want to plot the trait means.
  if(means_o==TRUE){
    #values per species
    taxon_mean<-aggregate(teble_trait_taxon[,1], list(teble_trait_taxon[,2]), mean, na.rm = TRUE)

    #make a column with repeted means for each specie
    table_all<-cbind(teble_trait_taxon, teble_trait_taxon[,2])
    for(i in 1:length(unique(taxonID))){
      name<-unique(table_all[,4])[i]
      T_F<-taxon_mean[,1]==name
      mean_name<-taxon_mean[,2][T_F]
      table_all[,4]<-gsub(name, mean_name, table_all[,4])
    }
    names(table_all)[names(table_all) == "teble_trait_taxon[, 2]"] <- "means_o"


  }




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


  ggplot_dist<-ggplot2::ggplot(table_all) +
    ggplot2::stat_density(adjust = adjust_o, ggplot2::aes(x = trait, group = taxonID, fill=taxonID), alpha = alpha_o, geom='polygon', position = 'identity') +
    ggplot2::facet_wrap(~ siteID, ncol=n_col, nrow = length(sites2use), scales = scale_o) +
    ggplot2::scale_fill_manual(values = colorvalues) +
    ggplot2::geom_text(ggplot2::aes(label = lab), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1) +
    ggplot2::scale_x_continuous(name = name_x, limits = limits_xo) +
    ggplot2::scale_y_continuous(name = name_y, expand = c(0,0))






  if (means_o==TRUE) {
    ggplot_means<-ggplot2::ggplot(table_all)+
      ggplot2::geom_vline(data=table_all, ggplot2::aes(xintercept=as.numeric(means_o),  colour=taxonID,  group=taxonID, alpha = alpha_o), size=0.5)+
      ggplot2::facet_wrap(~ siteID, ncol=n_col ,nrow = length(sites2use), scales = scale_o) +
      ggplot2::scale_color_manual(values = colorvalues) +
      ggplot2::geom_text(ggplot2::aes(label = lab), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1) +
      ggplot2::scale_x_continuous(name = name_x, limits = limits_xo) +
      ggplot2::scale_y_continuous(name = name_y, expand = c(0,0))


  }


  #final

  if (means_o==TRUE){
    gridExtra::grid.arrange(ggplot_dist, ggplot_means, ncol=2)

  } else {

    ggplot_dist

  }
}
