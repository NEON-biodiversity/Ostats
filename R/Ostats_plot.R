#'@title  Plotting Community Overlap
#'
#'@description This function plots the overlap of traits among
#'  species for each community.
#'
#'@param plots Site identity: a vector of names of each community.
#'@param sp Taxon identity: a vector of species or taxa names.
#'@param traits A vector of trait measurements for each individual.
#'@param overlap_dat an object containing the output of \code{\link{Ostats}}
#'@param sites2use a vector of sites to plot. If NULL, the function will plot all the sites.
#'@param n_col Number of columns for layout of individual panels. Default is 1.
#'@param colorvalues Vector of color values for the density polygons. Defaults to a rainbow palette if none provided.
#'@param alpha defines the transparency level for the density polygons. Default is 0.5
#'@param adjust multiplicate the bandwidth adjustment of the density polygons. The less, the tiny your density polygons will be. Default is 2.
#'@param limits_x the limits (min and max values) of the x axis. Default is \code{c(0.5*min(traits,na.rm=TRUE), 1.5*max(traits,na.rm=TRUE))}
#'@param scale If you want the scale of x, y or both x and y axis to be adjusted according to each site density probability set the argument to "free_x", "free_y" or "free" respectively. Default = "fixed" which uses the same scale across all sites.
#'@param name_x a character indicating the name of your x axis (i.e. the name of your trait). Default is 'trait value'
#'@param name_y a character indicating the name of your y axis. Default is 'probability density'
#'@param means if TRUE it plot traits means for each species in an additional plot column next to the traits distribution plots for each site. Default is FALSE, which make the function plot only the traits distribution for each site.
#'@return Density plots of species trait distribution plotted on the same graph
#'  for each community to show how they overlap each other.
#'  The overlap value obtained as output from \code{\link{Ostats}} is labelled on each community graph.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'# set the arguments:
#'plots <- small_mammal_data$siteID
#'sp <- small_mammal_data$taxonID
#'traits <- log10(small_mammal_data$weight)
#'
#'# to plot only selected sites:
#'sites2use <- c('BART','KONZ','JORN')
#'
#'
#'Ostats_plot(plots = plots, sp = sp, traits = traits,
#'            overlap_dat = small_mammal_Ostats,
#'            sites2use = sites2use, means = TRUE)


#'@export
#'
Ostats_plot<-function(plots,
                      sp,
                      traits,
                      overlap_dat,
                      sites2use = NULL,
                      n_col = 1,
                      scale = "fixed",
                      colorvalues = NULL,
                      alpha = 0.5,
                      adjust = 2,
                      limits_x = c(0.5*min(traits,na.rm=TRUE), 1.5*max(traits,na.rm=TRUE)),
                      name_x = 'trait value',
                      name_y = 'probability density',
                      means=FALSE) {



  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(sites2use)) {
    sites2use <- unique(plots)
  }

  #filter only for sites2use
  ostat_norm<-overlap_dat$overlaps_norm
  ostat_norm <- subset(ostat_norm, rownames(ostat_norm) %in% sites2use)

  traits <- subset(traits, plots %in% sites2use)
  sp<-subset(sp, plots %in% sites2use)
  plots<-subset(plots, plots %in% sites2use)

  plot_dat <- data.frame(traits = traits, sp = sp, plots = plots)

  # Calculate mean value by taxon.
  taxon_mean <- stats::aggregate(traits, list(sp, plots), mean, na.rm = TRUE)
  names(taxon_mean) <- c('sp', 'plots', 'means')

  # If a color vector is not provided, create a default palette.
  if (is.null(colorvalues)) {
    colorvalues <- sample(grDevices::rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,
                                  alpha, rev = FALSE), size = length(unique(sp)), replace = TRUE)
  }

  names(colorvalues) <- unique(sp)

  ggplot2::theme_set(
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                         axis.text = ggplot2::element_text(size = 12),
                                         axis.title = ggplot2::element_text(size = 12),
                                         axis.text.y = ggplot2::element_blank(),
                                         axis.ticks.y = ggplot2::element_blank(),
                                         legend.position = 'none',
                                         strip.background = ggplot2::element_blank()))

  overlap_labels <- data.frame(plots = row.names(ostat_norm),
                               lab = paste('Overlap =', round(ostat_norm[,1], 2)))

  ggplot_dist<-ggplot2::ggplot(plot_dat) +
    ggplot2::stat_density(adjust = adjust, ggplot2::aes(x = traits, group = sp, fill = sp), alpha = alpha, geom='polygon', position = 'identity') +
    ggplot2::facet_wrap(~ plots, ncol=n_col, nrow = length(sites2use), scales = scale) +
    ggplot2::scale_fill_manual(values = colorvalues) +
    ggplot2::geom_text(ggplot2::aes_string(label = 'lab'), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1) +
    ggplot2::scale_x_continuous(name = name_x, limits = limits_x) +
    ggplot2::scale_y_continuous(name = name_y, expand = c(0,0))


  if (means) {
    ggplot_means<-ggplot2::ggplot(taxon_mean) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=means,  colour=sp,  group=sp, alpha = alpha), size=0.5)+
      ggplot2::facet_wrap(~ plots, ncol = n_col ,nrow = length(sites2use), scales = scale) +
      ggplot2::scale_colour_manual(values = colorvalues) +
      ggplot2::scale_x_continuous(name = name_x, limits = limits_x) +
      ggplot2::scale_y_continuous(expand = c(0,0))
  }

  if (means){
    gridExtra::grid.arrange(ggplot_dist, ggplot_means, ncol=2)
  } else {
    ggplot_dist
  }
}
