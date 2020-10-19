#'@title  Plotting Community Overlap
#'
#'@description This function plots the overlap of traits among
#'  species for each community.
#'
#'@param indiv_dat Individual data: a data frame containing individual measurments for a certain
#'trait,its species identity, and the community identity it belongs to.
#'@param plots Site identity: a column in indiv_dat data frame that indicates the names of each community.
#'@param sp Taxon identity: a column in indiv_dat data frame that indicates species or taxa names.
#'@param traits you want to overlap among species: a column in indiv_dat data frame containing traits measurements for each individual.
#'@param overlap_dat an object containing the results from the Ostats function - see \code{\link{Ostats}}
#'@param sites2use a vector that select the sites you want to plot. If NULL, the function will plot all the sites.
#'@param n_col Number of columns for layout of individual panels. Default is 1.
#'@param colorvalues Vector of color values for the density polygons. Defaults to a viridis palette if none provided.
#'@param alpha_o defines the colors trasparency level for the density polygons. Default is 0.5
#'@param adjust multiplicate the bandwidth adjustment of the density polygons. The less, the tiny your density polygons will be. Default is 2.
#'@param limits_x the limits (min and max values) of the x axis. Default is \code{c(0.5*min(traits,na.rm=TRUE), 1.5*max(traits,na.rm=TRUE))}
#'@param scale If you want the scale of x, y or both x and y axis to be adjusted according to each site density probability set the argument to "free_x", "free_y" or "free" respectively. Default= "fixed" which uses the same scale across all sites.
#'@param name_x a character indicating the name of your x axis (i.e. the name of your trait). Default is 'traits value'
#'@param name_y a character indicating the name of your y axis. Default is 'Probability Density'
#'@param means if TRUE it plot traits means for each species in an additionlan plot column next to the traits distribution plots for each site. Default is FALSE, which make the function plot only the traits distribution for each site.
#'@return Density plots of species trait distribution plotted on the same graph
#'  for each community to show how they overlap each other. The overlap value obtained as output from the function \code{\link{Ostats}}, is labelled on each community graph.
#'
#'@seealso \code{\link{Ostats}} to Calculate O-statistics (community-level
#'  pairwise niche overlap statistics)
#'
#'@examples
#'#load data:
#'#indiv_dat <- read_csv('https://ndownloader.figshare.com/files/9167548', col_names = T)
#'
#'#set the arguments:
#'#overlap_dat <- Ostats_bysite2015 #your results from the Ostat function - see \code{\link{Ostats}}
#'#plots <- indiv_dat$siteID
#'#sp <- indiv_dat$taxonID
#'#traits <- indiv_dat$logweight
#'
#'#to plot only selected sites:
#'#sites2use<- c('BART','KONZ','JORN')
#'
#'#to plot all sites:
#'sites2use<- NULL
#'
#'
#'#Ostats_plot(indiv_dat = indiv_dat,
#'            plots = plots, sp = sp, traits = traits,
#'            overlap_dat = overlap_dat,
#'            sites2use = sites2use, means=T)


#'@export
#'
Ostats_plot<-function(indiv_dat,
                      plots,
                      sp,
                      traits,
                      overlap_dat,
                      sites2use = NULL,
                      n_col=1,
                      scale = "fixed",
                      colorvalues = NULL,
                      alpha = 0.5,
                      adjust = 2,
                      limits_x =c(0.5*min(traits,na.rm=TRUE), 1.5*max(traits,na.rm=TRUE)),
                      name_x = 'traits value',
                      name_y = 'Probability Density',
                      means=FALSE) {



  # Unless a subset of sites is provided, use all sites in dataset.
  if (is.null(sites2use)) {
    sites2use <- unique(plots)
  }

  #filter only for sites2use
  ostat_norm<-overlap_dat$overlaps_norm
  ostat_norm <- subset(ostat_norm, rownames(ostat_norm) %in% sites2use)

  traits <- subset(traits, plots %in% sites2use)
  indiv_dat <- subset(indiv_dat, plots %in% sites2use)
  sp<-subset(sp, plots %in% sites2use) #filter the taxons in the sites2use
  plots<-subset(plots, plots %in% sites2use)

  #organize data in a table
  table_traits_taxon<-data.frame(traits, sp, plots)


  # If the user want to plot the traits means.
  if(means==TRUE){
    #values per species
    taxon_mean<-aggregate(table_traits_taxon[,1], list(table_traits_taxon[,2]), mean, na.rm = TRUE)

    #make a column with repeted means for each specie
    table_all<-cbind(table_traits_taxon, table_traits_taxon[,2])
    table_all<-data.frame(table_all,stringsAsFactors = FALSE)

    n<-length(unique(taxon_mean$Group.1))
    for(i in 1:n){
      name<-unique(taxon_mean$Group.1)[i]
      T_F<-taxon_mean[,1]==name
      mean_name<-taxon_mean[,2][T_F]
      table_all <- within(table_all, table_traits_taxon...2.[sp==name] <- mean_name)
    }

    names(table_all)[names(table_all) == "table_traits_taxon...2."] <- "means"

  }


  # If a color vector is not provided, create a default palette.
  if (is.null(colorvalues)) {
    colorvalues <- sample(rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,
                                  alpha, rev = FALSE), size = length(unique(table_all$sp)), replace = TRUE)
  }

  names(colorvalues) <- unique(taxon_mean$Group.1)

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

  ggplot_dist<-ggplot2::ggplot(table_all) +
    ggplot2::stat_density(adjust = adjust, ggplot2::aes(x = traits, group = sp, fill=sp), alpha = alpha, geom='polygon', position = 'identity') +
    ggplot2::facet_wrap(~ plots, ncol=n_col, nrow = length(sites2use), scales = scale) +
    ggplot2::scale_fill_manual(values = colorvalues) +
    ggplot2::geom_text(ggplot2::aes(label = lab), data = overlap_labels, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1) +
    ggplot2::scale_x_continuous(name = name_x, limits = limits_x) +
    ggplot2::scale_y_continuous(name = name_y, expand = c(0,0))


  if (means==TRUE) {
    ggplot_means<-ggplot2::ggplot(table_all)+
      ggplot2::geom_vline(data=table_all, ggplot2::aes(xintercept=as.numeric(means),  colour=sp,  group=sp, alpha = alpha), size=0.5)+
      ggplot2::facet_wrap(~ plots, ncol=n_col ,nrow = length(sites2use), scales = scale) +
      ggplot2::scale_colour_manual(values = colorvalues) +
      ggplot2::scale_x_continuous(name = name_x, limits = limits_x) +
      ggplot2::scale_y_continuous(expand = c(0,0))


  }



  if (means==TRUE){
    gridExtra::grid.arrange(ggplot_dist, ggplot_means, ncol=2)

  } else {

    ggplot_dist

  }
}
