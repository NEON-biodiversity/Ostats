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
#'
#' #example use same dataset as ostats (need to change later)
#'
#'library(ggplot2)
#'library(dplyr)
#'
#'overlap_dat <- Ostats_bysite2015 #results of ostat function
#'final_mammal_data <- read.csv('final_NEON_mammal_data.csv', stringsAsFactors = FALSE) #raw data
#'
#'indiv_dat <-data.frame(final_mammal_data$siteID, final_mammal_data$taxonID, final_mammal_data$weight)
#'siteID <- indiv_dat$final_mammal_data.siteID
#'taxonID <- indiv_dat$final_mammal_data.taxonID
#'trait <- log(indiv_dat$final_mammal_data.weight)
#'
#'
#'Ostats_plot(indiv_dat = indiv_dat, siteID = siteID, taxonID = taxonID, trait = trait, overlap_dat = overlap_dat)
#'@export
#'
Ostats_plot<-function(indiv_dat, siteID, taxonID, trait, overlap_dat, sites2use=NULL){

  #put a default that sites2use is null


  ostat_norm<-overlap_dat$overlaps_norm
  #colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = length(unique(taxonID)), replace = TRUE) #it does not work with 26, it asks for 50, but I dont know why
  colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = 50, replace = TRUE)
  theme_set(
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 18),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       legend.position = 'none',
                       strip.background = element_blank()))



  rownames(ostat_norm)
  positions<-factor(match(siteID, rownames(ostat_norm)))

  M<-matrix(data = NA, nrow = length(siteID), ncol = 3)
  M[,1]<-siteID
  M[,2]<-factor(positions)
  M<-data.frame(M)

  for(i in 1:nlevels(positions)){

    M$X3[M$X2==i]<-ostat_norm[i,]
  }

  Ovalues1<-M[,3]
  Ovalues1[is.na(Ovalues1)] = 0
  Ovalues2<-round(x = as.numeric(Ovalues1),digits=4)
  Ovalues3<-cbind(siteID, Ovalues2)
  Ovalues<-paste(Ovalues3[,1], Ovalues3[,2], sep="_")
  indiv_dat$Ovalues<-Ovalues

  #if(is.null(sites2use)==TRUE){
  #PLOT THEM ALL

  ggplot(indiv_dat %>% mutate(siteID = factor(siteID))) +
    stat_density(adjust = 2, size = 1,aes(x = trait, group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
    facet_wrap(~factor(Ovalues), ncol = 3) +
    scale_fill_manual(values = colorvalues) +
    scale_x_continuous(name = colnames(ostat_norm), limits = c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))) +
    scale_y_continuous(name = 'Probability Density',expand = c(0,0))

  #geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 1.5, y = 8.5), color = 'black', data = data.frame(ostat_norm) %>% mutate(siteID = levels(factor(siteID))))


  #} else {
  #PLOT ONLY SOME OF THEM

  #ggplot(filter(indiv_dat, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites2use)))+
  #stat_density(adjust = 2, size = 1,aes(x = trait, group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
  #facet_wrap(~ siteID, ncol = 3) +
  #scale_fill_manual(values = colorvalues) +
  #scale_x_continuous(name = colnames(ostat_norm), limits = c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))) +
  #scale_y_continuous(name = 'Probability Density',expand = c(0,0)) +
  #geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 1.5, y = 8.5), color = 'black', data = overlap_dat %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites2use)))


}

