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
Ostats_plot<-function(indiv_dat, siteID, taxonID, trait, overlap_dat, sites2use=NULL){


  if(is.null(sites2use)==FALSE){

    ostat_norm<-overlap_dat$overlaps_norm
    ostat_norm<-subset(ostat_norm,rownames(ostat_norm)%in%sites2use)
    indiv_dat<-subset(indiv_dat, indiv_dat$siteID%in%sites2use)
    siteID <- indiv_dat$siteID
    taxonID <- indiv_dat$taxonID
    trait <- indiv_dat$logweight


    colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = length(unique(taxonID)), replace = TRUE)
    theme_set(
      theme_bw() + theme(panel.grid = element_blank(),
                         axis.text = element_text(size = 12),
                         axis.title = element_text(size = 18),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position = 'none',
                         strip.background = element_blank()))




    positions<-factor(match(indiv_dat$siteID, rownames(ostat_norm)))

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


    ggplot(indiv_dat %>% mutate(siteID = factor(siteID))) +
      stat_density(adjust = 2, size = 1,aes(x = trait, group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
      facet_wrap(~factor(Ovalues), ncol = 3) +
      scale_fill_manual(values = colorvalues) +
      scale_x_continuous(name = colnames(ostat_norm), limits = c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))) +
      scale_y_continuous(name = 'Probability Density',expand = c(0,0))





  } else {

    ostat_norm<-overlap_dat$overlaps_norm
    colorvalues <- sample(hcl.colors(10, palette = 'viridis'), size = length(unique(taxonID)), replace = TRUE)
    theme_set(
      theme_bw() + theme(panel.grid = element_blank(),
                         axis.text = element_text(size = 12),
                         axis.title = element_text(size = 18),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position = 'none',
                         strip.background = element_blank()))




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




    ggplot(indiv_dat %>% mutate(siteID = factor(siteID))) +
      stat_density(adjust = 2, size = 1,aes(x = trait, group = taxonID, fill=taxonID), alpha = 0.5, geom='polygon', position = 'identity') +
      facet_wrap(~factor(Ovalues), ncol = 3) +
      scale_fill_manual(values = colorvalues) +
      scale_x_continuous(name = colnames(ostat_norm), limits = c(0.5*min(trait,na.rm=TRUE), 1.5*max(trait,na.rm=TRUE))) +
      scale_y_continuous(name = 'Probability Density',expand = c(0,0))




  }

}
