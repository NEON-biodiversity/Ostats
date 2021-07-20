# Ostats: O-statistics, or Pairwise Community-Level Niche Cverlap Statistics

The `Ostats` package calculates overlap statistic to measure the degree of community-level trait overlap by fitting nonparametric kernel density functions to each species' trait distribution and calculating their areas of overlap (Mouillot et al. 2005, Geange et al. 2011, Read et al. 2018). For instance, the median pairwise overlap for a community is calculated by first determining the overlap of each species pair in trait space, and then taking the median overlap of each species pair in a community. This median overlap value is called the O-statistic (O for overlap). The `Ostats()` function calculates separate univariate overlap statistics for each trait, while the `Ostats_multivariate()` function calculates a single multivariate overlap statistic for all traits. O-statistics can be evaluated against null models to obtain standardized effect sizes. Grady et al. (2018) provide a teaching module that goes into detail about how to interpret the results presented in Read et al. (2018).

Ostats is part of the collaborative Macrosystems Biodiversity Project *Local- to continental-scale drivers of biodiversity across the National Ecological Observatory Network (NEON)*. For more information on this project, see the [Macrosystems Biodiversity Website](https://neon-biodiversity.github.io/).

## Authors

* Quentin D. Read
* Arya Yue
* Isadora E. Fluck
* Benjamin Baiser
* John M. Grady
* Phoebe L. Zarnetske
* Sydne Record

## How to install

Type the following into your R prompt.

```
remotes::install_github('NEON-biodiversity/Ostats')
```

## Getting started

See [the Ostats introduction vignette](https://neon-biodiversity.github.io/Ostats/articles/Ostats-introduction.html) for more information.  

## Funding 

National Science Foundation Division of Environmental Biology (Population & Community Ecology, MacroSysBIO & NEON-Enabled Science); Awards to P.L. Zarnetske (Michigan State University); S. Record (Bryn Mawr College); Ben Baiser (University of Florida); Angela Strecker (Western Washington University).

The National Ecological Observatory Network is a program sponsored by the National Science Foundation and operated under cooperative agreement by Battelle. This material is based in part upon work supported by the National Science Foundation through the NEON Program.

[![NSF-1926567](https://img.shields.io/badge/NSF-1926567-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1926567)
[![NSF-1926568](https://img.shields.io/badge/NSF-1926568-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1926568)
[![NSF-1926569](https://img.shields.io/badge/NSF-1926569-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1926569)
[![NSF-1926610](https://img.shields.io/badge/NSF-1926610-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1926610)

## References

Geange, S.W., S. Pledger, K.C. Burns, and J.S. Shima. 2011. A unified analysis of niche overlap incorporating data of different types. *Methods in Ecology and Evolution* 2(2):175-184. https://doi.org/10.1111/j.2041-210X.2010.00070.x

Grady, J.M., Q.D. Read, S. Record, P.L. Zarnetske, B. Baiser, K. Thorne, and J. Belmaker. 2018. Size, niches, and the latitudinal diversity gradient. *Teaching Issues and Experiments in Ecology* 14: Figure Set #1.

Mouillot, D., W. Stubbs, M. Faure, O. Dumay, J.A. Tomasini, J.B. Wilson, and T. Do Chi. 2005. Niche overlap estimated based on quantitative functional traits: A new family of non-parametric indices. *Oecologia* 145(3):345-353. https://doi.org/10.1007/s00442-005-0151-z

Read, Q.D., J.M. Grady, P.L. Zarnetske, S. Record, B. Baiser, J. Belmaker, M.-N. Tuanmu, A. Strecker, L. Beaudrot, and K.M. Thibault. 2018. Among-species overlap in rodent body size distributions predicts species richness along a temperature gradient. *Ecography* 41(10):1718-1727. https://doi.org/10.1111/ecog.03641

