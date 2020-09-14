# Ostats: O-statistics, or pairwise community-level niche overlap statistics

The Ostats package calculates overlap statistic to measure the degree of community-level trait overlap by fitting nonparametric kernel density functions to each species' trait distribution and calculating their areas of overlap (Mouillot et al. 2005, Geange et al. 2011, Read et al. 2018). For instance, the median pairwise overlap for a community is calculated by first determining the overlap of each species pair in trait space, and then taking the median overlap of each species pair in a community. This function will only work with univariate trait data. Grady et al. (2018) provides a teaching module that goes into detail about how to interpret the O-stats results presented in Read et al. (2018).


## Authors

* Quentin D. Read, SESYNC
* Arya Yue
* Isadora E. Fluck
* Benjaming Baiser
* John M. Grady
* Phoebe L. Zarnetske
* Sydne Record

## How to install

Type the following into your R prompt.

```
devtools::install_github('NEON-biodiversity/Ostats')
```

## References

Read, Q. D., J. M. Grady, P. L. Zarnetske, S. Record, B. Baiser, J. Belmaker, M.-N. Tuanmu, A. Strecker, L. Beaudrot, and K. M. Thibault. 2018. Among-species overlap in rodent body size distributions predicts species richness along a temperature gradient. Ecography. DOI:10.1111/ecog.03641
