# Prepare datasets for archiving with the package.

# Raw weight measurements
indiv_dat <- read.csv('https://ndownloader.figshare.com/files/9167548')
indiv_dat <- indiv_dat[, c('siteName', 'siteID', 'plotID', 'taxonID', 'weight', 'sex', 'lifeStage')]
small_mammal_data <- indiv_dat

# Precalculated overlap statistics
Ostats_bysite2015 <- Ostats(traits = log10(as.matrix(indiv_dat[,'weight',drop=FALSE])),
                            plots = factor(indiv_dat$siteID),
                            sp = factor(indiv_dat$taxonID),
                            random_seed = 111,
                            nperm = 999)
small_mammal_Ostats <- Ostats_bysite2015

save(small_mammal_data, file = 'data/small_mammal_data.rda')
save(small_mammal_Ostats, file = 'data/small_mammal_Ostats.rda')

# Pitcher traits (added by QDR 1 Apr 2021)
pitcher <- read.csv('~/Downloads/pitcher_traits.csv')

# Anonymize, get rid of long and lat, and get only a subset of the sites.
pitcher_traits <- pitcher[pitcher$continent %in% 'N America', !names(pitcher) %in% c('continent', 'latitude', 'longitude')]

# Right now I will not add any "noise" to data but could do so later if needed.
save(pitcher_traits, file = 'data/pitcher_traits.rda')
