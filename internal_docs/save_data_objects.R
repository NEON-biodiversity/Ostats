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

