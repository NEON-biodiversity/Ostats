library(Ostats)
context("Ostats")

# Set up data for test.
# Use small mammal data from the vignette.

dat <- small_mammal_data[small_mammal_data$siteID %in% c('HARV', 'JORN'),
                         c('siteID', 'taxonID', 'weight')]
dat <- dat[!is.na(dat$weight), ]
dat$log_weight <- log10(dat$weight)

# Test 1: one community with multiple species
result1 <- Ostats(traits = as.matrix(dat[, 'log_weight', drop = FALSE]), plots = factor(dat$siteID), sp = factor(dat$taxonID), nperm = 1, random_seed = 919)$overlaps_norm
expected1 <- structure(c(0.8946, 0.0183), .Dim = 2:1, .Dimnames = list(
  c("HARV", "JORN"), "log_weight"))

test_that (
  "Ostats returns expected output",
  {
    expect_equal(result1, expected1, tolerance = 0.001)
  }
)

# Test 2: multiple communities, one of which has only one species
dat2 <- dat[dat$siteID %in% 'HARV' | dat$taxonID %in% 'CHPE', ]
result2 <- Ostats(traits = as.matrix(dat2[, 'log_weight', drop = FALSE]), plots = factor(dat2$siteID), sp = factor(dat2$taxonID), nperm = 1, random_seed = 919)$overlaps_norm
expected2 <- structure(c(0.8946, NA), .Dim = 2:1, .Dimnames = list(
  c("HARV", "JORN"), "log_weight"))

test_that (
  "Ostats deals with communities with only one species",
  {
    expect_equal(result2, expected2, tolerance = 0.001)
  }
)

# Test 3: species with insufficient data for density estimate
# Edge case where every species in the community has only one individual. Should return NA.

dat3 <- do.call(rbind, lapply(split(dat, interaction(dat$siteID, dat$taxonID), drop = TRUE),
                              function(x) x[1,]))
result3 <- Ostats(traits = as.matrix(dat3[, 'log_weight', drop = FALSE]), plots = factor(dat3$siteID), sp = factor(dat3$taxonID), nperm = 1, random_seed = 919)$overlaps_norm
expected3 <- structure(c(NA, NA), .Dim = 2:1, .Dimnames = list(
  c("HARV", "JORN"), "log_weight"))


test_that (
  "Ostats deals with species with insufficient data",
  {
    expect_equal(result3, expected3, tolerance = 0.001)
  }
)
