library(Ostats)
context("Ostats")

# Set up data for test.
# Use small mammal data from the vignette.

dat <- small_mammal_data[small_mammal_data$siteID %in% c('HARV', 'JORN'),
                         c('siteID', 'taxonID', 'weight')]
dat <- dat[!is.na(dat$weight), ]
dat$log_weight <- log10(dat$weight)

# Test 1: one community with multiple species
result1 <- Ostats(traits = as.matrix(dat[, 'log_weight', drop = FALSE]), plots = factor(dat$siteID), sp = factor(dat$taxonID), run_null_model = FALSE)$overlaps_norm
expected1 <- matrix(c(0.8946, 0.0183), nrow = 2)

test_that (
  "Ostats returns expected output",
  {
    expect_equivalent(result1, expected1, tolerance = 0.001)
  }
)

# Test 2: multiple communities, one of which has only one species
dat2 <- dat[dat$siteID %in% 'HARV' | dat$taxonID %in% 'CHPE', ]
result2 <- Ostats(traits = as.matrix(dat2[, 'log_weight', drop = FALSE]), plots = factor(dat2$siteID), sp = factor(dat2$taxonID), run_null_model = FALSE)$overlaps_norm
expected2 <- matrix(c(0.8946, NA), nrow = 2)

test_that (
  "Ostats deals with communities with only one species",
  {
    expect_equivalent(result2, expected2, tolerance = 0.001)
  }
)

# Test 3: species with insufficient data for density estimate
# Edge case where every species in the community has only one individual. Should return NA.

dat3 <- do.call(rbind, lapply(split(dat, interaction(dat$siteID, dat$taxonID), drop = TRUE),
                              function(x) x[1,]))
result3 <- Ostats(traits = as.matrix(dat3[, 'log_weight', drop = FALSE]), plots = factor(dat3$siteID), sp = factor(dat3$taxonID), run_null_model = FALSE)$overlaps_norm
expected3 <- matrix(c(NA, NA), nrow = 2)


test_that (
  "Ostats deals with species with insufficient data",
  {
    expect_equivalent(result3, expected3, tolerance = 0.001)
  }
)

### Additional tests from vignette
# Test 4: Additional arguments to density, vignette line 111
result4 <- Ostats(traits = as.matrix(dat[,'log_weight', drop = FALSE]),
                  sp = factor(dat$taxonID),
                  plots = factor(dat$siteID),
                  run_null_model = FALSE,
                  density_args=list(bw = 'nrd0', adjust = 2, n=200))$overlaps_norm
expected4 <- matrix(c(0.895,0.0188), nrow = 2)

test_that (
  "Ostats handles the different density arguments correctly",
  {
    expect_equivalent(result4, expected4, tolerance = 0.001)
  }
)

# Test 5: Circular data example from vignette line 162
# Is the hourly circular data handled correctly?
result5 <- Ostats(traits = as.matrix(ant_data[, 'time', drop = FALSE]),
                  sp = factor(ant_data$species),
                  plots = factor(ant_data$chamber),
                  discrete = TRUE,
                  circular = TRUE,
                  unique_values = 0:23,
                  run_null_model = FALSE)$overlaps_norm

expected5 <- matrix(c(0.6803, 0.6348), nrow = 2)

test_that (
  "discrete circular data is handled correctly",
  {
    expect_equivalent(result5, expected5, tolerance = 0.001)
  }
)

# Test 6: Non-circular discrete example. Use fake data with 0.75 overlap.
result6 <- Ostats(traits = matrix(c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,5,5), ncol = 1),
                  sp = factor(rep(c('a','b'), c(8, 8))),
                  plots = factor(rep('a', 16)),
                  discrete = TRUE,
                  circular = FALSE,
                  run_null_model = FALSE)$overlaps_norm

expected6 <- 0.75

test_that (
  "discrete non-circular data is handled correctly",
  {
    expect_equivalent(result6, expected6, tolerance = 0.001)
  }
)
