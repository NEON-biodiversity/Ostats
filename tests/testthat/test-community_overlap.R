library(Ostats)
context("community_overlap")

# Set up data for test.
# Use small mammal data from the vignette.
# For this test we will only use HARV.

dat <- small_mammal_data[small_mammal_data$siteID %in% c('HARV', 'JORN'),
                         c('siteID', 'taxonID', 'weight')]
dat <- dat[!is.na(dat$weight), ]
dat$log_weight <- log10(dat$weight)

trait_harv <- as.matrix(dat[dat$siteID %in% 'HARV', 'log_weight', drop = FALSE])
spp_harv <- factor(dat$taxonID[dat$siteID %in% 'HARV'])

# Test 1. Verify correct output format if raw = FALSE
result1 <- community_overlap(traits = trait_harv, sp = spp_harv)
expected1 <- 0.8946

test_that (
  "community_overlap() returns expected output when raw = FALSE",
  {
    expect_equivalent(result1, expected1, tolerance = 0.001)
  }
)

# Test 2. Verify correct output format if raw = TRUE
result2 <- community_overlap(traits = trait_harv, sp = spp_harv, raw = TRUE)

test_that (
  "community_overlap() returns expected output when raw = TRUE",
  {
    expect_equivalent(result2[['value']], expected1, tolerance = 0.001)
    expect_s3_class(result2[['raw']], 'data.frame')
    expect_length(result2[['raw']][['overlap']], 10)
  }
)

