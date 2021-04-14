library(Ostats)
context("Ostats_multivariate")

# Set up iris data for test
iris_traits <- as.matrix(iris[,1:4])

# Test 1: iris

result1 <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE)
result1 <- as.numeric(result1$overlaps_norm)
expected1 <- c(0.62, 0.72)

test_that (
  "Ostats_multivariate returns expected output",
  {
    expect_equal(result1, expected1, tolerance = 0.01)
  }
)

# Test 2: run multivariate data on single trait
result2 <- Ostats_multivariate(traits = iris_traits[,1,drop=FALSE], plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE)
expected2 <-  Ostats(traits = iris_traits[,1,drop=FALSE], plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE)

# expect_equivalent ignores names.
test_that (
  "Ostats_multivariate returns the same output as Ostats if one trait is selected",
  {
    expect_equivalent(result2$overlaps_norm, expected2$overlaps_norm, tolerance = 0.01)
  }
)
