library(Ostats)
context("Ostats_multivariate")

# Set up iris data for test
iris_traits <- as.matrix(iris[,1:4])

# Test 1: Does iris return appropriate output?
test_that (
  "Ostats_multivariate returns expected output",
  {
    result1 <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE, hypervolume_args = list(method = 'box'), hypervolume_set_args = list(num.points.max = 1000))
    result1 <- as.numeric(result1$overlaps_norm)
    expected1 <- c(0.755, 0.783)
    expect_equal(result1, expected1, tolerance = 0.01)
  }
)

# Test 2: run multivariate data on single trait
# expect_equivalent ignores names.
test_that (
  "Ostats_multivariate returns the same output as Ostats if one trait is selected",
  {
    result2 <- Ostats_multivariate(traits = iris_traits[,1,drop=FALSE], plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE)
    expected2 <- Ostats(traits = iris_traits[,1,drop=FALSE], plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE)
    expect_equivalent(result2$overlaps_norm, expected2$overlaps_norm, tolerance = 0.01)
  }
)

# Test 3: Are the alternate arguments to hypervolume and hypervolume_set passed through correctly?
test_that (
  "Ostats_multivariate correctly handles arguments to hypervolume functions",
  {
    result3 <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = FALSE,
                                   hypervolume_args = list(method = 'box'), hypervolume_set_args = list(distance.factor = 2, num.points.max = 1000))
    result3 <- as.numeric(result3$overlaps_norm)
    expected3 <- c(0.93, 0.95)
    expect_equal(result3, expected3, tolerance = 0.01)
  }
)

# Test 4: Does null model return identical output if identical random seed is set?
test_that (
  "Ostats_multivariate null model returns identical output if identical random seed is set",
  {
    result4 <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = TRUE, nperm = 10, hypervolume_args = list(method = 'box'), hypervolume_set_args = list(num.points.max = 1000))
    result4b <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, run_null_model = TRUE, nperm = 10, hypervolume_args = list(method = 'box'), hypervolume_set_args = list(num.points.max = 1000))
    expect_equivalent(result4, result4b)
  }
)

# Test 5: Does Ostats_multivariate_plot return an object of class Ostats_plot_object when there are >2 traits?
test_that (
  "Ostats_multivariate_plot returns an object of class Ostats_plot_object when there are >2 traits",
  {
    plot5 <- Ostats_multivariate_plot(plots = rep(1, nrow(iris)), sp = iris$Species, traits = iris_traits)
    expect_s3_class(plot5, 'Ostats_plot_object')
  }
)

# Test 6: Does Ostats_multivariate_plot return an object of class Ostats_plot_object when there are exactly 2 traits?
test_that (
  "Ostats_multivariate_plot returns an object of class Ostats_plot_object when there are exactly 2 traits",
  {
    plot6 <- Ostats_multivariate_plot(plots = rep(1, nrow(iris)), sp = iris$Species, traits = iris_traits[, 1:2])
    expect_s3_class(plot6, 'Ostats_plot_object')
  }
)
