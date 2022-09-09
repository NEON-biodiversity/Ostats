library(Ostats)
context("Ostats_plot")

# Tests verify that objects of the correct class are created
# for continuous data, discrete data, and circular data.
# If a single trait is plotted, test inherits(p, "Ostats_plot_object")
# If multiple traits are plotted, test length(p) == number of traits
# and all(sapply(p, inherits, "Ostats_plot_object"))

# suppressWarnings() must be used to suppress a ggplot2 warning
# that arises when setting limits on geom_histogram().

# Set up data for test, making multiple trait columns from the ant data.
data(ant_data)

ant_data$time2 <- ant_data$time %% 12
ant_data$time3 <- ant_data$time %% 6

sp <- ant_data$species
plots <- ant_data$chamber
traits <- ant_data[, c('time', 'time2', 'time3')]

# Tests 1-4: single trait plots
test_that(
  "Single continuous trait plot with means panels returns valid output",
  {
    p1 <- Ostats_plot(plots = plots, sp = sp, traits = traits[, 1, drop = FALSE], discrete = FALSE, normalize = FALSE, means = TRUE)
    expect_true(inherits(p1, "Ostats_plot_object"))
  }
)

test_that(
  "Single discrete trait plot returns valid output",
  {
    p2 <- Ostats_plot(plots = plots, sp = sp, traits = traits[, 1, drop = FALSE], discrete = TRUE, normalize = FALSE)
    expect_true(inherits(p2, "Ostats_plot_object"))
  }
)

test_that(
  "Single continuous circular trait plot with means panels returns valid output",
  {
    p3 <- Ostats_plot(plots = plots, sp = sp, traits = traits[, 1, drop = FALSE], discrete = FALSE, normalize = FALSE, circular = TRUE, means = TRUE)
    expect_true(inherits(p3, "Ostats_plot_object"))
  }
)

test_that(
  "Single discrete circular trait plot returns valid output",
  {
    p4 <- Ostats_plot(plots = plots, sp = sp, traits = traits[, 1, drop = FALSE], discrete = TRUE, normalize = FALSE, circular = TRUE)
    expect_true(inherits(p4, "Ostats_plot_object"))
  }
)

# Tests 5-8: multiple trait plots
test_that(
  "Multiple continuous trait plot with means panels returns valid output",
  {
    plist1 <- Ostats_plot(plots = plots, sp = sp, traits = traits, discrete = FALSE, normalize = FALSE, means = TRUE)
    expect_true(length(plist1) == 3L & all(sapply(plist1, inherits, "Ostats_plot_object")))
  }
)

test_that(
  "Multiple discrete trait plot returns valid output",
  {
    plist2 <- Ostats_plot(plots = plots, sp = sp, traits = traits, discrete = TRUE, normalize = FALSE)
    expect_true(length(plist2) == 3L & all(sapply(plist2, inherits, "Ostats_plot_object")))
  }
)

test_that(
  "Multiple continuous circular trait plot with means panels returns valid output",
  {
    plist3 <- Ostats_plot(plots = plots, sp = sp, traits = traits, discrete = FALSE, normalize = FALSE, circular = TRUE, means = TRUE)
    expect_true(length(plist3) == 3L & all(sapply(plist3, inherits, "Ostats_plot_object")))
  }
)

test_that(
  "Multiple discrete circular trait plot returns valid output",
  {
    plist4 <- Ostats_plot(plots = plots, sp = sp, traits = traits, discrete = TRUE, normalize = FALSE, circular = TRUE)
    expect_true(length(plist4) == 3L & all(sapply(plist4, inherits, "Ostats_plot_object")))
  }
)
