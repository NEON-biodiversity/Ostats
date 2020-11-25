#vignette tests Isadora learning#

##
#example of how to use the testthat:
result1 <- Ostats_multivariate(traits = iris_traits, plots = factor(rep(c('a','b'),times=75)), sp = iris$Species, random_seed = 111, nperm = 1)

result1 <- as.numeric(result1$overlaps_norm)
expected1 <- c(0.62, 0.72)

test_that (
  "Ostats_multivariate returns expected output", #description of the test
  {
    expect_equal(result1, expected1, tolerance = 0.01) # what the function does and what it should be doing and I expect it to be equal
  }
)

#if you run this and nothing happens it means its true
##

library(testthat)
?test_that

###line 111#####
# Ostats handles the different density arguments correctly?
# Example of specifying the arguments for making density estimates
Ostats_example2 <- Ostats(traits = as.matrix(dat[,'log_weight', drop = FALSE]),
                          sp = factor(dat$taxonID),
                          plots = factor(dat$siteID),
                          density_args=list(bw = 'nrd0', adjust = 2, n=200))



result1<-as.numeric(Ostats_example2$overlaps_norm)
expected1<-c(0.89498411,0.01749386)

test_that (
  "Ostats_Ostats_example2_ handles the different density arguments correctly",
  {
    expect_equal(result1, expected1)
  }
)

##I cant get what this function does. It compares if the values are equal? or it can see the whole function behind it? If yes, how if the code is outside the function test_that?
##How do I know if I tested the different density arguments?

####line162####
#the hourly circular data is handled correctly?
circular_example <- Ostats(traits = as.matrix(ant_data[, 'time', drop = FALSE]),
                           sp = factor(ant_data$species),
                           plots = factor(ant_data$chamber),
                           data_type = "circular_discrete")


result2<-as.numeric(circular_example$overlaps_norm)
expected2<-c(0.6833794, 0.6363885)

test_that (
  "the hourly circular data is handled correctly",
  {
    expect_equal(result2, expected2, tolerance = 0.01)
  }
)

#if I dont specify the tolerance, it gets an error saying that the values are not the same. why?

#####line 217#####
Ostats_reg_example <- Ostats_regional(traits = as.matrix(HARV2018[, 'log_weight', drop = FALSE]),
                                      plots = factor(HARV2018$siteID),
                                      sp = factor(HARV2018$taxonID),
                                      reg_pool_traits = reg_pool_traits,
                                      reg_pool_sp = reg_pool_sp)

result3<-as.numeric(Ostats_reg_example$overlaps_reg)
expected3<-0.914048

test_that (
  "Ostats_reg_example",
  {
    expect_equal(result3, expected3, tolerance = 0.01)
  }
)

#without the tolerance there is an error. Thats because the numbers in the results are far more extensive that those displayed?
