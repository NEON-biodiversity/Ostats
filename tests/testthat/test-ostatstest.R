context("pairwise_overlap")
library(tidyverse)
library(Ostats)
library(sfsmisc)
#Ostats example
#how do I make it into a test?
#attempt with pairwise_overlap
 a <- mtcars %>%
   filter(cyl==6)
 a <- a[,1]
 b <- mtcars %>%
   filter(cyl==4)
 b <- b[,1]
 result<- pairwise_overlap(a,b)
 result

 test_that (
   "pairwise overlap works",
   {
     expect_equal((pairwise_overlap(a,b)),result)
   }
 )
