context("pairwise_overlap")
#Ostats example
#how do I make it into a test?
#attempt with pairwise_overlap

# Test 1: pairwise overlap for two distributions

n <- 1e7
set.seed(111)
a <- runif(n, min = 0, max = 2)
b <- runif(n, min = 1, max = 3)

result1 <- pairwise_overlap(a, b)
expected1 <- c(overlap = 0.5, overlap_a = 0.5, overlap_b = 0.5)

test_that (
"pairwise overlap returns expected output",
{
  expect_equal(result1, expected1, tolerance = 0.01)
}
)
