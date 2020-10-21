# Make some fake data to test our hypervolume overlap functions.

# Test on the classic iris dataset

hypervolume::hypervolume(data = subset(iris, Species %in% 'setosa')[,1:4], method = 'gaussian', verbose = TRUE)

iris_spp <- split(iris, iris$Species)

# Use MV
set.seed(1111)
pairwise_overlap_mv(iris_spp$setosa[,1:4], iris_spp$versicolor[,1:4])

set.seed(1112)
pairwise_overlap_mv(iris_spp$setosa[,1:4], iris_spp$virginica[,1:4])

set.seed(1113)
pairwise_overlap_mv(iris_spp$versicolor[,1:4], iris_spp$virginica[,1:4])

# Use newly created pairwise_overlap function that determines whether input is univariate or multivariate
set.seed(1111)
pairwise_overlap(iris_spp$setosa[,1:4], iris_spp$versicolor[,1:4])

set.seed(1112)
pairwise_overlap(iris_spp$setosa[,1:4], iris_spp$virginica[,1:4])

set.seed(1113)
pairwise_overlap(iris_spp$versicolor[,1:4], iris_spp$virginica[,1:4])

#### test Ostats_multivariate()!

iris_traits <- iris[,1:4]
iris_sp <- iris[,5]
iris_plots <- rep('foobar', nrow(iris))

O1 <- Ostats_multivariate(traits = as.matrix(iris_traits), sp = as.factor(iris_sp), plots = as.factor(iris_plots), nperm = 5, hypervolume_args = list(verbose = FALSE, method = 'box'))

traits = as.matrix(iris_traits)
plots = as.factor(iris_plots)
sp = as.factor(iris_sp)
output = "median"
weight_type= "hmean"
nperm = 99
nullqs = c(0.025, 0.975)
shuffle_weights = FALSE
swap_means = FALSE
hypervolume_args = list()

# Comparison with the non multivariate version
O1separate <- Ostats(traits = as.matrix(iris_traits), sp = as.factor(iris_sp), plots = as.factor(iris_plots), nperm = 5)
