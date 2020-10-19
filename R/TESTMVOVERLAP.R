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
