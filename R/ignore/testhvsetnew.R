# Benchmark original method versus new method for hypervolume overlap estimation

traits = as.matrix(iris[, 1:4])
sp = iris[,5]
plots = factor(rep(1, nrow(iris)))

# Clean input, removing missing values and species with <2 values.
sp <- as.character(sp)
dat <- cbind(as.data.frame(traits), sp = sp)
dat <- dat[stats::complete.cases(dat), ]
abunds <- table(dat$sp)
abunds <- abunds[abunds>1]
dat <- dat[dat$sp %in% names(abunds), ]
traitlist <- split(dat[, -ncol(dat)], dat$sp)
uniquespp <- unique(dat$sp)
nspp <- length(uniquespp)

# Define common grid limits so that all density functions and hypervolumes are estimated across the same domain.
grid_limits <- apply(dat[, -ncol(dat), drop = FALSE], 2, range)

# Add a multiplicative factor to the upper and lower end of the range so that the tails aren't cut off.
extend_grid <- c(-0.5, 0.5) %*% t(apply(grid_limits,2,diff))
grid_limits <- grid_limits + extend_grid

combs = utils::combn(1:nspp, 2)

overlaps = NULL

density_list <- lapply(traitlist, function(x) Ostats:::trait_density(x, grid_limits, normal = TRUE, discrete = FALSE, circular = FALSE, unique_values = NULL, density_args = list(), circular_args = list()))

for (idx in 1:ncol(combs)) {

  overlaps <- c(overlaps, Ostats:::pairwise_overlap(density_list[[combs[1, idx]]], density_list[[combs[2, idx]]], discrete = FALSE, density_args = list(), hypervolume_set_args = list()))

}

density_joined <- do.call(hypervolume::hypervolume_join, density_list)
overlaps2 <- hypervolume::hypervolume_set_n_intersection(density_joined)
