##
#
# Test the Symbolic PCA for intervals with the iris data.
#
##

source("./sym_pca.R")

data(iris)
set.seed(73) # for reproducibility

## Start by grouping the data inside each species to create some fake intervals.

# Create intervals out of approximately `grouping_size` observations.
grouping_size <- 5
# Group observations with the intervals indexed by this column.
interval_grouping <- iris[,5] != iris[,5]

for (factor in levels(iris[,5])) {
    mask <- iris[,5] == factor
    how_many <- sum(mask)
    groups <- floor(how_many/grouping_size)
    # Force all groupings to be non-empty.
    interval_grouping[mask] <- max(interval_grouping) + c(1:groups, sample(1:(groups), how_many-groups, replace = TRUE))
}

n_groups <- max(interval_grouping)
C <- matrix(0, nrow = n_groups, ncol = 4)
R <- matrix(0, nrow = n_groups, ncol = 4)
colnames(C) <- colnames(R) <- colnames(iris[, c(1, 2, 3, 4)])
species <- factor(character(length = n_groups), levels(iris[, 5]))
# Create the intervals
colMax <- function(data) { sapply(data, max, na.rm = TRUE) }
colMin <- function(data) { sapply(data, min, na.rm = TRUE) }
for (i in 1:n_groups) {
    mask <- interval_grouping == i
    subset <- iris[mask, c(1, 2, 3, 4)]
    subset_species <- iris[mask, 5]
    mins <- colMin(subset)
    maxs <- colMax(subset)
    C[i, ] <- (maxs + mins)/2
    R[i, ] <- maxs - mins
    if (length(unique(subset_species)) > 1) stop("There's something wrong.")
    species[i] <- subset_species[1]
}