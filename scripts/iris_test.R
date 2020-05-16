##
#
# Test the Symbolic PCA for intervals with the iris data.
#
##

source("./sym_pca.R")

library(ggplot2)

data(iris)
set.seed(73) # for reproducibility

## Set symbolic variables here

K <- 2                       # the symbolic model to be used.
algebra <- "vector space"    # the interval algebra to be used.
restriction <- "orthogonal"  # the type of restriction to enforce.

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

# Prepare plotting area
par(mfrow=c(2,2))

# Find the 2 columns with greatest variance and plot their intervals.
sigma.k <- cov.k(C, R, K)
vars.k <- diag(sigma.k)
cols <- sort.int(vars.k, decreasing = TRUE, index.return = TRUE)$ix[1:2]
print("The columns with greatest variance are")
print(colnames(C)[cols])
plot(
    c(
        min(C[, cols[1]]-R[, cols[1]])-2,
        max(C[, cols[1]]+R[, cols[1]])+2
    ),
    c(
        min(C[, cols[2]]-R[, cols[2]])-2,
        max(C[, cols[2]]+R[, cols[2]])+2
    ),
    xlab = colnames(C)[cols[1]],
    ylab = colnames(C)[cols[2]],
    main = paste("Intervals for", colnames(C)[cols[1]], "and", colnames(C)[cols[2]]),
)
rect(
    xleft = C[, cols[1]]-R[, cols[1]],
    ybottom = C[, cols[2]]-R[, cols[2]],
    xright = C[, cols[1]]+R[, cols[1]],
    ytop = C[, cols[2]]+R[, cols[2]],
    border = c("red", "blue", "green")[match(species, levels(species))]
)
# Now plot the intervals of the variable with greatest variance.
plot(
    c(
        min(C[, cols[1]]-R[, cols[1]])-2,
        max(C[, cols[1]]+R[, cols[1]])+2
    ),
    c(
        min(R[, cols[1]])-2,
        max(R[, cols[1]])+2
    ),
    xlab = colnames(C)[cols[1]],
    ylab = "Range",
    main = paste("Intervals for", colnames(C)[cols[1]], "with y depending on the range"),
)
rect(
    xleft = C[, cols[1]]-R[, cols[1]],
    ybottom = R[, cols[1]]-0.04,
    xright = C[, cols[1]]+R[, cols[1]],
    ytop = R[, cols[1]]+0.04,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)

## Call the symbolic PCAs here.

result <- sym.pca(C, R, K, interval.algebra = "vector space", restriction = "orthogonal")
newC <- C %*% result$vectors
newR <- R %*% result$vectors
absR <- abs(newR)
# Plot the rectangles for the first two PCAs
plot(
    c(
        min(newC[,1]-absR[,1])-2,
        max(newC[,1]+absR[,1])+2
    ),
    c(
        min(newC[,2]-absR[,2])-2,
        max(newC[,2]+absR[,2])+2
    ),
    xlab = "1st PCA",
    ylab = "2nd PCA",
    main = "Intervals of the 1st and 2nd PCAs.",
)
rect(
    xleft = newC[, 1]-newR[, 1],
    ybottom = newC[, 2]-newR[, 2],
    xright = newC[, 1]+newR[, 1],
    ytop = newC[, 2]+newR[, 2],
    border = c("red", "blue", "green")[match(species, levels(species))]
)
# Now only plot the 1st PCA with y depending on the range of the interval.
plot(
    c(
        min(newC[,1]-absR[,1])-2,
        max(newC[,1]+absR[,1])+2
    ),
    c(
        min(newR[,1]),
        max(newR[,1])
    ),
    xlab = "Intervals",
    ylab = "Ranges",
    main = "Intervals of the 1st PCA with y depending only on the range.",
)
rect(
    xleft = newC[, 1]-newR[, 1],
    ybottom = newR[, 1]-0.02,
    xright = newC[, 1]+newR[, 1],
    ytop = newR[, 1]+0.02,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)
