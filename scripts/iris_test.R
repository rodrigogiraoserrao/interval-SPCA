##
#
# Test the Symbolic PCA for intervals with the iris data.
# Plot intervals for two variables and for two PCs, as well as the intervals for those variables alone.
#
##

source("./sym_pca.R")

library(ggplot2)

data(iris)
set.seed(73) # for reproducibility

## Plot the conventional variables conv[1] and conv[2] when ordered by variance.
## Plot the PCs pc[1] and pc[2] when ordered by variance.
conv <- c(1, 3)
pc <- c(1, 3)

## Set symbolic variables here
K <- 5                         # the symbolic model to be used.
algebra <- "extended algebra"      # the interval algebra to be used.
restriction <- "orthogonal"  # the type of restriction to enforce.

# Helper function to set the axis of the plots
stretch.axis <- function(axis) {
    axis.range <- axis[2] - axis[1]
    delta <- 0.1*axis.range
    c(axis[1]-delta, axis[2]+delta)
}

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
par(mfrow = c(2,3))

# Find the 2 columns with greatest variance and plot their intervals.
sigma.k <- cov.k(C, R, K)
vars.k <- diag(sigma.k)
trace <- sum(vars.k)
cols <- sort.int(vars.k, decreasing = TRUE, index.return = TRUE)$ix[conv]
xaxis <- stretch.axis(c(
    min(C[, cols[1]]-R[, cols[1]]),
    max(C[, cols[1]]+R[, cols[1]])
))
yaxis <- stretch.axis(c(
    min(C[, cols[2]]-R[, cols[2]]),
    max(C[, cols[2]]+R[, cols[2]])
))
plot(
    xaxis,
    yaxis,
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
# Now plot the intervals of the variable conv[1].
plot_height <- max(R[, cols[1]]) - min(R[, cols[1]])
plot(
    xaxis,
    c(
        min(R[, cols[1]]),
        max(R[, cols[1]])
    ),
    xlab = colnames(C)[cols[1]],
    ylab = "Range",
    main = paste("Intervals for", colnames(C)[cols[1]], "with y = range"),
)
rect(
    xleft = C[, cols[1]]-R[, cols[1]],
    ybottom = R[, cols[1]]-plot_height/100,
    xright = C[, cols[1]]+R[, cols[1]],
    ytop = R[, cols[1]]+plot_height/100,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)
# Now plot the intervals of the variable conv[2].
plot_height <- max(R[, cols[2]]) - min(R[, cols[2]])
plot(
    yaxis,
    c(
        min(R[, cols[2]]),
        max(R[, cols[2]])
    ),
    xlab = colnames(C)[cols[2]],
    ylab = "Range",
    main = paste("Intervals for", colnames(C)[cols[2]], "with y = range"),
)
rect(
    xleft = C[, cols[2]]-R[, cols[2]],
    ybottom = R[, cols[2]]-plot_height/100,
    xright = C[, cols[2]]+R[, cols[2]],
    ytop = R[, cols[2]]+plot_height/100,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)

## Call the symbolic PCAs here.

result <- sym.pca(C, R, K, interval.algebra = algebra, restriction = restriction)
newC <- C %*% result$vectors
newR <- R %*% result$vectors
newVars <- var.k(newC, newR, K)
newTrace <- sum(newVars)
print(paste("abs(newTrace - trace) = ", abs(newTrace - trace)))
print(round(100*cumsum(newVars)/newTrace, 2))
absR <- abs(newR)
# Plot the rectangles for the PCs in pc
xaxis <- stretch.axis(c(
    min(newC[,pc[1]]-absR[,pc[1]]),
    max(newC[,pc[1]]+absR[,pc[1]])
))
yaxis <- stretch.axis(c(
    min(newC[,pc[2]]-absR[,pc[2]]),
    max(newC[,pc[2]]+absR[,pc[2]])
))
plot(
    xaxis,
    yaxis,
    xlab = "pc[1]",
    ylab = "pc[2]",
    main = paste("Intervals for PCs", pc[1], "and", pc[2]),
)
rect(
    xleft = newC[, pc[1]]-newR[, pc[1]],
    ybottom = newC[, pc[2]]-newR[, pc[2]],
    xright = newC[, pc[1]]+newR[, pc[1]],
    ytop = newC[, pc[2]]+newR[, pc[2]],
    border = c("red", "blue", "green")[match(species, levels(species))]
)
# Now only plot the pc[1] with y depending on the range of the interval.
plot_height <- max(newR[, pc[1]]) - min(newR[, pc[1]])
plot(
    xaxis,
    c(
        min(newR[,pc[1]]),
        max(newR[,pc[1]])
    ),
    xlab = "Intervals",
    ylab = "Ranges",
    main = paste("Intervals of PC", pc[1], "with y = range."),
)
rect(
    xleft = newC[, pc[1]]-newR[, pc[1]],
    ybottom = newR[, pc[1]]-plot_height/100,
    xright = newC[, pc[1]]+newR[, pc[1]],
    ytop = newR[, pc[1]]+plot_height/100,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)
# Now only plot the pc[2] with y = range.
plot_height <- max(newR[, pc[2]]) - min(newR[, pc[2]])
plot(
    yaxis,
    c(
        min(newR[,pc[2]]),
        max(newR[,pc[2]])
    ),
    xlab = "Intervals",
    ylab = "Ranges",
    main = paste("Intervals of PC", pc[2], "with y = range."),
)
rect(
    xleft = newC[, pc[2]]-newR[, pc[2]],
    ybottom = newR[, pc[2]]-plot_height/100,
    xright = newC[, pc[2]]+newR[, pc[2]],
    ytop = newR[, pc[2]]+plot_height/100,
    density = 50,
    col = c("red", "blue", "green")[match(species, levels(species))]
)