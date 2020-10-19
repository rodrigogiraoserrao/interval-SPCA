# This script shows an example of symbolic PCA applied to the RTT data.
# For this example, we wish to find uncorrelated PCs and a diagonal model and we compare the three interval algebras.

# author: Rodrigo Girão Serrão
# email: rodrigogiraoserrao@gmail.com
# date: 15/09/2020

source("./sym.R")
source("./sym_pca.R")
source("./utils.R")
source("./RTTutils.R")

library(tidyverse)


# Load and prepare the data -----------------------------------------------

load("data/dataRTT_all.RData")

data <- data.T1.all.noNA %>% rtt.to.symbolic  # dataset to use
# Delete heavy variables from memory.
remove(
    list=c(
        paste0("data.T", 1:4, ".all"),
        paste0("data.T", 1:4, ".all.noNA")
    )
)



# Global variables --------------------------------------------------------

# Plot blue rectangles from regular observations and red rectangles for anomalies.
colours <- c("chartreuse", paste0("indianred", 1:4))
clevels <- levels(data$relayname)

makePlots <- TRUE
saveToPDF <- "scripts/mpca_imgs/"      # on saving plots: https://www.stat.berkeley.edu/~s133/saving.html
conv <- c(1, 2)  # plot conv[1]th VS conv[2]th conventional variables
pc <- c(1, 2)    # plot pc[1]th VS pc[2]th principal components
# Symbolic global variables
K <- 2
METHOD <- "moore"
RESTRICTION <- "orthogonal"
NORMALIZE <- FALSE

C <- as.matrix(data %>% select(!!paste0("C", 1:12)))
colnames(C) <- probes
R <- as.matrix(data %>% select(!!paste0("R", 1:12)))
colnames(R) <- probes


# Plot data without using PCA ---------------------------------------------

# Find the 2 columns with greatest variance and plot their intervals.
Cov <- estimate.cov.k(C, R, K)
vars <- diag(Cov)
trace <- sum(vars)

if (makePlots) {
    # Sort the conventional variables by variance.
    conv_ <- sort.int(vars, decreasing = TRUE, index.return = TRUE)$ix[conv]
    plot.VS(C, R, conv_, border = colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
    # Now plot the intervals of the variable conv[1].
    plot.interval(C, R, conv_[1], colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
    # Now plot the intervals of the variable conv[2].
    plot.interval(C, R, conv_[2], colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
}
print("Variance explained pre-PCA")
print(unname(100*cumsum(sort(vars, decreasing = TRUE))/sum(vars)))


# Apply PCA with the 3 interval algebras ----------------------------------

# Call the symbolic PCAs here.
result <- sym.pca(C, R, K, interval.algebra = METHOD, restriction = RESTRICTION, normalize = NORMALIZE)
rownames(result$vectors) <- probes
newC <- result$score.C(C)
colnames(newC) <- paste0(METHOD, "_pc_", 1:ncol(newC))
newR <- result$score.R(R)
colnames(newR) <- paste0(METHOD, "_pc_", 1:ncol(newR))
newCov <- estimate.cov.k(newC, newR, K)
newVars <- diag(newCov)
newTrace <- sum(newVars)
print(paste("abs(newTrace - trace)/trace = ", abs(newTrace - trace)/trace))
print(unname(round(100*cumsum(newVars)/newTrace, 2)))
absR <- abs(newR)

if (makePlots) {
    # Plot the rectangles for the PCs in pc
    plot.VS(newC, newR, pc, border = colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
    # Now only plot the pc[1] with y depending on the range of the interval.
    plot.interval(newC, newR, pc[1], col = colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
    # Now only plot the pc[2].
    plot.interval(newC, newR, pc[2], col = colours[match(data$relayname, clevels)], saveToPDF = saveToPDF)
}