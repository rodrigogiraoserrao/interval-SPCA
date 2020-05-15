##
#
# Define functions relevant to symbolic analysis of interval data.
#
##

# Define all the symbolic methods by their \delta_k constant.
DELTAS <- c(0, 0.25, 1/12, 0.25, 1/12, 0.125, 1/24, 1/36 - dnorm(3)/(6*(2*dnorm(3)-1)))

# Retrieve the \delta_k constant relative to a symbolic model.
delta.k <- function(k = 1) {
    if (k > length(DELTAS) || k <= 0) stop(paste("k must be between 1 and", length(DELTAS)))
    DELTAS[k]
}

# Computes the symbolic variance of a vector of intervals.
var.k <- function(C, R, k = 1) {}

# Computes the covariance matrix for a vector of intervals.
cov.k <- function(C, R, k = 1) {}