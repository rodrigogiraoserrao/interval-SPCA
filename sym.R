##
#
# Define functions relevant to symbolic analysis of interval data.
#
##

# We need some of the symbolic utility functions.
source("./utils.R")

# Define all the symbolic methods by their \delta_k constant.
DELTAS <- c(0, 0.25, 1/12, 0.25, 1/12, 0.125, 1/24, 1/36 - dnorm(3)/(6*(2*dnorm(3)-1)))

# Checks if this integer k makes sense as the index for a symbolic model.
.is.legal.k <- function(k) {
    if (!as.integer(k) == k) stop("k should be an integer.")
    if (k < 0 || k > length(DELTAS)) stop(paste("k should be between 1 and", k))
}

# Checks if a given index k corresponds to a full symbolic model.
is.full.k <- function(k) {
    .is.legal.k(k)
    k %in% c(2, 3)
}

# Checks if a given index k corresponds to a diagonal symbolic model.
is.diagonal.k <- function(k) {
    .is.legal.k(k)
    k %in% c(1, 4, 5, 6, 7, 8)
}

# Retrieve the \delta_k constant relative to a symbolic model.
delta.k <- function(k = 1) {
    .is.legal.k(k)
    DELTAS[k]
}

# Computes the symbolic variance of a vector of intervals.
var.k <- function(C, R, k = 1) {
    .is.legal.k(k)
    diag(var(C)) + delta.k(k)*diag(t(R)%*%R)/nrow(R)
}

# Computes the covariance matrix for a vector of intervals.
cov.k <- function(C, R, k = 1) {
    .is.legal.k(k)
    sigma.cc <- cov(C)
    e.rr <- (t(R)%*%R)/nrow(R)
    delta <- delta.k(k)
    
    if (is.full.k(k)) sigma.cc + delta*e.rr
    else sigma.cc + delta*diag.matrix(e.rr)
}

# Normalize the symbolic variables to centre 0 and variance 1.
normalize.k <- function(C, R, k = 1) {
    .is.legal.k(k)
    # recenter
    mu.C <- colMeans(C)
    C_ <- sweep(C, 2, mu.C, "-")
    sd <- sqrt(var.k(C_, R, k))
    C_ <- sweep(C_, 2, sd, "/")
    R_ <- sweep(R, 2, sd, "/")
    return(list(C=C_, R=R_))
}