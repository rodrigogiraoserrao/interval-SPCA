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

# Computes the symbolic covariance matrix.
cov.k <- function(sigma.cc, mu.r, sigma.rr, k = 1) {
    .is.legal.k(k)
    
    delta <- delta.k(k)
    mu.r <- to.column(mu.r)
    e.rr <- sigma.rr + mu.r%*%t(mu.r)
    if (is.full.k(k)) sigma.cc + delta*e.rr
    else sigma.cc + delta*diag.matrix(e.rr)
}

# Estimates the symbolic covariance matrix from some interval-valued observations.
estimate.cov.k <- function(C, R, k = 1) {
    .is.legal.k(k)
    sigma.cc <- cov(C)
    mu.r <- to.column(colMeans(R))
    sigma.rr <- cov(R)
    cov.k(sigma.cc, mu.r, sigma.rr, k)
}

# Estimates the symbolic covariance matrix in a robust way.
robust.cov.k <- function(C, R, k = 1, ...) {
    .is.legal.k(k)
    r <- robustbase::covMcd(cbind(C, R), ...)
    p <- ncol(C)
    r$mu.c <- to.column(r$center[1:p])
    r$sigma.cc <- r$cov[1:p, 1:p]
    r$mu.r <- to.column(r$center[(p+1):(2*p)])
    r$sigma.rr <- r$cov[(p+1):(2*p), (p+1):(2*p)]
    r$e.rr <- r$sigma.rr + r$mu.r%*%t(r$mu.r)
    
    if (is.full.k(k)) r$cov.k <- r$sigma.cc + delta.k(k)*r$e.rr
    else r$cov.k <- r$sigma.cc + delta.k(k)*diag.matrix(r$e.rr)
    return(r)
}
