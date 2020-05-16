##
#
# This module implements Principal Component Analysis (PCA) for interval valued variables.
#
##

source("./utils.R")
source("./sym.R")

# Computes the symbolic principal components.
# The columns in the centers and ranges input matrices denote different variables and rows define measurements.
sym.pca <- function(C, R, k,
                    interval.algebra = c("vector space", "extended algebra", "moore algebra"),
                    restriction = c("orthogonal", "uncorrelated")) {
    .is.legal.k(k)
    if (nrow(C) != nrow(R)) stop("C and R matrices have different numbers of rows.")
    if (ncol(C) != ncol(R)) stop("C and R matrices have different numbers of columns.")
    interval.algebra <- match.arg(interval.algebra)
    restriction <- match.arg(restriction)
    
    if (interval.algebra == "vector space") {
        
        sigma.k <- cov.k(C, R, k)
        if (restriction == "orthogonal" || (is.full.k(k) && restriction = "uncorrelated")) {
            eigen(sigma.k)
        }
        
    } else {
        stop(paste("interval.algebra", interval.algebra, "not implemented yet."))
    }
}