##
#
# This module implements Principal Component Analysis (PCA) for interval valued variables.
#
##

source("./utils.R")
source("./sym.R")

# Computes the symbolic principal components.
# The columns in the centers and ranges input matrices denote different variables and rows define measurements.
# Returns a list with $values and $vectors, much like the eigen() function.
sym.pca <- function(C, R, k,
                    interval.algebra = c("vector space", "extended algebra", "moore algebra"),
                    restriction = c("orthogonal", "uncorrelated")
) {
    .is.legal.k(k)
    if (nrow(C) != nrow(R)) stop("C and R matrices have different numbers of rows.")
    if (ncol(C) != ncol(R)) stop("C and R matrices have different numbers of columns.")
    interval.algebra <- match.arg(interval.algebra)
    restriction <- match.arg(restriction)
    
    p <- ncol(C)
    n <- nrow(C)
    
    result <- list(
        values = rep(0, p),
        vectors = matrix(0, nrow = p, ncol = p),
        interval.algebra = interval.algebra,
        restriction = restriction,
        k = k
    )
    
    if (interval.algebra == "vector space") {
        
        sigma.k <- cov.k(C, R, k)
        partial.result <- eigen(sigma.k, symmetric = TRUE)
        result$values <- partial.result$values
        result$vectors <- partial.result$vectors
        
        if (restriction == "orthogonal" || (restriction == "uncorrelated" && is.full.k(k))) {
            return(result)
            
        } else if (restriction == "uncorrelated" && is.diagonal.k(k)) {
            sigma.cc <- cov(C)
            for (i in 2:p) {
                # Build the projection matrix out of the previous PCs
                orthogonalized <- gram.schmidt(result$vectors[, 1:(i-1)])
                projection.matrix <- orthogonal.projection.matrix(orthogonalized)
                modified.sigma.k <- t(projection.matrix)%*%sigma.k%*%projection.matrix
                partial.result <- eigen(modified.sigma.k)
                result$values[i] <- partial.result$values[1]
                result$vectors[, i] <- partial.result$vectors[, 1]
            }
            return(result)
        }
        
    } else {
        stop(paste("interval.algebra", interval.algebra, "not implemented yet."))
    }
}