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
    
    if (interval.algebra == "vector space" && restriction == "orthogonal") {
        return(.sym.pca.conventional.case(result, cov.k(C, R, k)))
    } else if (interval.algebra == "vector space" && restriction == "uncorrelated" && is.full.k(k)) {
        return(.sym.pca.conventional.case(result, cov.k(C, R, k)))
    } else if (interval.algebra == "vector space" && restriction == "uncorrelated" && is.diagonal.k(k)) {
        return(.sym.pca.projection.case(
            result,
            conventional.matrix = cov.k(C, R, k),
            orthogonality.matrix = cov(C)
        ))
    } else if (interval.algebra == "extended algebra" && restriction == "orthogonal") {
        return(.sym.pca.conventional.case(result, cov.k(C, R, k)))
    } else if (interval.algebra == "extended algebra" && restriction == "uncorrelated" && is.diagonal.k(k)) {
        return(.sym.pca.projection.case(
            result,
            conventional.matrix = cov.k(C, R, k),
            orthogonality.matrix = cov(C)
        ))
    } else {
        stop(paste(interval.algebra, restriction, k, "not implemented yet."))
    }
}

.sym.pca.conventional.case <- function(result, conventional.matrix) {
    partial.result <- eigen(conventional.matrix, symmetric = TRUE)
    result$vectors <- partial.result$vectors
    result$values <- partial.result$values
    return(result)
}

.sym.pca.projection.case <- function(result, conventional.matrix, orthogonality.matrix = diag(ncol(result$vectors))) {
    
    # For the first vector the orthogonality restrictions are empty.
    partial.result <- eigen(conventional.matrix, symmetric = TRUE)
    result$values[1] <- partial.result$values[1]
    result$vectors[, 1] <- partial.result$vectors[, 1]
    
    # For each subsequent vector, find the projection matrix and solve the modified problem.
    for (i in safe.colon(2, ncol(conventional.matrix))) {
        # Build the orthogonal projection matrix.
        modified.vectors <- orthogonality.matrix %*% result$vectors[, 1:(i-1)]
        orthogonalized <- gram.schmidt(modified.vectors)
        projection.matrix <- orthogonal.projection.matrix(orthogonalized)
        # Modify the matrix inducing the quadratic form.
        modified.conventional.matrix <- t(projection.matrix) %*% conventional.matrix %*% projection.matrix
        # Find the largest eigenvalue (and corresp. eigenvector) and save it.
        partial.result <- eigen(modified.conventional.matrix, symmetric = TRUE)
        result$values[i] <- partial.result$values[1]
        result$vectors[, i] <- partial.result$vectors[, 1]
    }
    
    return(result)
}