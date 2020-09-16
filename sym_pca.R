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
                    interval.algebra = c("vector space", "extended", "moore"),
                    restriction = c("orthogonal", "uncorrelated"),
                    normalize = FALSE
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
        restriction = restriction
    )
    
    if (!isTRUE(all.equal(as.vector(colMeans(C)), rep(0, p)))) {
        warning("C is not centred")
        result$muc <- colMeans(C)
        C <- t(t(C) - result$muc)
    } else {
        result$muc <- rep(0, p)
    }
    
    if (normalize) {
        l <- normalize.k(C, R, k)
        C <- l$C
        R <- l$R
    }
    
    if (interval.algebra == "vector space" && restriction == "orthogonal") {
        result <- .sym.pca.conventional.case(result, cov.k(C, R, k))
        
    } else if (interval.algebra == "vector space" && restriction == "uncorrelated" && is.full.k(k)) {
        result <- .sym.pca.conventional.case(result, cov.k(C, R, k))
        
    } else if (interval.algebra == "vector space" && restriction == "uncorrelated" && is.diagonal.k(k)) {
        result <- .sym.pca.projection.case(
            result,
            conventional.matrix = cov.k(C, R, k),
            orthogonality.matrix = cov(C)
        )
        
    } else if (interval.algebra == "extended" && restriction == "orthogonal") {
        result <- .sym.pca.conventional.case(result, cov.k(C, R, k))
        
    } else if (interval.algebra == "extended" && restriction == "uncorrelated" && is.diagonal.k(k)) {
        result <- .sym.pca.projection.case(
            result,
            conventional.matrix = cov.k(C, R, k),
            orthogonality.matrix = cov(C)
        )
        
    } else if (interval.algebra == "extended" && restriction == "uncorrelated" && is.full.k(k)) {
        # this instance is special-cased, solve directly
        sigma.cc <- cov(C)
        quadratic.matrix <- sigma.cc + delta.k(k) * t(R)%*%R/nrow(R)
        for (i in safe.colon(1, p)) {
            modified.vectors <- sigma.cc%*%result$vectors[, safe.colon(1, i-1)] + delta.k(k) * (t(R)%*%abs(R%*%result$vectors[, safe.colon(1, i-1)]))/nrow(R)
            orthogonalized <- gram.schmidt(modified.vectors)
            projection.matrix <- orthogonal.projection.matrix(orthogonalized)
            modified.matrix <- t(projection.matrix) %*% quadratic.matrix %*% projection.matrix
            partial.result <- eigen(modified.matrix, symmetric = TRUE)
            result$values[i] <- partial.result$values[1]
            eigenvector <- partial.result$vectors[, 1]
            # we may need to modify the sign of `eigenvector` depending on the sign of `t(eigenvector)%*%R`
            # we only have observations of R so estimate it
            result$vectors[, i] <- sgn(colMeans(R)%*%eigenvector)[1]*eigenvector
        }
        
    } else if (interval.algebra == "moore" && restriction == "orthogonal") {
        result <- .sym.pca.modified.quadratic.by.absolute.values(
            result,
            conv.quadratic.matrix = cov(C),
            abs.quadratic.matrix = delta.k(k) * t(R)%*%R/nrow(R)
        )
        
    } else if (interval.algebra == "moore" && restriction == "uncorrelated" && is.diagonal.k(k)) {
        sigma.cc <- cov(C)
        result <- .sym.pca.modified.quadratic.by.absolute.values(
            result,
            conv.quadratic.matrix = sigma.cc,
            abs.quadratic.matrix = delta.k(k) * t(R)%*%R/nrow(R),
            orthogonality.matrix = sigma.cc
        )
        
    } else if (interval.algebra == "moore" && restriction == "uncorrelated" && is.full.k(k)) {
        # this instance is special-cased, solve directly
        sigma.cc <- cov(C)
        e.rr <- t(R)%*%R/nrow(R)
        if (!all.equal(e.rr, abs(e.rr))) stop("E[RR^T] matrix should be non-negative.")
        Ms <- list()
        Ds <- list()
        d.matrix <- d(p)
        for (s in safe.colon(1, 2^p)) {
            Ds[[s]] <- diag.matrix(d.matrix[, s])
            Ms[[s]] <- sigma.cc + delta.k(k) * Ds[[s]]%*%e.rr%*%Ds[[s]]
        }
        for (i in safe.colon(1, p)) {
            lambda.s <- rep(0, 2^p)
            alpha.s <- matrix(0, nrow = p, ncol = ncol(d.matrix))
            same.s <- rep(FALSE, 2^p)
            prev.vectors <- result$vectors[, safe.colon(1, i-1)]
            for (s in safe.colon(1, length(Ds))) {
                modified.vectors <- sigma.cc%*%prev.vectors + delta.k(k)*Ds[[s]]%*%e.rr%*%abs(prev.vectors)
                orthogonalized <- gram.schmidt(modified.vectors)
                if (!isTRUE(all.equal(
                    diag(ncol(modified.vectors)),
                    t(orthogonalized)%*%orthogonalized
                ))) stop("Gram-Scmidt failed.")
                projection.matrix <- orthogonal.projection.matrix(orthogonalized)
                quadratic.matrix <- t(projection.matrix)%*%Ms[[s]]%*%projection.matrix
                partial.result <- eigen(quadratic.matrix, symmetric = TRUE)
                lambda.s[s] <- partial.result$values[1]
                eig <- t(t(partial.result$vectors[, 1]))
                alpha.s[, s] <- eig
                same.s[s] <- isTRUE(all.equal(
                    t(eig)%*%sigma.cc%*%eig + delta.k(k)*t(eig)%*%Ds[[s]]%*%e.rr%*%Ds[[s]]%*%eig,
                    t(eig)%*%sigma.cc%*%eig + delta.k(k)*  abs(t(eig))   %*%e.rr%*%  abs(eig)
                ))
                if (!isTRUE(all.equal(eig, projection.matrix%*%eig))) stop("Projection failed.")
                if (!isTRUE(all.equal(
                    rep(0, i-1),
                    as.vector(t(eig)%*%(sigma.cc%*%prev.vectors + delta.k(k)*Ds[[s]]%*%e.rr%*%abs(prev.vectors)))
                ))) stop("Eigenvector isn't uncorrelated")
                # print(isTRUE(all.equal(
                #     t(eig)%*%sigma.cc%*%prev.vectors + delta.k(k)*t(eig)%*%Ds[[s]]%*%e.rr%*%abs(prev.vectors),
                #     t(eig)%*%sigma.cc%*%prev.vectors + delta.k(k)*abs(t(eig))%*%e.rr%*%abs(prev.vectors)
                # )))
                # if (eig[1]*d.matrix[1, s] < 0) eig <- -eig
                # print(isTRUE(all.equal(
                #     t(eig)%*%sigma.cc%*%prev.vectors + delta.k(k)*t(eig)%*%Ds[[s]]%*%e.rr%*%abs(prev.vectors),
                #     t(eig)%*%sigma.cc%*%prev.vectors + delta.k(k)*abs(t(eig))%*%e.rr%*%abs(prev.vectors)
                # )))
            }
            if (!any(same.s)) stop(paste("No subproblem for i =", i, "worked."))
            s.star <- which.max(lambda.s)
            result$values[i] <- lambda.s[s.star]
            # adjust the signs of `eigenvector` if needed
            eigenvector <- alpha.s[, s.star]
            print("============")
            # print(eigenvector)
            # print(sgn(eigenvector))
            # print(d.matrix[, s.star])
            # print(lambda.s[s.star])
            # print(t(eigenvector)%*%sigma.cc%*%eigenvector + delta.k(k)*(t(eigenvector)%*%Ds[[s.star]]%*%e.rr%*%Ds[[s.star]]%*%eigenvector))
            # print(t(eigenvector)%*%sigma.cc%*%eigenvector + delta.k(k)*(abs(t(eigenvector))%*%e.rr%*%abs(eigenvector)))
            print(t(eigenvector)%*%sigma.cc%*%prev.vectors + delta.k(k)*t(eigenvector)%*%Ds[[s.star]]%*%e.rr%*%abs(prev.vectors))
            
            #eigenvector <- diag.matrix(sgn(Ds[[s.star]]%*%eigenvector))%*%eigenvector
            if (eigenvector[1]*d.matrix[1, s.star] < 0) eigenvector <- -eigenvector
            
            print(t(eigenvector)%*%sigma.cc%*%prev.vectors + delta.k(k)*t(eigenvector)%*%Ds[[s.star]]%*%e.rr%*%abs(prev.vectors))
            print(t(eigenvector)%*%sigma.cc%*%prev.vectors + delta.k(k)*abs(t(eigenvector))%*%e.rr%*%abs(prev.vectors))
            
            result$vectors[, i] <- eigenvector
        }
        
    } else {
        stop(paste(interval.algebra, restriction, k, "not implemented."))
    }
    
    
    result$score.C <- function(C) {
        return(C %*% result$vectors)
    }
    result$score.R <- function(R) {
        if (interval.algebra == "vector space") return(R %*% result$vectors)
        if (interval.algebra ==     "extended") return(abs(R %*% result$vectors))
        if (interval.algebra ==        "moore") return(R %*% abs(result$vectors))
    }
    result$score <- function(C, R) {
        return(list(
            C = result$score.C(C),
            R = result$score.R(R)
        ))
    }
    return(result)
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

.sym.pca.modified.quadratic.by.absolute.values <- function(
    result, conv.quadratic.matrix, abs.quadratic.matrix, orthogonality.matrix = diag(ncol(result$vectors))
) {
    p <- ncol(result$vectors)
    Ms <- list()
    d.matrix <- d(p)
    # Create all the matrices modified by signs.
    for (s in safe.colon(1, 2^p)) {
        Ds <- diag.matrix(d.matrix[, s])
        Ms[[s]] <- conv.quadratic.matrix + Ds %*% abs.quadratic.matrix %*% Ds
    }
    
    alpha.s <- matrix(0, nrow = p, ncol = ncol(d.matrix))
    lambda.s <- rep(0, ncol(alpha.s))
    for (i in safe.colon(1, p)) {
        modified.vectors <- orthogonality.matrix %*% result$vectors[, safe.colon(1, i-1)]
        orthogonalized <- gram.schmidt(modified.vectors)
        projection.matrix <- orthogonal.projection.matrix(orthogonalized)
        # Go over each of the subproblems
        for (s in safe.colon(1, ncol(alpha.s))) {
            modified.matrix <- t(projection.matrix) %*% Ms[[s]] %*% projection.matrix
            partial.result <- eigen(modified.matrix, symmetric = TRUE)
            lambda.s[s] <- partial.result$values[1]
            alpha.s[, s] <- partial.result$vectors[, 1]
        }
        s.star <- which.max(lambda.s)
        result$values[i] <- lambda.s[s.star]
        result$vectors[, i] <- alpha.s[, s.star]
    }
    
    return(result)
}
