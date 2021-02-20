##
#
# This module implements Principal Component Analysis (PCA) for interval valued variables.
#
##

source("./utils.R")
source("./sym.R")

# Applies SPCA to the data given.
# Data can be a combination of the raw interval-valued observations and estimations of their variances/covariances.
# C and R are matrices where columns encode variables and rows encode observations.
# sigma.cc and sigma.rr are pxp square matrices where p=ncol(C)=ncol(R) (if C/R are present).
# mu.c and mu.r are px1 column vectors where p=ncol(C)=ncol(R) (if C/R are present).
# Using sigma.xx and mu.xx takes precedence over estimating parameters directly from the data in C or R.
# k, interval.algebra and restriction configure the SPCA setting.
spca <- function(
    interval.algebra = c("extended", "moore", "vector space"),
    restriction = c("orthogonal", "centre-uncorrelated"),
    model = c("diagonal", "full"),
    delta,
    normalise = FALSE,
    C = NULL,
    R = NULL,
    mu.c = NULL,
    sigma.cc = NULL,
    mu.r = NULL,
    sigma.rr = NULL
) {
    # INPUT VALIDATION - symbolic model checking.
    interval.algebra <- match.arg(interval.algebra)
    restriction <- match.arg(restriction)
    model <- match.arg(model)
    
    # INPUT VALIDATION - ensure the value for delta makes sense.
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")
    
    # Some of the following checks are redundant but are in place to try and the provide the most helpful error messages.
    # INPUT VALIDATION - dimension checking for C/R.
    if (!is.null(C)) {
        p.c <- ncol(C)
        n.c <- nrow(C)
    } else p.c <- n.c <- -1
    if (!is.null(R)) {
        p.r <- ncol(R)
        n.r <- nrow(R)
    } else p.r <- n.r <- -1
    if (p.c > -1 && p.r > -1) {
        if (!all.equal(dim(C), dim(R))) stop(
            paste0("C (", paste(dim(C), collapse = ', '), ") and R (", paste(dim(R), collapse = ', '), " have mismatched dims.")            
        )
    }
    
    # INPUT VALIDATION - check mu.c and mu.r are available or can be estimated.
    if (is.null(mu.c)) {
        if(is.null(C)) stop("Cannot infer E[C] without C.")
        mu.c <- colMeans(C)
    }
    mu.c <- to.column(mu.c)
    p.mu.c <- nrow(mu.c)
    if (is.null(mu.r)) {
        if(is.null(R)) stop("Cannot infer E[R] without R.")
        mu.r <- colMeans(R)
    }
    mu.r <- to.column(mu.r)
    p.mu.r <- nrow(mu.r)
    # Check the two mean column vectors have the same dim.
    if (!all.equal(dim(mu.c), dim(mu.r))) stop(paste0("E[C] (", p.mu.c, ") and E[R] (", p.mu.r, ") have mismatched dims."))
    
    # INPUT VALIDATION - check sigma.cc and sigma.rr are available or can be estimated.
    if (is.null(sigma.cc)) {
        if(is.null(C)) stop("Cannot infer Var(C) without C.")
        sigma.cc <- cov(C)
    }
    if (ncol(sigma.cc) != nrow(sigma.cc)) stop("Var(C) is not square.")
    p.sigma.cc <- ncol(sigma.cc)
    if (is.null(sigma.rr)) {
        if(is.null(R)) stop("Cannot infer Var(R) wihtout R.")
        sigma.rr <- cov(R)
    }
    if (ncol(sigma.rr) != nrow(sigma.rr)) stop("Var(R) is not square.")
    p.sigma.rr <- ncol(sigma.rr)
    # Check if the two covariance matrices have the same dim.
    if (!all.equal(dim(sigma.cc), dim(sigma.rr))) {
        paste0("Var(C) (", paste(dim(sigma.cc), collapse = ', '), ") and Var(R) (", paste(dim(sigma.rr), collapse = ', '), " have mismatched dims.")
    }
    
    # INPUT VALIDATION - there are inconsistencies in the value of `n` or `p`.
    ns <- c(n.c, n.r)
    if (!all((ns == max(ns)) | (-1 == ns))) stop("Number of observations cannot be inferred correctly.")
    ps <- c(p.c, p.r, p.mu.c, p.mu.r, p.sigma.cc, p.sigma.rr)
    if (!all((ps == max(ps)) | (-1 == ps))) stop("Number of variables cannot be inferred correctly.")
    p <- max(ps)
    
    e.rr <- sigma.rr + mu.r%*%t(mu.r)
    
    if (normalise) {
        A <- diag.matrix(cov.k(sigma.cc, e.rr, model, delta))^-.5
        sigma.cc <- A%*%sigma.cc%*%A
        e.rr <- A%*%e.rr%*%A
    }
    
    # Initialise the return variable.
    result <- list(
        values = rep(0, p),
        vectors = matrix(0, nrow = p, ncol = p),
        interval.algebra = interval.algebra, # Keep these values for reference. ↓
        restriction = restriction,
        model = model,
        delta = delta,
        normalise = normalise,
        sigma.cc = sigma.cc,
        e.rr = e.rr
    )
    
    if (model == "diagonal" && !is.null(C) && !is.null(R)) {
        warn.about.linear.dependences(C, R)
    }
    
    ### Extended Algebra
    if (interval.algebra == "extended" && restriction == "orthogonal") {
        result <- .spca.conventional.case(result, cov.k(sigma.cc, e.rr, model, delta))
        
    } else if (interval.algebra == "extended" && restriction == "centre-uncorrelated") {
        result <- .spca.projection.case(
            result,
            conventional.matrix = cov.k(sigma.cc, e.rr, model, delta),
            orthogonality.matrix = sigma.cc
        )
    
    ### Moore Algebra
    } else if (interval.algebra == "moore" && restriction == "orthogonal") {
        result <- .spca.modified.quadratic.by.absolute.values(
            result,
            conv.quadratic.matrix = sigma.cc,
            abs.quadratic.matrix = delta * e.rr
        )
        
    } else if (interval.algebra == "moore" && restriction == "centre-uncorrelated") {
        result <- .spca.modified.quadratic.by.absolute.values(
            result,
            conv.quadratic.matrix = sigma.cc,
            abs.quadratic.matrix = delta * e.rr,
            orthogonality.matrix = sigma.cc
        )
        
    ### Vector Space Algebra
    } else if (interval.algebra == "vector space" && restriction == "orthogonal") {
        result <- .spca.conventional.case(result, cov.k(sigma.cc, e.rr, model, delta))
        
    } else if (interval.algebra == "vector space" && restriction == "centre-uncorrelated") {
        result <- .spca.projection.case(
            result,
            conventional.matrix = cov.k(sigma.cc, e.rr, model, delta),
            orthogonality.matrix = sigma.cc
        )
        
    ### Error out on unknown algebras/restrictions.
    } else {
        stop(paste(interval.algebra, restriction, delta, model, "not implemented."))
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


# Computes the robust symbolic principal components.
# The columns in the centers and ranges input matrices denote different variables and rows define observations.
rspca <- function(
    interval.algebra = c("vector space", "extended", "moore"),
    restriction = c("orthogonal", "uncorrelated"),
    model = c("diagonal", "full"),
    delta,
    normalise = FALSE,
    C,
    R
) {
    interval.algebra <- match.arg(interval.algebra)
    restriction <- match.arg(restriction)
    if (nrow(C) != nrow(R)) stop("C and R matrices have different numbers of rows.")
    if (ncol(C) != ncol(R)) stop("C and R matrices have different numbers of columns.")
    model <- match.arg(model)
    
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")
    
    result <- list(
        robust = robust.cov.k(C, R, model, delta)
    )
    
    spca <- spca(
        interval.algebra, restriction,
        model, delta, normalise,
        C = C,
        R = R,
        mu.c = result$robust$mu.c,
        sigma.cc = result$robust$sigma.cc,
        mu.r = result$robust$mu.r,
        sigma.rr = result$robust$sigma.rr
    )
    
    # As seen in https://stackoverflow.com/a/37856431/2828287 to merge lists and update existing values.
    return(utils::modifyList(result, spca))
}

.spca.conventional.case <- function(result, conventional.matrix) {
    partial.result <- eigen(conventional.matrix, symmetric = TRUE)
    result$vectors <- partial.result$vectors
    result$values <- partial.result$values
    return(result)
}

.spca.projection.case <- function(result, conventional.matrix, orthogonality.matrix = diag(ncol(result$vectors))) {
    
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

.spca.modified.quadratic.by.absolute.values <- function(
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
        if (lambda.s[s.star] > 0) {
            result$vectors[, i] <- alpha.s[, s.star]
        } else {
            # Find any vector that is orthogonal to all the columns in `orthogonalized`.
            # Iterate over the vectors in the canonical basis and pick the first that is linearly independent.
            col <- ncol(orthogonalized) + 1
            larger.orthogonalized <- matrix(orthogonalized, ncol = col, nrow = nrow(orthogonalized))
            for (j in safe.colon(1, p)) {
                larger.orthogonalized[, col] <- rep(0, p)
                larger.orthogonalized[j, col] <- 1
                new.orthogonalized <- gram.schmidt(larger.orthogonalized)
                candidate <- new.orthogonalized[, col]
                if (sum(candidate*candidate) > 10^(-6)) {
                    result$vectors[, i] <- candidate
                    break
                }
            }
        }
    }
    
    return(result)
}
