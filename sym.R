##
#
# Define functions relevant to symbolic analysis of interval data.
#
##

# (TODO) rename functions that are prefixed with .k to something more meaningful,
# e.g. var.k and cov.k are terrible function names.

# We need some of the symbolic utility functions.
source("./utils.R")


#' Symbolic variance of a vector of intervals.
#'
#'Estimates the symbolic variance of interval-valued random variables based on centre and range samples.
#'
#' @param C The data matrix with the population centres of the intervals. Each row represents an observation and columns represent variables.
#' @param R The data matrix with the population ranges of the intervals. Each row represents an observation and columns represent variables.
#' @param delta The non-negative factor that is multiplied with the matrix of the second moments of the ranges when computing the symbolic covariance matrix.
#' 
#' @return \code{var.k} Returns a vector with the estimated symbolic variances.
#' 
#' @references Rodrigo Girão Serrão, M. Rosário Oliveira, and Lina Oliveira (20XX), Theoretical Derivation of Interval Principal Component Analysis
#'@seealso \code{\link{spca}}, \code{\link{cov.k}}, \code{\link{estimate.cov.k}}.
#' @export
var.k <- function(C, R, delta) {
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")

    return(diag(var(C)) + delta*diag(t(R)%*%R)/nrow(R))
}

#' Symbolic covariance matrix.
#'
#'Computes the symbolic covariance of interval-valued random variables based on the covariance of the centres and the second moments of the ranges.
#'
#' @param sigma.cc The covariance matrix of the centres.
#' @param e.rr The second moments matrix of the ranges.
#' @param model Whether the symbolic covariance matrix uses the full matrix of the second moments of the ranges of the intervals, or just the diagonal part.
#' @param delta The non-negative factor that is multiplied with the matrix of the second moments of the ranges when computing the symbolic covariance matrix.
#' 
#' @return \code{spca} Returns a matrix with the estimated symbolic covariance.
#' 
#' @references Rodrigo Girão Serrão, M. Rosário Oliveira, and Lina Oliveira (20XX), Theoretical Derivation of Interval Principal Component Analysis
#'@seealso \code{\link{spca}}, \code{\link{var.k}}, \code{\link{estimate.cov.k}}.
#' @export
cov.k <- function(sigma.cc, e.rr, model = c("diagonal", "full"), delta) {
    model <- match.arg(model)
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")

    if (model == "full") return(sigma.cc + delta*e.rr)
    else return(sigma.cc + delta*diag.matrix(e.rr))
}


#' Symbolic covariance matrix estimation.
#'
#'Estimates the symbolic covariance of interval-valued random variables based on centre and range samples.
#'
#' @param C The data matrix with the population centres of the intervals. Each row represents an observation and columns represent variables.
#' @param R The data matrix with the population ranges of the intervals. Each row represents an observation and columns represent variables.
#' @param model Whether the symbolic covariance matrix uses the full matrix of the second moments of the ranges of the intervals, or just the diagonal part.
#' @param delta The non-negative factor that is multiplied with the matrix of the second moments of the ranges when computing the symbolic covariance matrix.
#' 
#' @return \code{spca} Returns a matrix with the estimated symbolic covariance.
#' 
#' @references Rodrigo Girão Serrão, M. Rosário Oliveira, and Lina Oliveira (20XX), Theoretical Derivation of Interval Principal Component Analysis
#'@seealso \code{\link{spca}}, \code{\link{cov.k}}, \code{\link{var.k}}.
#' @export
estimate.cov.k <- function(C, R, model = c("diagonal", "full"), delta) {
    model <- match.arg(model)
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")
    
    sigma.cc <- cov(C)
    mu.r <- to.column(colMeans(R))
    sigma.rr <- cov(R)
    return(cov.k(sigma.cc, sigma.rr + mu.r%*%t(mu.r), model, delta))
}

# Estimates the symbolic covariance matrix in a robust way.
robust.cov.k <- function(C, R, model = c("diagonal", "full"), delta, ...) {
    model <- match.arg(model)
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")
    
    r <- robustbase::covMcd(cbind(C, R), ...)
    p <- ncol(C)
    r$mu.c <- to.column(r$center[1:p])
    r$sigma.cc <- r$cov[1:p, 1:p]
    r$mu.r <- to.column(r$center[(p+1):(2*p)])
    r$sigma.rr <- r$cov[(p+1):(2*p), (p+1):(2*p)]
    r$e.rr <- r$sigma.rr + r$mu.r%*%t(r$mu.r)
    
    if (model == "full") r$cov.k <- r$sigma.cc + delta*r$e.rr
    else r$cov.k <- r$sigma.cc + delta*diag.matrix(r$e.rr)
    
    return(r)
}

# Estimates the correlation between interval valued variables.
# The variables of the first (C, R) pair are along the rows and the variables of the second (C, R) pair are along the columns.
estimate.corr.k <- function(C1, R1, C2, R2, model = c("diagonal", "full"), delta) {
    model <- match.arg(model)
    if (!is.atomic(delta) || length(delta) > 1) stop("Delta should be a single number.")
    if (!is.numeric(delta)) stop("Delta should be a non-negative real number.")
    if (!delta >= 0) stop("Delta should be non-negative.")
    
    covs.c <- cov(C1, C2)
    result <- covs.c
    
    if (model == "full") result <- result + delta*t(R1)%*%R2/nrow(R1)
    
    sds1 <- sqrt(diag(estimate.cov.k(C1, R1, model, delta)))
    sds2 <- sqrt(diag(estimate.cov.k(C2, R2, model, delta)))
    
    # Divide by the standard deviations.
    result <- result / sds1
    result <- sweep(result, 2, sds2, `/`)
    
    return(result)
}

warn.about.linear.dependences <- function(C, R) {
    # Issue warnings for the pairs (i, j) of interval-valued random variables
    # that are very similar, in the sense that one is a linear function of the other.
    # To determine these pairs of variables, we perform linear regression on the centres
    # of the variables and for the pairs of centres that exhibit linear dependence we
    # use the same regression coefficients on the corresponding ranges.
    # If those coefficients also represent a good linear regression for the ranges, then
    # we take that to mean the interval valued variables are essentially the same
    # and we issue a warning.
    
    R.means <- colMeans(R)
    for (i in 1:(ncol(C)-1)) {
        for (j in (i+1):ncol(C)) {
            obj <- lm(C[,i] ~ C[,j])
            r.squared <- suppressWarnings(summary(obj)$r.squared)
            if (r.squared < 0.95) next
            
            # Try to apply the same coefficients to the ranges.
            coefs <- obj$coefficients
            Ri <- coefs[2]*R[,j] + coefs[1]
            res <- R[,i] - Ri
            r.squared <- 1 - sum(res^2)/sum((R[,i]-R.means[i])^2)
            if (r.squared < 0.95) next
            
            warning(paste("Interval valued variables", i, "and", j, "are very similar."))
        }
    }
}