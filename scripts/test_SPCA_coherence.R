### This script performs a basic check to see if the covariance of the PCs matches
## the theoretical formula I found.

source("./sym.R")
source("./sym_pca.R")
source("./utils.R")
source("./RTTutils.R")

library(tidyverse)


# Load and prepare the data -----------------------------------------------

load("data/dataRTT_all.RData")

data <- data.T1.all.noNA %>% rtt.to.horizontal.symbolic  # dataset to use
# Delete heavy variables from memory.
remove(
    list=c(
        paste0("data.T", 1:4, ".all"),
        paste0("data.T", 1:4, ".all.noNA")
    )
)
C <- as.matrix(data %>% select(!!paste0("C", 1:12)))
colnames(C) <- probes
R <- as.matrix(data %>% select(!!paste0("R", 1:12)))
colnames(R) <- probes

GRID.CHECK <- TRUE # set this to enable/disable thorough checking

sigma.cc <- cov(C)
e.rr <- cov(R) + to.column(colMeans(R))%*%t(colMeans(R))

METHODS <- c("vector space", "extended", "moore")
RESTRICTIONS <- c("orthogonal", "uncorrelated")
KS <- 1:8
# Check if we want to go through everything or not
if (!GRID.CHECK) {
    METHODS <- c(readline(prompt = "method >> "))
    RESTRICTIONS <- c(readline(prompt = "restriction >> "))
}

sapply(METHODS, function(METHOD) {
    sapply(RESTRICTIONS, function(RESTRICTION) {
        checks <- sapply(KS, function(K) {
            checks <- c(FALSE, FALSE)
            names(checks) <- c("covs match", "valid restriction")
            
            # Apply PCA
            result <- sym.pca(K, interval.algebra = METHOD, restriction = RESTRICTION, C, R)
            newC <- result$score.C(C)
            newR <- result$score.R(R)
            newCov <- estimate.cov.k(newC, newR, K)
            
            #### SANITY CHECKS START HERE
            
            # Check if the covariances computed both ways match.
            if (is.diagonal.k(K)) {
                if (METHOD == "moore") cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*diag.matrix(abs(t(result$vectors)) %*% e.rr %*% abs(result$vectors))
                else cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*diag.matrix(t(result$vectors) %*% e.rr %*% result$vectors)
            } else if (METHOD == "vector space") {
                cov <- t(result$vectors) %*% (sigma.cc + delta.k(K)*e.rr) %*% result$vectors
            } else if (METHOD == "extended") {
                a.r <- abs(R %*% result$vectors)
                a.rr <- cov(a.r) + to.column(colMeans(a.r))%*%t(colMeans(a.r))
                cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*a.rr
            } else if (METHOD == "moore") {
                cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*(abs(t(result$vectors)) %*% e.rr %*% abs(result$vectors))
            } else {
                stop("You don't know what you are doing!")
            }
            checks[1] <- isTRUE(all.equal(cov, newCov, tolerance = 10^(-6)))
            
            # Check if the PCs satisfy the restriction they should
            if (RESTRICTION == "orthogonal")  {
                valid.restriction <- isTRUE(all.equal(t(result$vectors) %*% result$vectors, diag(nrow(newCov)), tolerance = 10^(-6)))
            } else if (RESTRICTION == "uncorrelated") {
                cov.newC <- cov(newC)
                valid.restriction <- isTRUE(all.equal(cov.newC, diag.matrix(cov.newC), tolerance = 10^(-6)))
            }
            checks[2] <- valid.restriction
            if (!GRID.CHECK) print(checks)
            return(checks)
        })
        
        # check which values of K passed all tests
        test.results <- apply(checks, 2, all)
        print(
            paste(
                ifelse(all(test.results), "V", "X"),
                str_pad(METHOD, 3+max(map_int(METHODS, nchar)), side = "right"),
                str_pad(RESTRICTION, 3+max(map_int(RESTRICTIONS, nchar)), side = "right"),
                paste(ifelse(test.results, ".", "X"), collapse = "")
            )
        )
    })
})