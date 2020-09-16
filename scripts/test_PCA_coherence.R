### This script performs a basic check to see if the covariance of the PCs matches
## the theoretical formula I found.

# author: Rodrigo Girão Serrão
# email: rodrigogiraoserrao@gmail.com
# date: 16/09/2020

source("./sym.R")
source("./sym_pca.R")
source("./utils.R")
source("./RTTutils.R")

library(tidyverse)


# Load and prepare the data -----------------------------------------------

load("data/dataRTT_all.RData")

data <- data.T1.all.noNA %>% rtt.to.symbolic  # dataset to use
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

# Symbolic global variables
K <- 3
METHOD <- "vector space"
RESTRICTION <- "uncorrelated"
NORMALIZE <- TRUE

sigma.cc <- cov(C)
e.rr <- (t(R)%*%R)/nrow(R)

# Apply PCA
result <- sym.pca(C, R, K, interval.algebra = METHOD, restriction = RESTRICTION, normalize = NORMALIZE)
newC <- result$score.C(C)
newR <- result$score.R(R)
newCov <- cov.k(C, R, K)

if (is.diagonal.k(K)) {
    cov <- t(result$vectors) %*% sigma.cc %*% result$vectors
} else if (METHOD == "vector space") {
    cov <- t(result$vectors) %*% (sigma.cc + delta.k(K)*e.rr) %*% result$vectors
} else if (METHOD == "extended") {
    a.r <- abs(t(R) %*% result$vectors)
    cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*((t(a.r)%*%a.r)/nrow(a.r))
} else if (METHOD == "moore") {
    cov <- (t(result$vectors) %*% sigma.cc %*% result$vectors) + delta.k(K)*(abs(t(result$vectors)) %*% e.rr %*% abs(result$vectors))
} else {
    stop("You don't know what you are doing!")
}

print(isTRUE(all.equal(newCov, cov)))