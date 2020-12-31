### This script performs basic checks to see if SPCA preservers total variability and generalised variance.

source("./sym.R")
source("./sym_pca.R")
source("./utils.R")
source("./RTTutils.R")

library(tidyverse)


# Load and prepare the data -----------------------------------------------

load("data/dataRTT_all.RData")

data <- data.T3.all.noNA %>% rtt.to.horizontal.symbolic  # dataset to use
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

base.dets <- rep(0, length(KS))
base.trs <- rep(0, length(KS))
first.values <- rep(0, length(KS))
for (k in KS) {
    svar <- estimate.cov.k(C, R, k)
    det <- determinant(svar)
    base.dets[k] <- det$sign * as.numeric(det$modulus)
    base.trs[k] <- sum(diag(svar))
    first.values[k] <- max(diag(svar))
}

print("Results are (original measure) * (1+value) = (PC measure)")

r <- sapply(METHODS, function(METHOD) {
    sapply(RESTRICTIONS, function(RESTRICTION) {
        checks <- sapply(KS, function(K) {
            values <- c(0, 0)
            names(values) <- c("tr", "log det")
            
            # Apply PCA
            result <- sym.pca(K, interval.algebra = METHOD, restriction = RESTRICTION, C, R)
            newC <- result$score.C(C)
            newR <- result$score.R(R)
            newCov <- estimate.cov.k(newC, newR, K)
            
            #### SANITY CHECKS START HERE
            
            tr <- sum(diag(newCov))
            values[1] <- round(100*(tr - base.trs[K])/base.trs[K], 2)
            
            det <- determinant(newCov)
            det <- det$sign * as.numeric(det$modulus)
            values[2] <- round(100*(det - base.dets[K])/base.dets[K], 2)
            
            if (!GRID.CHECK) print(c(tr, det))
            return(values)
        })
        
        print(
            paste(
                str_pad(METHOD, 3+max(map_int(METHODS, nchar)), side = "right"),
                str_pad(RESTRICTION, 3+max(map_int(RESTRICTIONS, nchar)), side = "right")
            )
        )
        print(checks)
        print(to.tex.numeric.table(checks))
    })
})

print("Results are (original measure) * (1+value) = (PC measure)")
