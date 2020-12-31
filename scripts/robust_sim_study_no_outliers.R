# Simulation study for ROBUST SPCA
# Comparison with classical SPCA in the absence of outliers.

source("./utils.R")
source("./sym.R")
source("./sym_pca.R")

# Set the seed for reproducibility
set.seed(73)

# Generate N data points
N <- 100
M <- 500
MINP <- 4
MAXP <- 10

RESTRICTION <- "orthogonal"

MU.C <- 5*(1:MAXP)
SIGMA.CC <- diag((1:MAXP)^2/25)
MU.R <- rep(10, MAXP) + ((1:MAXP)-1)/2
SIGMA.RR <- diag((1:MAXP)^2/100)

result <- sapply(MINP:MAXP, function(P) {
    mu.c <- to.column(MU.C[1:P])
    sigma.cc <- SIGMA.CC[1:P, 1:P]
    mu.r <- to.column(MU.R[1:P])
    sigma.rr <- SIGMA.RR[1:P, 1:P]
    print(P)
    
    sd.c <- sqrt(diag(sigma.cc))
    sd.r <- sqrt(diag(sigma.rr))
    
    result <- sapply(c("extended", "moore"), function(algebra) {
        result <- sapply(c(3, 5, 8), function(k) {
            baseline <- sym.pca(
                k, interval.algebra = algebra, restriction = RESTRICTION,
                mu.c = mu.c, sigma.cc = sigma.cc, mu.r = mu.r, sigma.rr = sigma.rr
            )
            baseline.lengths <- diag(t(baseline$vectors)%*%baseline$vectors)
            
            result <- replicate(M, {
                C <- matrix(rnorm(P*N, mean = c(mu.c), sd = sd.c), nrow = N, ncol = P, byrow = TRUE)
                R <- matrix(rnorm(P*N, mean = c(mu.r), sd = sd.r), nrow = N, ncol = P, byrow = TRUE)
                # Ensure we did not generate negative ranges
                assertthat::assert_that(sum(R <= 0) == 0)                
                
                robust.r <- robust.sym.pca(k, interval.algebra = algebra, restriction = RESTRICTION, C = C, R = R)
                robust.MRE <- abs(baseline$values - robust.r$values)/baseline$values
                robust.cumvar <- var.explained(robust.r$values)
                robust.ACV <- abs(diag(t(abs(baseline$vectors))%*%abs(robust.r$vectors))/(baseline.lengths*diag(t(robust.r$vectors)%*%robust.r$vectors)))
                
                classical.r <- sym.pca(k, interval.algebra = algebra, restriction = RESTRICTION, C = C, R = R)
                classical.MRE <- abs(baseline$values - classical.r$values)/baseline$values
                classical.cumvar <- var.explained(classical.r$values)
                classical.ACV <- abs(diag(t(abs(baseline$vectors))%*%abs(classical.r$vectors))/(baseline.lengths*diag(t(classical.r$vectors)%*%classical.r$vectors)))
                
                return(list(
                    robust.MRE = robust.MRE,
                    robust.cumvar = robust.cumvar,
                    robust.ACV = robust.ACV,
                    classical.MRE = classical.MRE,
                    classical.cumvar = classical.cumvar,
                    classical.ACV = classical.ACV
                ))
            })
            
            # Helper function to summarise across the M replications.
            summarise <- function(result) {
                classical.MRE <- classical.cumvar <- classical.ACV <- matrix(0, nrow = M, ncol = P)
                robust.MRE <- robust.cumvar <- robust.ACV <- matrix(0, nrow = M, ncol = P)
                sapply(1:M, function (i) {
                    robust.MRE[i, ] <<- result[, i]$robust.MRE
                    robust.cumvar[i, ] <<- result[, i]$robust.cumvar
                    robust.ACV[i, ] <<- result[, i]$robust.ACV
                    classical.MRE[i, ] <<- result[, i]$classical.MRE
                    classical.cumvar[i, ] <<- result[, i]$classical.cumvar
                    classical.ACV[i, ] <<- result[, i]$classical.ACV
                })
                return(list(
                    baseline.cumvar = var.explained(baseline$values),
                    robust.MRE = colMeans(robust.MRE),
                    robust.cumvar = colMeans(robust.cumvar),
                    robust.ACV = colMeans(robust.ACV),
                    classical.MRE = colMeans(classical.MRE),
                    classical.cumvar = colMeans(classical.cumvar),
                    classical.ACV = colMeans(classical.ACV)
                ))
            }
            result <- summarise(result)
            result$baseline <- baseline
            return(result)
        })
        return(list(
            k3 = result[, 1],
            k5 = result[, 2],
            k8 = result[, 3]
        ))
    })
    return(list(
        extended_k3 = result[, 1]$k3,
        extended_k5 = result[, 1]$k5,
        extended_k8 = result[, 1]$k8,
        moore_k3 = result[, 2]$k3,
        moore_k5 = result[, 2]$k5,
        moore_k8 = result[, 2]$k8
    ))
})

values <- list(
    extended_k3 = result[1, ],
    extended_k5 = result[2, ],
    extended_k8 = result[3, ],
    moore_k3 = result[4, ],
    moore_k5 = result[5, ],
    moore_k8 = result[6, ]
)

save(list = "values", file = "./scripts/robust_sim_study_no_outliers.RData")
