# Simulation study for ROBUST SPCA
# Shifts means of the distribution of the centres.

source("./utils.R")
source("./sym.R")
source("./sym_pca.R")

# Set the seed for reproducibility
set.seed(73)

# Generate N data points
N <- 100
M <- 500
P <- 4

RESTRICTION <- "orthogonal"

mu.c <- 5*(1:P)
sigma.cc <- diag((1:P)^2/25)
sd.c <- sqrt(diag(sigma.cc))
mu.r <- rep(10, P) + ((1:P)-1)/2
sigma.rr <- diag((1:P)^2/100)
sd.r <- sqrt(diag(sigma.rr))

contaminated.mu.c <- mu.c
contaminated.mu.r <- mu.r

contaminated.sd.c <- sd.c
contaminated.sd.r <- sd.r

## UNCOMMENT the contaminations wanted.
#contaminated.mu.c[P] <- contaminated.mu.c[P] + mu.c[1]
#contaminated.mu.r[P] <- contaminated.mu.r[P] + mu.r[1]
#contaminated.sd.c <- 4*contaminated.sd.c
#contaminated.sd.r <- 4*contaminated.sd.r

#contaminated.mu.c <- 10*contaminated.mu.c
#contaminated.mu.r <- 10*contaminated.mu.r

# Iterate over contamination levels
result <- sapply(c(5, 10, 20), function(epsilon) {
    print(paste0("epsilon = ", epsilon))
    result <- sapply(c("extended", "moore"), function(algebra) {
        print(algebra)
        result <- sapply(c(3, 5, 8), function(k) {
            print(k)
            baseline <- sym.pca(
                k, interval.algebra = algebra, restriction = RESTRICTION,
                mu.c = mu.c, sigma.cc = sigma.cc, mu.r = mu.r, sigma.rr = sigma.rr
            )
            baseline.lengths <- diag(t(baseline$vectors)%*%baseline$vectors)
            
            result <- replicate(M, {
                C <- matrix(rnorm(P*N, mean = c(mu.c), sd = sd.c), nrow = N, ncol = P, byrow = TRUE)
                R <- matrix(rnorm(P*N, mean = c(mu.r), sd = sd.r), nrow = N, ncol = P, byrow = TRUE)
                
                # Contaminate the data
                for (i in 1:round(epsilon*N/100)) {
                    C[N+1-i, ] <- rnorm(P, mean = c(contaminated.mu.c), sd = contaminated.sd.c)
                    R[N+1-i, ] <- rnorm(P, mean = c(contaminated.mu.r), sd = contaminated.sd.r)
                }
                
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

filename <- readline(prompt = "Save to what .RData file? >> ")
if (nchar(filename) > 0) {
    filepath <- paste0("./scripts/", filename, ".RData")
    save(list = "values", file = filepath)
}
