source("spca.R")

# Generate N data points for P variables, all independent.
# The centres follow normal distributions and the ranges follow exponential distributions.

N <- 10
P <- 5

C <- matrix(rnorm(N*P, 1:P), nrow = N, ncol = P, byrow = TRUE)
R <- matrix(rexp(N*P), nrow = N, ncol = P, byrow = TRUE)

# Compute the symbolic covariance matrix.
sigma.cc <- cov(C)
mu.r <- to.column(colMeans(R))
e.rr <- cov(R) + mu.r%*%t(mu.r)
scov.matrix <- cov.k(sigma.cc, e.rr, "diagonal", 0.5)
print(scov.matrix)

# Normalise the data by hand.
A <- diag.matrix(diag(scov.matrix)^-.5)
norm.sigma.cc <- A%*%sigma.cc%*%A
norm.e.rr <- A%*%e.rr%*%A
print(cov.k(norm.sigma.cc, norm.e.rr, "diagonal", 0.5))
