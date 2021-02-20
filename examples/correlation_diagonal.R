# Explore what happens in a diagonal model when two variables are the same,
# or very, very similar.

source("spca.R")

# Generate N data points for P variables, all independent.
# The centres follow normal distributions and the ranges follow exponential distributions.

N <- 100
P <- 5

C <- matrix(rnorm(N*P, 1:P), nrow = N, ncol = P, byrow = TRUE)
R <- matrix(rexp(N*P), nrow = N, ncol = P, byrow = TRUE)

# Coefficients to transform some columns
a <- 3.5
b <- 0.5
# 6th column is exact linear function of 1st column,
# and 7th column is close linear function of 2nd column.
C <- cbind(C, a*C[,1]+b, a*C[,2]+b+rnorm(N, sd = 0.1))
R <- cbind(R, a*R[,1]+b, a*R[,2]+b+rnorm(N, sd = 0.1))

result <- spca("vector space", "centre-uncorrelated", "diagonal", 0.5, FALSE, C, R)
