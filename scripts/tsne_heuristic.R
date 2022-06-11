# Based off of https://github.com/Xuxl2020/t-SNE-classifier/blob/4e08403c4f46de0bdddbd9df428fd759240c1a3b/tsne%20with%20AD.R
# Implementation of out-of-sample t-SNE embedding:
# A t-SNE Based Classification Approach to Compositional Microbiome Data
# Xu Xueli, Xie Zhongming, Yang Zhenyu, Li Dongfang, Xu Ximing
# 10.3389/fgene.2020.620143

library(Rtsne)
library(FNN)
library(dplyr)
library(ggplot2)

source("RTTutils.R")


load("data/dataRTT_all.RData")

base.data <- data.T3.all.noNA

data <- rtt.to.horizontal.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01"))
matrices <- extract.CR(data)
# Use 24 variables (centres + ranges) or just centres
# matrix.data <- cbind(matrices$C, matrices$R)
matrix.data <- matrices$C

# Distances
dists <- as.matrix(dist(matrix.data))

# Split in training and testing
split <- 0.6
n <- nrow(matrix.data)
train.idx <- 1:round(n*split)
test.idx <- (round(n * split) + 1):n
train <- matrix.data[train.idx, ]
test <- matrix.data[test.idx, ]

# # Make a nice plot with the projected data.
# projected <- data.frame(
#     timestamp = data$timestamp[1:split.to],
#     value = Rtsne(train, dims = 1, num_threads = 0)$Y,
#     anomalyQ = data$anomalyQ[1:split.to]
# )
# 
# projected %>% ggplot(aes(x = value.1, y = value.2, color = anomalyQ)) + geom_point()

k <- 10; per <- 480; dim <- 1; iter <- 2000

# Find K nearest neighbours based on Euclidean distance
test.dist <- dists[test.idx, train.idx]
test.neighbours.index <- t(apply(test.dist, 1, function(x) {order(x)[1:k]}))
test.neighbours.dist <- t(apply(test.dist, 1, function(x) {sort(x)[1:k]}))

# Model fit
tsne.fit = as.matrix(
    Rtsne(train, dims = dim, max_iter = iter, num_threads = 0)$Y
)

# Project test data.
projections <- matrix(nrow = nrow(test), ncol = dim)
for (i in 1:nrow(test)) {
    # Compute weights
    w0 = exp(-test.neighbours.dist[i,])
    w0[w0 < 1e-12] = 1e-12
    w = w0 / sum(w0)
    # Get approximate projection
    projections[i, ] <- w %*% tsne.fit[test.neighbours.index[i,],]
}

# Find low-dim nearest neighbours.
neighbs <- get.knnx(tsne.fit, projections, k = k)$nn.index
# Are the low-dim neighbours anomalies or not?
neighbs.anomalies <- matrix(
    as.integer(as.character(data$anomalyQ))[neighbs],
    nrow = nrow(neighbs),
    ncol = ncol(neighbs),
)  # Do a sort of 2D indexing
is.anomaly <- apply(neighbs.anomalies, 1, sum) >= (k / 2)  # Sum each row and check for majority voting

anomalies <- data$anomalyQ[test.idx] == 1
Re <- mean(is.anomaly[anomalies])
Pr <- mean(anomalies[is.anomaly])
F1 <- (2*Re*Pr)/(Re + Pr)
FPR <- mean(is.anomaly[!anomalies])
print(c(Re, Pr, F1, FPR))

# data.to.plot <- data.frame(
#     timestamp = data$timestamp[test.idx],
#     value = projections,
#     anomalyQ = is.anomaly
# )
# data.to.plot %>% ggplot(aes(x = timestamp, y = value, color = anomalyQ)) + geom_point()
