source("./sym.R")
source("./sym_pca.R")

library(tidyverse)

# Take an RTT data set and build the symbolic variables.
rtt.to.symbolic <- function(df) {
    # Reshape the data into symbolic data.
    data <- df
    for (i in 1:12) {
        c <- paste0("C", i)
        r <- paste0("R", i)
        m <- paste0("mRTT", i)
        M <- paste0("MRTT", i)
        data <- data %>%
            mutate(
                !!c := (.data[[M]] + .data[[m]])/2,
                !!r := (.data[[M]] - .data[[m]])
            ) %>%
            select(-!!c(m, M, paste0(c("avgRTT", "MedRTT", "StdDev"), i)))
    }
    return(data)
}

# Helper function to set the axis of the plots
stretch.axis <- function(axis) {
    axis.range <- axis[2] - axis[1]
    delta <- 0.1*axis.range
    c(axis[1]-delta, axis[2]+delta)
}

# Helper function to plot two variables, one VS the other
plot.VS <- function(C, R, vars, border = NULL) {
    # Plots the vars[1] variable VS vars[2].
    xaxis <- stretch.axis(c(
        min(C[, vars[1]]-abs(R[, vars[1]])),
        max(C[, vars[1]]+abs(R[, vars[1]]))
    ))
    yaxis <- stretch.axis(c(
        min(C[, vars[2]]-abs(R[, vars[2]])),
        max(C[, vars[2]]+abs(R[, vars[2]]))
    ))
    plot(
        xaxis,
        yaxis,
        xlab = colnames(C)[vars[1]],
        ylab = colnames(C)[vars[2]],
        main = paste("Intervals for", colnames(C)[vars[1]], "and", colnames(C)[vars[2]]),
    )
    rect(
        xleft = C[, vars[1]]-R[, vars[1]],
        ybottom = C[, vars[2]]-R[, vars[2]],
        xright = C[, vars[1]]+R[, vars[1]],
        ytop = C[, vars[2]]+R[, vars[2]],
        border = border
    )
    legend("topright", legend=as.character(clevels), fill=colours, title="relay")
}

# Helper function to plot interval variables
plot.interval <- function(C, R, var, col = NULL) {
    # Plots the intervals of the variable var.
    xaxis <- stretch.axis(c(
        min(C[, var]-abs(R[, var])),
        max(C[, var]+abs(R[, var]))
    ))
    plot_height <- max(R[, var]) - min(R[, var])
    plot(
        xaxis,
        c(
            min(R[, var]),
            max(R[, var])
        ),
        xlab = colnames(C)[var],
        ylab = "Range",
        main = paste("Intervals for", colnames(C)[var], "with y = range"),
    )
    rect(
        xleft = C[, var]-R[, var],
        ybottom = R[, var]-plot_height/100,
        xright = C[, var]+R[, var],
        ytop = R[, var]+plot_height/100,
        density = 50,
        col = col
    )
}

load("data/dataRTT_all.RData")

data <- data.T1.all.noNA %>% rtt.to.symbolic  # dataset to use
# Delete heavy variables from memory.
remove(
    list=c(
        paste0("data.T", 1:4, ".all"),
        paste0("data.T", 1:4, ".all.noNA")
    )
)
# Plot blue rectangles from regular observations and red rectangles for anomalies.
colours <- c("chartreuse", paste0("indianred", 1:4))
clevels <- levels(data$relayname)

makePlots <- FALSE
conv <- c(1, 2)  # plot conv[1]th VS conv[2]th conventional variables
pc <- c(1, 2)    # plot pc[1]th VS pc[2]th principal components
# Symbolic global variables
K <- 5
algebra <- "moore"
restriction <- "orthogonal"

C <- as.matrix(data %>% select(!!paste0("C", 1:12)))
colnames(C) <- probes
R <- as.matrix(data %>% select(!!paste0("R", 1:12)))
colnames(R) <- probes

# Prepare plotting area
if (makePlots) par(mfrow=c(2, 3))

# Find the 2 columns with greatest variance and plot their intervals.
Cov <- cov.k(C, R, K)
vars <- diag(Cov)
trace <- sum(vars)

if (makePlots) {
    # Sort the conventional variables by variance.
    conv_ <- sort.int(vars, decreasing = TRUE, index.return = TRUE)$ix[conv]
    plot.VS(C, R, conv_, border = colours[match(data$relayname, clevels)])
    # Now plot the intervals of the variable conv[1].
    plot.interval(C, R, conv_[1], colours[match(data$relayname, clevels)])
    # Now plot the intervals of the variable conv[2].
    plot.interval(C, R, conv_[2], colours[match(data$relayname, clevels)])
}
    
# Call the symbolic PCAs here.

result <- sym.pca(C, R, K, interval.algebra = algebra, restriction = restriction, normalize = TRUE)
newC <- C %*% result$vectors
colnames(newC) <- paste("pc", 1:ncol(newC))
newR <- R %*% result$vectors
newCov <- cov.k(newC, newR, K)
newVars <- diag(newCov)
newTrace <- sum(newVars)
print(paste("abs(newTrace - trace) = ", abs(newTrace - trace)))
print(round(100*cumsum(newVars)/newTrace, 2))
absR <- abs(newR)

if (makePlots) {
    # Plot the rectangles for the PCs in pc
    plot.VS(newC, newR, pc, border = colours[match(data$relayname, clevels)])
    # Now only plot the pc[1] with y depending on the range of the interval.
    plot.interval(newC, newR, pc[1], col = colours[match(data$relayname, clevels)])
    # Now only plot the pc[2].
    plot.interval(newC, newR, pc[2], col = colours[match(data$relayname, clevels)])
}

# # Project the PCs into a subspace.
# keep <- 10
# 
# projC <- 0*newC
# projC[, 1:keep] <- newC[, 1:keep]
# projR <- 0*newR
# projR[, 1:keep] <- newR[, 1:keep]
# 
# diff <- projC - newC
# distsC <- sqrt(diag(diff %*% t(diff)))
# diff <- projR - newR
# distsR <- sqrt(diag(diff %*% t(diff)))
# 
# data %>%
#     select(anomalyQ, relayname) %>%
#     add_column(dC = distsC, dR = distsR) %>%
#     ggplot(aes(x = dC, y = dR, color = relayname)) +
#     scale_color_manual(values = colours) +
#     geom_point(alpha = 0.6) +
#     labs(color = "relay name")