rm(list = ls())
load("./scripts/robust_sim_study_no_outliers.RData")

library(dplyr)
library(ggplot2)

Ps <- length(values[[1]])
MINP <- length(values[[1]][[1]]$robust.MRE)
MAXP <- length(values[[1]][[Ps]]$robust.MRE)

# Transform the data into a data frame with the means ----

# Columns are (model, maxp, p, classical/robust, mre)
col.names <- c("model", "maxp", "p", "measure", "robustness", "value")
data <- matrix(NA, nrow = 0, ncol = length(col.names))
colnames(data) <- col.names
for (model in 1:length(values)) {
    model_name <- names(values)[model]
    for (p_ in 1:Ps) {
        maxp <- length(values[[model]][[p_]]$robust.MRE)
        sub <- matrix(NA, nrow = 4*maxp, ncol = length(col.names))
        sub[, 1] <- model_name
        sub[, 2] <- maxp
        sub[, 3] <- rep(1:maxp, 4)
        # Add the MRE information
        sub[1:(2*maxp),  4] <- "MRE"
        sub[1:maxp,      5] <- "robust"
        sub[1:maxp,      6] <- values[[model]][[p_]]$robust.MRE
        sub[maxp+1:maxp, 5] <- "classical"
        sub[maxp+1:maxp, 6] <- values[[model]][[p_]]$classical.MRE
        # Add the ACV information
        o <- 2*maxp
        sub[o + 1:(2*maxp), 4] <- "ACV"
        sub[o + 1:maxp,      5] <- "robust"
        sub[o + 1:maxp,      6] <- values[[model]][[p_]]$robust.ACV
        sub[o + maxp+1:maxp, 5] <- "classical"
        sub[o + maxp+1:maxp, 6] <- values[[model]][[p_]]$classical.ACV
        
        data <- rbind(data, sub)
    }
}
data <- data.frame(data) %>% transform(
    maxp = as.numeric(as.character(maxp)),
    p = as.numeric(as.character(p)),
    value = as.numeric(as.character(value))
)

# Plot mean MRE across models ----

data %>%
    filter(measure == "ACV") %>%
    group_by(model, robustness, maxp) %>%
    summarise(mean = mean(value)) %>%
    ggplot(aes(x = maxp, y = 1-mean, colour = model), alpha = 0.8) +
    geom_line(size=0.75) +
    geom_point(size=1.5) +
    facet_grid(rows = vars(robustness)) +
    labs(
        x = "p",
        y = "global (1-ACV)",
        color = "Config."
        #title = "Mean MRE for each value of p and varying models."
    ) +
    scale_color_discrete(labels = c("EA, k = 3", "EA, k = 5", "EA, k = 8", "MA, k = 3", "MA, k = 5", "MA, k = 8"))

# Old plots ----

##############################################################################################
##############################################################################################
# Compare mean MRE for the various algebras/ks, as P varies.

robust.meanMRE <- matrix(0, nrow = length(values), ncol = Ps)
classical.meanMRE <- matrix(0, nrow = length(values), ncol = Ps)
MREdiffs <- matrix(0, nrow = length(values), ncol = Ps)
rownames(classical.meanMRE) <- rownames(robust.meanMRE) <- names(values)
for (model in 1:length(values)) {
    for (p in 1:Ps) {
        robust.meanMRE[model, p] <- mean(values[[model]][[p]]$robust.MRE)
        classical.meanMRE[model, p] <- mean(values[[model]][[p]]$classical.MRE)
        MREdiffs[model, p] <- robust.meanMRE[model, p] - classical.meanMRE[model, p]
    }
}
r <- max(abs(MREdiffs))*1.2
par(xpd=TRUE)
plot(range(MINP:MAXP), range(c(0, r)), type = "n", xlab = "p", ylab = "mean MRE difference")
colors <- hcl.colors(nrow(MREdiffs), "Blue-Red 3")
for (model in 1:nrow(MREdiffs)) {
    lines(MINP - 1 + 1:Ps, MREdiffs[model, ], col = colors[model], lty = model, type = "b")
}
#legend("top", legend = c("E k3", "E k5", "E k8", "M k3", "M k5", "M k8"), col = colors, lty = 1:MINP, title = "model (Extended or Moore)", horiz = TRUE, inset = c(0, -0.2))
legend("bottom", inset = c(0, 0.3), legend = c("E k3", "E k5", "E k8", "M k3", "M k5", "M k8"), col = colors, lty = 1:MINP, title = "model (Extended or Moore)", horiz = TRUE)

print(robust.meanMRE[1, ])
print(classical.meanMRE[1, ])

##############################################################################################
# Print the robust mean MRE for good measure.

r <- max(abs(robust.meanMRE))*1.2
par(xpd=TRUE)
plot(range(MINP:MAXP), range(c(0, r)), type = "n", xlab = "p", ylab = "mean robust MRE")
colors <- hcl.colors(nrow(robust.meanMRE), "Blue-Red 3")
for (model in 1:nrow(robust.meanMRE)) {
    lines(MINP - 1 + 1:Ps, robust.meanMRE[model, ], col = colors[model], lty = model, type = "b")
}
#legend("top", legend = c("E k3", "E k5", "E k8", "M k3", "M k5", "M k8"), col = colors, lty = 1:MINP, title = "model (Extended or Moore)", horiz = TRUE, inset = c(0, -0.2))
legend("bottom", inset = c(0, 0.3), legend = c("E k3", "E k5", "E k8", "M k3", "M k5", "M k8"), col = colors, lty = 1:MINP, title = "model (Extended or Moore)", horiz = TRUE)

##############################################################################################
##############################################################################################


##############################################################################################
##############################################################################################
# Compare mean ACV for the various algebras/ks, as P varies.
# ACV is not a suitable measure for results with the Moore algebra.

robust.meanACV <- matrix(0, nrow = length(values), ncol = Ps)
classical.meanACV <- matrix(0, nrow = length(values), ncol = Ps)
ACVdiffs <- matrix(0, nrow = length(values), ncol = Ps)
rownames(robust.meanACV) <- names(values)
for (model in 1:length(values)) {
    for (p in 1:Ps) {
        robust.meanACV[model, p] <- 1 - mean(values[[model]][[p]]$robust.ACV)
        classical.meanACV[model, p] <- 1 - mean(values[[model]][[p]]$classical.ACV)
        ACVdiffs[model, p] <- robust.meanACV[model, p] - classical.meanACV[model, p]
    }
}
r <- max(abs(ACVdiffs))*1.2
plot(range(MINP:MAXP), range(c(min(ACVdiffs), r)), type = "n", xlab = "p", ylab = "mean ACV difference")
colors <- hcl.colors(nrow(ACVdiffs), "Blue-Red 3")
for (model in 1:(nrow(ACVdiffs)/2)) {
    lines(MINP - 1 + 1:Ps, ACVdiffs[model, ], col = colors[model], lty = model, type = "b")
}
legend("top", legend = c("E k3", "E k5", "E k8"), col = colors, lty = 1:MINP, title = "model (Extended or Moore)", horiz = TRUE)

print(robust.meanACV[1, ])
print(classical.meanACV[1, ])

##############################################################################################
##############################################################################################