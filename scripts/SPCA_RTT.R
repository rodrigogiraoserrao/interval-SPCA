# Thoroughly apply (r)SPCA to the RTT datasets.

library(tidyverse)
library(dplyr)
library(ggplot2)

source("RTTutils.R")
source("sym_pca.R")
source("sym.R")
source("utils.R")

# Prepare the data ----

load("data/dataRTT_all.RData")

base.data <- data.T1.all.noNA
sorted.probes <- sort(probes)

vdata <- rtt.to.vertical.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01")) %>%
    apply.heuristic(sorted.probes)

data <- rtt.to.horizontal.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01"))

# Apply SPCA ----

K <- 3
IA <- "extended"
RESTRICTION <- "uncorrelated"

p <- length(sorted.probes)
n <- nrow(data)
matrices <- extract.CR(data)

# Uncomment the two lines below if using the robust.sym.pca method.
# matrices$C <- matrices$C[, c(1, 3:12)]
# matrices$R <- matrices$R[, c(1, 3:12)]
r <- robust.sym.pca(k = K, interval.algebra = IA, restriction = RESTRICTION,
             C = matrices$C, R = matrices$R)

# Symbolic correlation between Johannesburg 1 and Johannesburg 2 is 0.9566833
# cov.matrix <- estimate.cov.k(matrices$C, matrices$R, K)
# vars <- diag(cov.matrix)
# (cov.matrix[6,7]/sqrt(vars[6]))/sqrt(vars[7])

print("Coefficient vectors.")
writeLines(to.tex.numeric.table(round(r$vectors, 3)))

writeLines(to.tex.numeric.table.row(round(100*var.explained(r$values), 2), 6))

newC <- r$score.C(matrices$C)
newR <- r$score.R(matrices$R)

print("Correlations between vars and pcs.")
writeLines(to.tex.numeric.table(round(estimate.corr.k(matrices$C, matrices$R, newC, newR, K), 3)))



# Prepare to plot stuff ----

shadow.attacks <- function(plot, annotations) {
    return(
        plot + geom_rect(
            data = annotations,
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lbl),
            alpha = 0.2,
            inherit.aes = FALSE
        )
    )
}

basic.plot <- function(data, annotations) {
    shadow.attacks(
        data %>% ggplot(aes(x = timestamp, ymin = C-R, ymax = C+R)) +
            scale_x_datetime() +
            labs(x = "date") +
            geom_ribbon(), annotations
    )
}

## Find the attacks to annotate plots ----

s <- 1
xmins <- c()
xmaxs <- c()
lbls <- c()
while (s < nrow(data)) {
    if (data[s,]$relayname != "reg") {
        lbl <- data[s,]$relayname
        end <- s
        while (data[end,]$relayname == lbl) {
            end <- end + 1
        }
        xmins <- c(xmins, data[s,]$timestamp)
        xmaxs <- c(xmaxs, data[end,]$timestamp)
        lbls <- c(lbls, lbl)
        s <- end + 1
    } else s <- s + 1
}
annotations <- data.frame(
    ymin = -Inf,
    ymax = Inf,
    xmin = xmins,
    xmax = xmaxs,
    lbl = lbls
) %>% mutate(xmin = as.POSIXct(xmin, origin = "1970-01-01"),
             xmax = as.POSIXct(xmax, origin = "1970-01-01"),
             lbl = factor(relays[lbl - 1]))


pc.vdata <- CR.to.vertical(newC, newR, data$timestamp)

pc.vdata %>% group_by(probe) %>% filter((abs(C - mean(C))/sd(C) < 2.5) & (abs(R - mean(R))/sd(R) < 2.5)) %>% ungroup() %>%
    basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free")

# pc.vdata %>% filter(probe <= 4) %>%
#     basic.plot(annotations) +
#     facet_grid(rows = vars(probe), scales = "free")

# filtered <- pc.vdata %>% filter(probe <= 4) %>% group_by(probe) %>% filter((abs(C - mean(C))/sd(C) < 2.5) & (abs(R - mean(R))/sd(R) < 2.5)) %>% ungroup()
# print(nrow(filtered))
# print(nrow(pc.vdata %>% filter(probe <= 4)))
# 
# filtered %>%
#     basic.plot(annotations) +
#     facet_grid(rows = vars(probe), scales = "free")
