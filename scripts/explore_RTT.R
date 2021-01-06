library(dplyr)
library(ggplot2)

source("RTTutils.R")
source("sym_pca.R")
source("sym.R")

# Auxiliary functions ----

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

# Prepare data ----

load("data/dataRTT_all.RData")

base.data <- data.T4.all.noNA
sorted.probes <- sort(probes)

vdata <- rtt.to.vertical.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01")) %>%
    apply.heuristic(sorted.probes)

print(nrow(vdata))
# vdata <- vdata %>% filter(C + R < 800)
# print(nrow(vdata))

data <- rtt.to.horizontal.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01"))


# Find the attacks to annotate later plots ----

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


# Exploratory plot (cf. Subtil) ----

print(nrow(vdata) - nrow(vdata %>% filter(C + R < 800)))
vdata %>% filter(C + R < 800) %>%
    basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free") +
    labs(fill = "Relay")


# SMALL exploratory plot (cf. Subtil) for extended abstract ----

start_at <- which.max(6 == lubridate::day(data$timestamp)) # Plot from 6th of June forward.
vdata %>% filter(C + R < 800) %>% filter(probe %in% c("Amsterdam","VdM")) %>% filter(timestamp >= data$timestamp[start_at]) %>%
    basic.plot(annotations %>% filter(xmin > data$timestamp[start_at])) +
    facet_grid(rows = vars(probe), scales = "free") +
    labs(fill = "Relay")


# Exploratory plot (cf. Subtil) with heuristic ----

shadow.attacks(
    vdata %>% filter(C + R < 800, !is.na(heuristic)) %>%
        ggplot(aes(x = timestamp, y = C, color = heuristic)) +
        geom_point(size = 0.5) +
        facet_grid(rows = vars(probe), scales = "free"), annotations
)


# Aggregated vote with heuristic ----

shadow.attacks(
    vdata %>% group_by(timestamp) %>%
        summarise(heuristic = sum(heuristic)) %>%
        mutate(is.attack = heuristic >= 6) %>%
        ggplot(aes(x = timestamp, y = heuristic, color = is.attack)) +
        geom_point(size = 0.5), annotations
)


# Plot Chicago for T1 or Frankfurt for T2 because it looks interesting ----

# For T1, the target is Chicago1 and so the probe Chicago2 is (probably) in the
## physical machine, meaning the RTTs are 0 except when there's relays.
# Same thing for T2, but with Frankfurt.

probe.name <- "Frankfurt2"
vdata %>% filter(probe == probe.name) %>%
    basic.plot(annotations) +
    ggtitle(paste("Probe", probe.name)) +
    geom_line(mapping = aes(x = timestamp, y = C))


# Plot original variable with higher variance ----

matrices <- extract.CR(data)
cov <- estimate.cov.k(matrices$C, matrices$R, 3)
i <- which.max(diag(cov))

# Original plot
data.frame(timestamp = data$timestamp, C = matrices$C[, i], R = matrices$R[, i]) %>%
    basic.plot(annotations) +
    ggtitle(paste0(sorted.probes[i], " -- probe with greatest variability"))

# Remove average C before plotting
avg <- mean(matrices$C[, i])
shadow.attacks(
    data.frame(timestamp = data$timestamp, C = matrices$C[, i]-avg, R=matrices$R[, i]) %>%
        ggplot(aes(x = timestamp, ymin = C-R, ymax = C+R)) +
        geom_ribbon() +
        ggtitle(paste0(sorted.probes[i], " -- greatest var, recentred")), annotations
)


# Prepare for SPCA ----

K <- 3
IA <- "extended"
RESTRICTION <- "uncorrelated"

p <- length(sorted.probes)
n <- nrow(data)

# Apply SPCA and rSPCA and plot first variable ----

r <- sym.pca(k = K, interval.algebra = IA, restriction = RESTRICTION,
             C = matrices$C, R = matrices$R)
newC <- r$score.C(matrices$C)
newR <- r$score.R(matrices$R)

# Regular plot
shadow.attacks(
    data.frame(timestamp = data$timestamp, C = newC[, 1], R = newR[, 1]) %>%
    ggplot(aes(x = timestamp, ymin = C-R, ymax = C+R)) +
    geom_ribbon() +
    ggtitle("First PC variable"), annotations
)
# Remove the average
avg <- mean(newC[, 1])
shadow.attacks(
    data.frame(timestamp = data$timestamp, C = newC[, 1]-avg, R = newR[, 1]) %>%
        ggplot(aes(x = timestamp, ymin = C-R, ymax = C+R)) +
        geom_ribbon() +
        ggtitle("First PC variable, without mean"), annotations
)
# Remove the average in a robust way
r <- robust.cov.k(newC, newR, K)
shadow.attacks(
    data.frame(timestamp = data$timestamp, C = newC[, 1]-r$mu.c[1], R = newR[, 1]) %>%
        ggplot(aes(x = timestamp, ymin = C-R, ymax = C+R)) +
        geom_ribbon() +
        ggtitle("First PC variable, without robust mean"), annotations
)


# Vertical data from the PCs ----

pc.vdata <- CR.to.vertical(newC, newR, data$timestamp) %>% apply.heuristic(1:length(sorted.probes))
pc1.vdata <- pc.vdata %>% filter(probe == 1)


# Plot all PC variables ----

pc.vdata %>% basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free") +
    ggtitle("All PC variables.")


# Plot first 4 PC variables ----

pc.vdata %>% filter(probe <= 4) %>%
    basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free") +
    ggtitle("First 4 PC variables.")


# Plot mid 4 PC variables ----

pc.vdata %>% filter(4 < probe & probe <= 8) %>%
    basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free") +
    ggtitle("Mid 4 PC variables.")


# Plot last 4 PC variables ----

pc.vdata %>% filter(probe > 8) %>%
    basic.plot(annotations) +
    facet_grid(rows = vars(probe), scales = "free") +
    ggtitle("Last 4 PC variables.")


# Plot heuristic for first PC ----

shadow.attacks(
    pc.vdata %>%
        filter(probe == 1) %>%
        filter(!is.na(heuristic)) %>%
        ggplot(aes(x = timestamp, y = C, color = heuristic)) +
        geom_point(size = 0.5), annotations
)
# and compare performance with benchmark data
valid <- !is.na(pc1.vdata$heuristic)
anomalies <- data[, ]$anomalyQ == 1

mean(pc1.vdata[valid& anomalies, ]$heuristic == anomalies[valid& anomalies])
mean(pc1.vdata[valid&!anomalies, ]$heuristic == anomalies[valid&!anomalies])
