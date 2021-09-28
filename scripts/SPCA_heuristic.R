# Use SPCA to apply the heuristic in Salvador & Nogueira (2014), Customer-Side Detection of Internet-Scale Traffic Redirection

library(dplyr)
library(ggplot2)
library(tidyverse)

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
            labs(x = "date", fill = "relay") +
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

# Apply SPCA and compute performance measures ----

GRID.CHECK <- TRUE
# Whether to print results that are ready to be pasted in a TeX tabular environment (using the booktabs package)
# or to print something nicer for the R session.
PRINT.TEX <- TRUE
# Use these many PCs for the heuristic.
KEEP.PCS <- 4
# Use `sym.pca` or `robust.sym.pca`
PCA.METHOD <- robust.sym.pca

METHODS <- c("extended", "moore", "vector space")
RESTRICTIONS <- c("orthogonal", "uncorrelated")
KS <- 1:8
# Check if we want to go through everything or not
if (!GRID.CHECK) {
    METHODS <- c(readline(prompt = "interval algebra >> "))
    RESTRICTIONS <- c(readline(prompt = "restriction >> "))
}

matrices <- extract.CR(data)

# # Remove the probe ``Chicago 2'' for the target ``Chicago 1''.
# matrices$C <- matrices$C[, c(1, 3:12)]
# matrices$R <- matrices$R[, c(1, 3:12)]

# # Remove the probe ``Frankfurt 2'' for the target ``Frankfurt 1''.
# matrices$C <- matrices$C[, c(1, 2, 4:12)]
# matrices$R <- matrices$R[, c(1, 2, 4:12)]

split <- 0.6
n <- nrow(matrices$C)
split.to <- round(n*split)
resume.at <- split.to + 1
trainC <- matrices$C[1:split.to, ]
trainR <- matrices$R[1:split.to, ]
testC <- matrices$C[resume.at:n, ]
testR <- matrices$R[resume.at:n, ]

print("CHECK IF USING SPCA OR ROBUST SPCA. CHECK IF ANY VARIABLES WERE REMOVED.")
r <- sapply(METHODS, function(METHOD) {
    sapply(RESTRICTIONS, function(RESTRICTION) {
        values <- sapply(KS, function(K) {
            
            # Apply PCA
            result <- PCA.METHOD(K, interval.algebra = METHOD, restriction = RESTRICTION, trainC, trainR)
            # Keep only the first PCs.
            newC <- result$score.C(testC)[, 1:KEEP.PCS, drop=FALSE] # cf. https://stackoverflow.com/a/7352287/2828287
                                                                    # ^ prevents single column to be demoted to vector.
            # Flip the centres if needed.
            for (i in 1:ncol(newC)) {
                if (sum(newC[, i] < 0) > 0.5*nrow(newC)) {
                    newC[, i] <- -newC[, i]
                    result$vectors[, i] <- -result$vectors[, i]
                }
            }
            # Compute ranges with the correct vectors.
            newR <- result$score.R(testR)[, 1:KEEP.PCS, drop=FALSE] 
            
            # Apply the heuristic
            results <- CR.to.vertical(newC, newR, data$timestamp[resume.at:n]) %>%
                        apply.heuristic(1:KEEP.PCS) %>%
                        group_by(timestamp) %>%
                        summarise(is.attack = sum(heuristic) >= ceiling(0.5*KEEP.PCS))
            subdata <- data[resume.at:n, ]
            
            valid <- !is.na(results$is.attack)
            anomalies <- subdata$anomalyQ == 1
            
            # Compute precision metrics.
            Re <- mean(results$is.attack[valid & anomalies])
            Pr <- mean(anomalies[valid & results$is.attack])
            F1 <- (2*Re*Pr)/(Re + Pr)
            FPR <- mean(results$is.attack[valid & !anomalies])
            values <- c(Re, Pr, F1, FPR)
            names(values) <- c("Re", "Pr", "F1", "FPR")
            
            if (!GRID.CHECK) print(values)
            return(values)
        })
        

        if (PRINT.TEX) {
            if (RESTRICTION == "orthogonal") res <- "orthogonal"
            else res <- "centre-uncorrelated"
            
            if (METHOD == "extended") ia <- "ea"
            else if (METHOD == "moore") ia <- "ma"
            else ia <- "vsa"
            
            writeLines("\\midrule")
            writeLines(paste0("\\multicolumn{9}{c}{\\Gls{", res, "-pcs} \\& \\", ia, "} \\\\"))
            writeLines(paste("   \\recall &", to.tex.numeric.table.row(round(values[1, ], 3), pad.with = 6), "\\\\"))
            writeLines(paste("\\precision &", to.tex.numeric.table.row(round(values[2, ], 3), pad.with = 6), "\\\\"))
            writeLines(paste(" \\fmeasure &", to.tex.numeric.table.row(round(values[3, ], 3), pad.with = 6), "\\\\"))
            writeLines(paste("      \\fpr &", to.tex.numeric.table.row(round(values[4, ], 3), pad.with = 6), "\\\\"))
        } else {
            print(
                paste(
                    str_pad(METHOD, 3+max(map_int(METHODS, nchar)), side = "right"),
                    str_pad(RESTRICTION, 3+max(map_int(RESTRICTIONS, nchar)), side = "right")
                )
            )
            print(values)
        }
    })
})

stop("Stopping here to prevent wasteful generation of plots.")

# Plot the PCs for a specific configuration.
algebra <- "moore"
restriction <- "uncorrelated"
k <- 5

# Apply PCA
result <- PCA.METHOD(k, interval.algebra = algebra, restriction = restriction, trainC, trainR)
# Keep only the first PCs.
newC <- result$score.C(testC)[, 1:KEEP.PCS, drop=FALSE] # cf. https://stackoverflow.com/a/7352287/2828287
                                                        # ^ prevents single column to be demoted to vector.
# Flip the centres if needed.
for (i in 1:ncol(newC)) {
    if (sum(newC[, i] < 0) > 0.5*nrow(newC)) {
        newC[, i] <- -newC[, i]
        result$vectors[, i] <- -result$vectors[, i]
    }
}
# Compute ranges with the correct vectors.
newR <- result$score.R(testR)[, 1:KEEP.PCS, drop=FALSE]

# Visualise the PCs
CR.to.vertical(newC, newR, data$timestamp[resume.at:n]) %>%
    # Filter points that are too far away from the mean.
    group_by(probe) %>% filter((abs(C - mean(C))/sd(C) < 2.5) & (abs(R - mean(R))/sd(R) < 2.5)) %>% ungroup() %>%
    # Remove annotations with timestamps that have nothing to do with the test set.
    basic.plot(annotations %>% filter(xmax > data$timestamp[resume.at])) +
    facet_grid(rows = vars(probe), scales = "free")

# Apply the heuristic.
results <- CR.to.vertical(newC, newR, data$timestamp[resume.at:n]) %>%
    apply.heuristic(1:KEEP.PCS) %>%
    group_by(timestamp) %>%
    summarise(votes = sum(heuristic), is.attack = votes >= ceiling(0.5*KEEP.PCS))
# Determine what type of classification is given to each point.
subdata <- data[resume.at:n, ]
valid <- !is.na(results$is.attack)
anomalies <- subdata$anomalyQ == 1

# Type 0 - regular, classified as regular
# Type 1 - regular, classified as attack
# Type 2 - attack, classified as attack
# Type 3 - attack, classified as regular
results$type <- "false negative" # Type 3 by default
results$type[!valid] <- "true regular"
results$type[valid & !results$is.attack & !anomalies] <- "true regular"
results$type[valid & results$is.attack & !anomalies] <- "false positive"
results$type[valid & results$is.attack & anomalies] <- "true attack"
results$type <- factor(results$type)

shadow.attacks(
    results %>%
        filter(!is.na(is.attack)) %>%
        ggplot(aes(x = timestamp, y = votes, color = type)) +
        # Force the y axis to always go up to the maximum voting number.
        ylim(0, KEEP.PCS) +
        labs(x = "date", fill = "relay") +
        geom_point(), (annotations %>% filter(xmax > min(results$timestamp)))
)

# Compute precision metrics.
Re <- mean(results$is.attack[valid & anomalies])
Pr <- mean(anomalies[valid & results$is.attack])
F1 <- (2*Re*Pr)/(Re + Pr)
FPR <- mean(results$is.attack[valid & !anomalies])
values <- c(Re, Pr, F1, FPR)
names(values) <- c("Re", "Pr", "F1", "FPR")
print(values)
