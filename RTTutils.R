library(dplyr)

# Take an RTT data set and build the symbolic variables,
# each probe corresponds to two columns Cx and Rx.
rtt.to.horizontal.symbolic <- function(df) {
    # Reshape the data into symbolic data.
    data <- df
    sorted <- sort(probes)
    for (j in 1:12) {
        i <- which(probes == sorted[j])
        c <- paste0("C", j)
        r <- paste0("R", j)
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

# Take an RTT data set and build the symbolic variables,
# each row has a probe and a timestamp, so that there's
# only two columns: C and R.
rtt.to.vertical.symbolic <- function(df) {
    # Reshape the data into symbolic data.
    data <- data.frame()
    sorted <- sort(probes)
    for (j in 1:12) {
        i <- which(probes == sorted[j])
        m <- paste0("mRTT", i)
        M <- paste0("MRTT", i)
        data <- bind_rows(data,
            df %>% mutate(
                C := (.data[[M]] + .data[[m]])/2,
                R := (.data[[M]] - .data[[m]]),
                probe = probes[i]
            ) %>%
                select(timestamp, anomalyQ, relay, relayname, C, R, probe)
        )
    }
    data
}

# Takes the symbolic RTT data and makes it vertical.
make.vertical <- function(df) {
    data <- data.frame()
    sorted <- sort(probes)
    for (j in 1:12) {
        i <- which(probes == sorted[j])
        c <- paste0("C", j)
        r <- paste0("R", j)
        data <- bind_rows(data,
                          df %>% mutate(
                              C := .data[[c]],
                              R := .data[[r]],
                              probe = probes[i]
                          ) %>%
                              select(timestamp, anomalyQ, relay, relayname, C, R, probe)
        )
    }
    data
}

# Takes the horizontal data and extracts a matrix for the centres and a matrix for the ranges
extract.CR <- function(hdata) {
    n <- nrow(hdata)
    nam <- names(hdata)
    p <- min(sum(startsWith(nam, "C")), sum(startsWith(nam, "R")))
    C <- matrix(0, nrow = n, ncol = p)
    R <- matrix(0, nrow = n, ncol = p)
    for (i in 1:p) {
        C[, i] <- data[[paste0("C", i)]]
        R[, i] <- data[[paste0("R", i)]]
    }
    return(list(C = C, R = R))
}

# Takes the C and R matrices and builds a vertical data set
CR.to.vertical <- function(C, R, timestamps) {
    df <- data.frame()
    for (i in 1:ncol(C)) {
        df <- bind_rows(df,
                        data.frame(timestamp = timestamps, C = C[, i], R = R[, i], probe = i)
        )
    }
    return(df)
}

# Simple Moving Average
sma <- function(x, n){stats::filter(x, rep(1/n, n), sides=1)}

# Takes a vertical RTT dataset and classifies attacks according to the heuristic seen in [1]
# [1] Salvador, P. and Nogueira, A. (2014). Customer-side detection of Internet-scale traffic redirection.
## In 2014 16th International Telecommunications Network Strategy and Planning Symposium (Networks),
## pages 1–5.
apply.heuristic <- function(vdata, probes) {
    vdata$heuristic <- NA
    N <- 480 # trailing length for moving average for the heuristic
    K <- 10  # how many consecutive times do we need to exceed the threshold?
    for (i in 1:length(probes)) {
        idx <- which(vdata$probe == probes[i])
        Cs <- vdata[idx, ]$C
        values <- rep(0, length(Cs))
        values[1:N] <- Cs[1:N]
        s <- 1
        threshold <- 1.2*mean(values[s:N])
        above <- rep(NA, length(Cs))
        j <- N+1
        while (j <= length(Cs)) {
            if (Cs[j] > threshold) {
                above[j] <- TRUE
            } else {
                above[j] <- FALSE
                s <- s + 1
                values[s+N] <- Cs[j]
                threshold <- 1.2*mean(values[s:(s+N)])
            }
            j <- j + 1
        }
        vdata[idx, ]$heuristic <- K == stats::filter(above, rep(1, K), sides = 1)
    }
    return(vdata)
}
