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