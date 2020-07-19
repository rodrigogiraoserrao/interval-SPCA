##
#
# This module implements distances for interval valued variables.
#
##

# Sum the mahalanobis distances of the lower and upper bounds of the intervals.
Dist.ULM <- function(C, R) {
    # cf. Dynamic Cluster Methods for Interval Data based on Mahalanobis Distances, section 2
    p <- ncol(C)
    QL <- cov(C-R)
    ML <- det(QL)^(1/p)*solve(QL)
    QU <- cov(C+R)
    MU <- det(QU)^(1/p)*solve(QU)
    dist.ULM <- function(XC, XR, YC, YR) {
        XL <- XC - XR
        XU <- XC + XR
        YL <- YC - YR
        YU <- YC + YR
        return(sqrt(
            t(XL-YL)%*%ML%*%(XL-YL) + t(XU-YU)%*%MU%*%(XU-YU)
        ))
    }
    return(dist.ULM)
}