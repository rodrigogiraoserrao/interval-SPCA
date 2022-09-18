#' Symbolic PCA based on the CIPCA method.
#'
#' Principal component analysis for symbolic objects described by interval variables - CIPCA (Complete Information PCA) method.
#'
#' @inheritParams CPCA
#' @return \code{CIPCA} returns a list containing the following components:
#' \item{centers}{A matrix of centers.}
#' \item{pc}{A 3-dimensional array of the interval principal components.}
#' \item{mat}{A covariance (or correlation) matrix.}
#' \item{loadings}{A matrix whose columns contain the eigenvectors.}
#' \item{percvar}{A Table of eigenvalues and percentage of variance explained.}
#' \item{cjv}{The correlations between the original interval variables and the interval principal components.}

#' @seealso \code{\link{SimpleFormPCA}}, \code{\link{CIPCAcovMcd_B}}, \code{\link{CIPCAcovMcd_C}}.
#' @export
CIPCA <- function (data,npc=2,Cor=FALSE) 
{
    n<-dim(data)[1]
    p<-dim(data)[2]
    centr <- as.matrix(data[, , 1] + data[, , 2])/2
    rownames(centr)<-c(1:n)
    colnames(centr)<-c(paste("Var", 1:p, sep = "_"))
    Mcentr<-apply(centr,2,mean)
    Mmean<-matrix(rep(Mcentr,n),n,p,byrow=T)
    Centt<-data
    Centt[,,1]<-data[, , 1]-Mmean
    Centt[,,2]<-data[, , 2]-Mmean
    
    ya<-matrix(0,n,npc)
    yb<-matrix(0,n,npc)
    
    var<-matrix(0,3,npc)
    colnames(var) <- c(paste("PC", 1:npc, sep = "_"))
    rownames(var) <- c("lambda i","Proportion Var","Cumul. Var")
    Cjv<-matrix(0,p,npc)
    colnames(Cjv) <- c(paste("PC", 1:npc, sep = "_"))
    rownames(Cjv) <- c(paste("Var", 1:p, sep = "_"))
    
    #inner product
    INP<-function(data,j,k)
    {
        res<-0
        if(j!=k){
            for(i in 1:dim(data)[1])
            {
                res<-res + 0.25*sum(data[i,j,])*sum(data[i,k,])
            }
        }
        else{
            for(i in 1:dim(data)[1])
            {
                res<-res + (data[i,j,1]^2+data[i,j,1]*data[i,j,2]+data[i,j,2]^2)/3
            }
        }
        res
    }
    
    
    
    if(Cor==TRUE){
        
        Standt<-data
        #       for(j in 1:p)
        #       {
        #         VaXj<-sqrt((1/n)*INP(Centt,j,j))
        #         for(i in 1:n)
        #         {
        #         Standt[i,j,1]<-(data[i,j,1]-Mcentr[j])/VaXj
        #         Standt[i,j,2]<-(data[i,j,2]-Mcentr[j])/VaXj
        #       }
        #     }
        sdX<-rep(0,p)
        for(j in 1:p)
        {
            sdX[j]<-sqrt((1/n)*INP(Centt,j,j))
        }
        for(i in 1:n)
        {
            for(j in 1:p)
            {
                Standt[i,j,1]<-(data[i,j,1]-Mcentr[j])/sdX[j]
                Standt[i,j,2]<-(data[i,j,2]-Mcentr[j])/sdX[j]
            }
        }
        
        CorM<-matrix(NA,p,p)
        for(i in 1:p)
        {
            for(j in 1:p)
            {
                CorM[i,j]<-(1/n)*INP(Standt,i,j)
            }
        }
        
        mat<-CorM
        aux<-eigen(mat)
        vec<-aux$vectors
        lambd<-aux$values
        
        xa<-matrix(0,n,p)
        xb<-matrix(0,n,p)
        for(pci in 1:npc) {
            for(i in 1:n) {
                for(j in 1:p) {
                    if(vec[j,pci]>0){
                        xa[i,j]<-Standt[i,j,1]
                        xb[i,j]<-Standt[i,j,2]
                    }
                    else{
                        xa[i,j]<-Standt[i,j,2]
                        xb[i,j]<-Standt[i,j,1]
                    }
                }
            }
            ya[,pci]<-xa%*%vec[,pci]
            yb[,pci]<-xb%*%vec[,pci]
            
            var[1,pci]<-lambd[pci]
            var[2,pci]<-lambd[pci]/sum(lambd)
            
            if(pci==1){
                var[3,pci]<-var[2,pci]
            }
            else{
                var[3,pci]<-var[3,(pci-1)]+var[2,pci]
            }
            for(pci in 1:npc) {
                for(j in 1:p) {
                    Cjv[j,pci]<-vec[j,pci]*sqrt(lambd[pci])
                }
            }
        }
    }
    #if Cor==FALSE 
    else{
        
        CovM<-matrix(NA,p,p)
        for(i in 1:p)
        {
            for(j in 1:p)
            {
                CovM[i,j]<-(1/n)*INP(Centt,i,j)
            }
        }
        
        mat<-CovM
        aux<-eigen(mat)
        vec<-aux$vectors
        lambd<-aux$values
        
        xa<-matrix(0,n,p)
        xb<-matrix(0,n,p)
        for(pci in 1:npc) {
            for(i in 1:n) {
                for(j in 1:p) {
                    if(vec[j,pci]>0){
                        xa[i,j]<-Centt[i,j,1]
                        xb[i,j]<-Centt[i,j,2]
                    }
                    else{
                        xa[i,j]<-Centt[i,j,2]
                        xb[i,j]<-Centt[i,j,1]
                    }
                }
            }
            ya[,pci]<-xa%*%vec[,pci]
            yb[,pci]<-xb%*%vec[,pci]
            
            var[1,pci]<-lambd[pci]
            var[2,pci]<-lambd[pci]/sum(lambd)
            
            if(pci==1){
                var[3,pci]<-var[2,pci]
            }
            else{
                var[3,pci]<-var[3,(pci-1)]+var[2,pci]
            }
            for(pci in 1:npc) {
                for(j in 1:p) {
                    Cjv[j,pci]<-vec[j,pci]*sqrt(lambd[pci]/CovM[j,j])
                }
            }
        }
    }
    
    pc <-array(NA,c(n,npc,2))
    pc[,,1]<-ya
    pc[,,2]<-yb
    dimnames(pc)<-list(c(1:n),c(paste("PC", 1:npc, sep = "_")),c("Min","Max"))
    colnames(mat)<-c(paste("Var", 1:p, sep = "_"))
    rownames(mat)<-c(paste("Var", 1:p, sep = "_"))
    colnames(vec)<-c(paste("E", 1:p, sep = "_"))
    list(centers=centr,pc = pc,mat=mat,loadings=vec,percvar=var,cjv=Cjv)
}

##############################

source("RTTutils.R")

load("data/dataRTT_all.RData")

PCS <- 1  # How many principal components to use?

base.data <- data.T4.all.noNA
data <- rtt.to.horizontal.symbolic(base.data) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin="1970-01-01"))
matrices <- extract.CR(data)
C <- matrices$C
Rad <- matrices$R / 2 # Radius is half the range.

cipca.data <- array(dim = c(nrow(C), ncol(C), 2))
cipca.data[,, 1] <- C - Rad
cipca.data[,, 2] <- C + Rad

# Split in training and testing
split <- 0.6
n <- nrow(C)
train.idx <- 1:round(n * split)
test.idx <- (round(n * split) + 1):n
train <- cipca.data[train.idx,, ]
test <- cipca.data[test.idx,, ]

fit <- CIPCA(train, npc=PCS)

# Project the interval endpoints according to the PCs.
new_data <- array(dim=c(length(test.idx), PCS, 2))
new_data[,, 1] <- test[,, 1] %*% fit$loadings[1:12, 1:PCS, drop=FALSE]
new_data[,, 2] <- test[,, 2] %*% fit$loadings[1:12, 1:PCS, drop=FALSE]

# Recompute centres and ranges.
newC <- as.matrix((new_data[,, 1] + new_data[,, 2]) / 2, nrow = length(test.idx), ncol = PCS)
newR <- as.matrix(new_data[,, 2] - new_data[,, 1], nrow = length(test.idx), ncol = PCS)
# Flip the centres and ranges if needed.
for (i in 1:ncol(newC)) {
    if (sum(newC[, i] < 0) > 0.5*nrow(newC)) {
        print("Flipping")
        newC[, i] <- -newC[, i]
        newR[, i] <- -newR[, i]
    }
}
print(sum(newR < 0))

# Apply the heuristic
results <- CR.to.vertical(newC, newR, data$timestamp[test.idx]) %>%
    apply.heuristic(1:PCS) %>%
    group_by(timestamp) %>%
    summarise(is.attack = sum(heuristic) >= ceiling(0.5*PCS))
subdata <- data[test.idx, ]

valid <- !is.na(results$is.attack)
anomalies <- subdata$anomalyQ == 1

# Compute precision metrics.
Re <- mean(results$is.attack[valid & anomalies])
Pr <- mean(anomalies[valid & results$is.attack])
F1 <- (2*Re*Pr)/(Re + Pr)
FPR <- mean(results$is.attack[valid & !anomalies])
values <- c(Re, Pr, F1, FPR)
names(values) <- c("Re", "Pr", "F1", "FPR")
print(values)
