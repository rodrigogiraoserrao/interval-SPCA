d <- function(I, p) {
    # Returns the I-th sign vector of size p.
    (-1)^as.numeric(intToBits(I))[1:p]
}

d.factory <- function(p) {
    # Returns a function that can generate the sign vectors of size p.
    d.curried <- function(I) { d(I, p) }
    d.curried
}

sgn <- function(vector) {
    # Implements the sign(v) function for a vector v.
    sign(0.5 + sign(vector))
}

gram.schmidt <- function(vector.matrix) {
    # Implements the Gram-Schmidt orthogonalization process over the columns of vector.matrix.
    result <- matrix(normalize(vector.matrix[,1]), ncol=1)
    for (i in 2:ncol(vector.matrix)) {
        print("The result is")
        print(result)
        # Use the modified Gram-Schmidt process to reduce floating-point inaccuracies
        v <- vector.matrix[,i]
        print("v is")
        print(v)
        for (k in 1:ncol(result)) {
            v <- v - (v%*%result[,k])[1,1]*result[,k]
        }
        result <- cbind(result, normalize(v))
    }
    result
}

normalize <- function(v) {
    # Normalize a vector to norm 1. If the vector has a norm very close to 0 return the 0 vector.
    norm <- sqrt(v%*%v)[1,1]
    if (round(norm, digits=5) > 0) {
        v/norm
    } else {
        0*v
    }
}