# Function to generate the sign vectors of length p.
# If I is 0 we return a matrix with all ds in columns.
d <- function(p, I = 0) {
    if (I == 0) {
        (-1)^matrix(as.numeric(intToBits(1:2^p - 1)), nrow = 32)[1:p,]
    } else {
        (-1)^as.numeric(intToBits(I-1))[1:p]
    }
}

# Curries the parameter p into the d function.
d.factory <- function(p) {
    d.curried <- function(I) { d(p, I) }
    d.curried
}

# Implements the vectorized sign function so that |v|_i = v_i * sign(v)_i.
sgn <- function(vector) {
    sign(0.5 + sign(vector))
}

# Function that builds a diagonal matrix out of an input vector v or an input matrix M.
## If the input is a vector, creates a diagonal matrix with v as main diagonal.
## If the input is a matrix, creates a diagonal matrix with the same diagonal as M.
diag.matrix <- function(vM) {
    if (is.matrix(vM)) vM <- diag(vM)
    diag(vM)
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