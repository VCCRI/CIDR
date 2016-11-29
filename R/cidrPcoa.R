#' @title Principal Coordinates Analysis
#'
#' @rdname cidrPcoa
#' @name cidrPcoa
#'
#' @description performs principal coordinates analysis on a dissimilarity matrix.
#'
#' @details
#' This function is used internally in the method \code{scPCA};
#' see the help page for \code{scPCA} for more details.
#'
#'
#' @param D a dissimilarity matrix.
#'
#' @return This function is used internally in the method \code{scPCA};
#' see the help page for \code{scPCA} for more details.
#'
#'
#' @export cidrPcoa
#' @import RcppEigen
#'

cidrPcoa <- function (D) {
    centre <- function(D, n) {
        One <- matrix(1, n, n)
        mat <- diag(n) - One/n
        ## use RcppEigen to calculate product: (mat * D * mat)
        mat.cen <- eigen_centre(mat, D)$centre
    }
    D <- as.matrix(D)
    n <- nrow(D)
    epsilon <- sqrt(.Machine$double.eps)
    names <- rownames(D)
    
    delta1 <- centre((-0.5 * D^2), n)
    ## getEigenSpace is c++ code that uses Rcpp & RcppEigen
    D.eig <- getEigenSpace(delta1)
    ## getEigenSpace sorts from smallest to largest
    min.eig <- D.eig$values[1]
    zero.eig <- which(abs(D.eig$values) < epsilon)
    D.eig$values[zero.eig] <- 0
    ## reverse the eigenspace - so that it is sorted descending
    eig_values <- D.eig$values[length(D.eig$values):1]
    eig_vectors <- D.eig$vectors[,ncol(D.eig$vectors):1]
    ## only interested in positive eigenvalues
    k <- length(which(eig_values > epsilon))
    vectors <- sweep(eig_vectors[, 1:k], 2, sqrt(eig_values[1:k]), FUN = "*")
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
    ## note returning eig_values for values - which will include 
    ## all eigenvalues - not just the positive ones
    out <- (list(values = eig_values, vectors = vectors))
    return(out)
}
