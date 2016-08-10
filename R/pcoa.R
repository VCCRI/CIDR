#' @title Principal Coordinates Analysis
#'
#' @rdname cidrPcoa
#' @name cidrPcoa
#'
#' @description Performs principal coordinates analysis on a dissimilarity matrix.
#'
#' @details
#' This function is used internally in the method \code{scPCA};
#' see the help page of \code{scPCA} for more details.
#'
#' This function is a modified version of the \emph{pcoa} function from the CRAN
#' package \emph{ape (version 3.5)}.  The function was modified on 6th August 2016.
#'
#'
#' @param D a dissimilarity matrix
#'
#' @return This function is used internally in the method \code{scPCA};
#' see the help page of \code{scPCA} for more details.
#'
#' @references
#' Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
#'
#' @export cidrPcoa
#'

cidrPcoa <- function (D) {
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  bstick.def <- function(n, tot.var = 1, ...) {
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    return(res)
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  names <- rownames(D)

  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1)
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0

  if (min.eig > -epsilon) {
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k]/trace
    vectors <- sweep(D.eig$vectors[, 1:k], 2, sqrt(eig[1:k]), FUN = "*")
    res <- data.frame(eig[1:k], rel.eig)
    colnames(res) <- c("Eigenvalues", "Relative_eig")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
    out <- (list(values = res, vectors = vectors, trace = trace))
  }
  else {
    k <- n
    eig <- D.eig$values
    rel.eig <- eig/trace
    k2 <- length(which(eig > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2], 2, sqrt(eig[1:k2]), FUN = "*")
    res <- data.frame(eig[1:k], rel.eig)
    colnames(res) <- c("Eigenvalues", "Relative_eig")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
    out <- (list(values = res, vectors = vectors, trace = trace))
    }
  return(out)
}
