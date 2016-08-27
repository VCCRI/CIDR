#' @title Clustering through Imputation and Dimensionality Reduction
#'
#' @description Fast and accurate clustering through imputation and dimensionality
#' reduction for single cell RNA-Seq data.
#'
#' @author Peijie Lin <P.Lin@victorchang.edu.au>, Michael Troup
#'
#' @docType package
#' @name cidr-package
#' @aliases cidr
#' @useDynLib cidr
#'
#' @examples
#' par(ask=FALSE)
#' ## Generate simulated single cell RNA-Seq tags.
#' N=3 ## 3 cell types
#' k=50 ## 50 cells per cell type
#' sData <- scSimulator(N=N, k=k)
#' ## tags - the tag matrix
#' tags <- as.matrix(sData$tags)
#' cols <- c(rep("RED",k), rep("BLUE",k), rep("GREEN",k))
#' ## Standard principal component analysis.
#' ltpm <- log2(t(t(tags)/colSums(tags))*1000000+1)
#' pca <- prcomp(t(ltpm))
#' plot(pca$x[,c(1,2)],col=cols,pch=1,xlab="PC1",ylab="PC2",main="prcomp")
#' ## Use cidr to analyse the simulated dataset.
#' ## The input for cidr should be a tag matrix.
#' sData <- scDataConstructor(tags)
#' sData <- determineDropoutCandidates(sData)
#' sData <- wThreshold(sData)
#' sData <- scDissim(sData)
#' sData <- scPCA(sData)
#' nCluster(sData)
#' sData <- scCluster(sData)
#' ## Two dimensional visualization: different colors denote different cell types,
#' ## while different plotting symbols denote the clusters output by cidr.
#' plot(sData@PC[,c(1,2)], col=cols,
#'      pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")
#' ## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
#' adjustedRandIndex(sData@clusters,cols)
#' ## 0.79
#' ## Alter the number of PCs used in the clustering.
#' sData <- scCluster(sData, nPC=2)
#' plot(sData@PC[,c(1,2)], col=cols,
#'      pch=sData@clusters,main="CIDR",xlab="PC1", ylab="PC2")
#' adjustedRandIndex(sData@clusters,cols)
#' ## 0.92
NULL
#' @title Adjusted Rand Index
#' @author Chris Fraley, Adrian Raftery, Luca Scrucca.
#' @rdname adjustedRandIndex
#' @name adjustedRandIndex
#' @description Calculates the Adjusted Rand Index which meansures the accuracy of clustering when the ground truth is known.
#' @importFrom mclust adjustedRandIndex
#' @export adjustedRandIndex
#' @details Imported from the package \emph{mclust}; see \code{adjustedRandIndex::mclust} help page for more details.
#' @references
#' Chris Fraley, Adrian E. Raftery, T. Brendan Murphy, and Luca Scrucca (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation Technical Report No. 597, Department of Statistics, University of Washington
#'
#' Chris Fraley and Adrian E. Raftery (2002) Model-based Clustering, Discriminant Analysis and Density Estimation Journal of the American Statistical Association 97:611-631

## class scData - single cell data object with properties releveant to
## clustering through imputation and dimensionality reduction
setClass("scData", representation(tags="matrix",
                                  sampleSize="numeric",
                                  librarySizes="vector",
                                  nData="matrix",
                                  priorTPM="numeric",
                                  dThreshold="vector",
                                  wThreshold="numeric",
                                  dropoutCandidates="matrix",
                                  PC="matrix",
                                  variation="vector",
                                  eigenvalues="vector",
                                  dissim="matrix",
                                  nCluster="numeric",
                                  clusters="vector",
                                  nPC ="numeric",
                                  cMethod="character",
                                  correction="character"),
         prototype(cMethod="ward.D2", nCluster=0, correction="none", priorTPM=1))

#' @title scData Constructor
#'
#' @rdname scDataConstructor
#'
#' @description
#' \code{scDataConstructor} creates a new scData class object from a tag table.
#'
#' @details
#' Creates an object whose type represents a single cell dataset class.
#' Properties of the class include scalar, vector and matrix
#' data types necessary for the CIDR analysis - such as gene counts, library
#' sizes, thresholds, dropouts, and clustering data.  Raw counts are converted
#' to log2 per million, normalised by dividing each gene count for a
#' particular cell by the total count for all genes in that cell.
#'
#' @param tags a matrix of tags where the rows crrespond to features (genes, transcripts, etc) and the columns correspond to cells
#' @export
#' @return an scData class object
#' @examples
#' ## Generate simulated single cell RNA-Seq tags.
#' N=3 ## 3 cell types
#' k=50 ## 50 cells per cell type
#' sData <- scSimulator(N=N, k=k)
#' ## The input for cidr should be a tag matrix. 
#' tags <- as.matrix(sData$tags)
#' ## create a new scData object
#' sData <- scDataConstructor(tags)
#' ## print the first 5 library sizes
#' sData@librarySizes[1:5]
#' ## print a portion of the data matrix of the class - contains raw tags
#' sData@tags[1:5, 30:34]
#' ## print part of the data matrix of the class - log tag per million
#' sData@nData[1:5, 30:34]
scDataConstructor <- function(tags){
    tags <- tags[rowSums(tags)>0,]
    object <- new("scData", tags=tags)
    object@sampleSize <- ncol(tags)
    object@librarySizes <- colSums(tags)
    object@nData <- log2(t(t(tags)/object@librarySizes)*1000000+object@priorTPM)
    return(object)
}

setGeneric("determineDropoutCandidates", function(object, min1=3, min2=8, N=2000, alpha=0.1, fast=TRUE, zerosOnly=FALSE){
    standardGeneric("determineDropoutCandidates")
})

#' @title Determine Dropout Candidates
#'
#' @rdname determineDropoutCandidates
#' @name determineDropoutCandidates
#'
#' @description
#' determines which entries in a single cell RNA-Seq dataset are potentially dropouts.
#'
#' @details
#' populates a truth matrix with the same dimension as the nData.
#' The value is \code{True} for an individual entry if it
#' is a dropout candidate, otherwise the value is \code{false}.
#'
#' @param object the scData class object
#' @param min1,min2 technical parameters used in estimating the minimum point between the first two modes of the density curve of logTPM for each cell
#' @param alpha a cutoff quantile in the range [0,1]. Thresholds outside this will be adjusted to the quantile boundary.
#' @param N number of cells to consider when determining the threshold value for dropout candidates; used in conjunction with the \code{fast} parameter
#' @param fast boolean; if \code{TRUE} (default), implements a fast version for datasets with a sample size greater than N
#' @param zerosOnly boolean; if \code{TRUE}, only zeros are considered as dropout candidates; by default \code{FALSE}.
#' @export
#' @return an updated scData class object with the following attributes updated
#'
#' \item{dThreshold}{a vector of library dependent dropout candidate thresholds}
#' \item{dropoutCandidates}{a matrix with the same dimension as nData. The value is \code{True} for an individual entry if it
#' is a dropout candidate, otherwise the value is \code{false}.}
#' @examples
#' example(cidr)
setMethod("determineDropoutCandidates", "scData", function(object, min1, min2, N, alpha, fast, zerosOnly){
    if(zerosOnly){
        object@dThreshold <- log2(rep(1, object@sampleSize)/object@librarySizes*1000000+object@priorTPM)
        object@dropoutCandidates <- (object@tags==0)
    } else {
        topLibraries <- 1:object@sampleSize
        if(fast & (object@sampleSize>N)){
            topLibraries <- order(object@librarySizes,decreasing = TRUE)[1:N]
        } else {
            N <- object@sampleSize
        }
        dTs <- rep(0,N)
        LT1 <- log2(rep(min1, N)/object@librarySizes[topLibraries]*1000000+object@priorTPM)
        LT2 <- log2(rep(min2, N)/object@librarySizes[topLibraries]*1000000+object@priorTPM)
        object@dropoutCandidates <- array(NA, dim=dim(object@nData))
        for(i in 1:N){
            dfn_m <- density(object@nData[, topLibraries[i]], kernel="epanechnikov", n=1024, from=LT2[i])
            dfn_max <- dfn_m$x[which.max(dfn_m$y)]
            dfn <- density(object@nData[, topLibraries[i]], kernel="epanechnikov", n=1024, from=LT1[i], to=dfn_max)
            dTs[i] <- dfn$x[which.min(dfn$y)]
        }
        if(fast & (object@sampleSize>N)){
            object@dThreshold <- rep(median(dTs),object@sampleSize)
        } else{
            limits <- quantile(dTs,c(alpha, 1-alpha))
            dTs[dTs<limits[1]] <- limits[1]
            dTs[dTs>limits[2]] <- limits[2]
            object@dThreshold <- dTs
        }
        object@dropoutCandidates <- t(t(object@nData) < object@dThreshold)
    }    
    return(object)
})

setGeneric("wThreshold", function(object, cutoff=0.5){
    standardGeneric("wThreshold")
})

#' @title Imputation Weighting Threshold
#'
#' @description
#' determines the imputation weighting threshold.
#'
#' @rdname wThreshold
#' @name wThreshold
#'
#' @param object the scData class object
#' @param cutoff parameter in the range (0,1), used in the calculation of imputation weighting threshold
#'
#' @importFrom minpack.lm nlsLM
#' @export
#'
#' @return an updated scData class object with the following attribute updated
#'
#' \item{wThreshold}{imputation weighting threshold}
#'
#' @examples
#' example(cidr)
setMethod("wThreshold", "scData", function(object, cutoff){
    delete <- which(rowSums(object@dropoutCandidates)==object@sampleSize)
    if(length(delete)>0){
        nData <- object@nData[-delete,]
        dropoutCandidates <- object@dropoutCandidates[-delete,]
    } else {
        nData <- object@nData
        dropoutCandidates <- object@dropoutCandidates
    }
    
    N <- object@sampleSize
    dropoutRates <- rowSums(dropoutCandidates)/N
    nzmean <- function(x){mean(x[x!=0])}
    averLcpm <- apply(nData*!dropoutCandidates, 1, nzmean)
    qu <- nlsLM(dropoutRates ~ 1/(1+exp(a*(averLcpm-b))),
           start=list(a=1, b=round(median(averLcpm))), trace=FALSE)
    a <- coef(qu)[1]
    b <- coef(qu)[2]
    threshold <- 1/a*log(1/cutoff-1)+b
    x <- 1:10000/10000*max(averLcpm)
    y <- 1/(1+exp(a*(x-b)))
    smoothScatter(averLcpm, dropoutRates, nrpoints=0,
                  xlab="Average of Expressed Entries (logTPM)",
                  ylab="Empirical Dropout Rate",
                  main="Tornado Plot",bty="l")
    points(x, y, col="RED", type="l")
    points(threshold, cutoff, col="RED", pch=16)
    object@wThreshold <- threshold
    return(object)
})

setGeneric("scDissim", function(object, correction=FALSE) {
    standardGeneric("scDissim")
})

#' @title CIDR Dissimilarity Matrix
#'
#' @description
#' computes the \emph{CIDR} dissimilarity matrix.
#'
#' @rdname scDissim
#' @name scDissim
#'
#' @param object the scData class object
#' @param correction boolean; if \code{TRUE} apply Cailliez correction
#'
#' @importFrom Rcpp evalCpp
#' @importFrom ade4 cailliez
#' @export
#'
#' @return an updated scData class object with the following attribute updated
#'
#' \item{dissim}{the \emph{CIDR} dissimilarity matrix}
#'
#' @examples
#' example(cidr)
setMethod("scDissim", "scData", function(object, correction){
      N <- ncol(object@nData)
      Dist <- array(0, dim=c(N, N))
      D <- cpp_dist(Dist, object@dropoutCandidates, object@nData, N, object@wThreshold)
      D <- sqrt(D)
      D <- D+t(D)
      if(correction){
        D <- as.matrix(cailliez(as.dist(D)))
        object@correction <- "Cailliez"
      }
      object@dissim <- D
      return(object)
})

setGeneric("scPCA", function(object) {
    standardGeneric("scPCA")
})

#' @title Single Cell Principal Coordinates Analysis
#'
#' @description
#' performs principal coordinates analysis on the \emph{CIDR} dissimilarity matrix
#'
#' @rdname scPCA
#' @name scPCA
#'
#' @param object the scData class object
#'
#' @export
#'
#' @return an updated scData class object with the following attributes updated
#'
#' \item{eigenvalues}{all eigenvalues (positive and negative) output by the principal coordinates analysis}
#' \item{PC}{principal coordinates}
#' \item{variation}{proportion of variation explained by each of the principal coordinates}
#'
#' @examples
#' example(cidr)
setMethod("scPCA", "scData", function(object){
    y <- cidrPcoa(object@dissim)
    variation <- y$values$Eigenvalues
    object@eigenvalues <- variation
    variation <- variation[variation>0]
    object@PC <- y$vectors[, 1:length(variation)]
    variation <- variation/sum(variation)
    object@variation <- variation
    plot(object@variation, xlab="PC", ylab="Proportion", main="Proportion of Variation")
    return(object)
})

setGeneric("nCluster", function(object, n=NULL, nPC=4) {
    standardGeneric("nCluster")
})

#' @title nCluster
#'
#' @description
#' outputs the Calinski-Harabasz Index versus number of clusters plot.
#'
#' @rdname nCluster
#' @name nCluster
#'
#' @details
#' this method outputs the Calinski-Harabasz Index versus number of clusters plot,
#' useful for deciding the number of clusters.
#' Since \emph{CIDR} automatically decides the number of clusters,
#' this method is optional and only used if the user wants to alter the default number of clusters.
#'
#' @param object the scData class object
#' @param nPC the number of PCs used in clustering; by default 4
#' @param n maximum number of clusters; if \code{NULL} (default), it is set to be nPC*2+2
#'
#' @importFrom clusterCrit intCriteria
#' @export
#'
#' @examples
#' example(cidr)
setMethod("nCluster", "scData", function(object, n, nPC){
    if(is.null(n)){
        n <- nPC*2+2
    }
    exp_clustering <- object@PC[, c(1:nPC)]
    y <- hclust(dist(exp_clustering), method=object@cMethod)
    CH <-rep(NA, n)
    CH[1] <- 0
    for (k in 2:n){
        clusters <- cutree(y, k=k)
        CH[k] <- intCriteria(as.matrix(exp_clustering), clusters, "Calinski_Harabasz")[[1]]
    }
    plot(2:n, CH[2:n], type="b",
         xlab="Number of Clusters", ylab="Calinski-Harabasz Index",
         bty="l")
})

setGeneric("scCluster", function(object, n=NULL, nCluster=NULL, nPC=4) {
    standardGeneric("scCluster")
})

#' @title Sincle Cell Clustering
#'
#' @description
#' performs heirarchical clustering based on \emph{CIDR} principal coordinates.
#'
#' @rdname scCluster
#' @name scCluster
#'
#' @param object the scData class object
#' @param nPC the number of PCs used in clustering; by default 4
#' @param n Calinski-Harabasz Index is used to decide which number between 2 and n is optimal as the number of clusters; if \code{NULL} (default), it is set to be nPC*2+2. User should not assign both n and nCluster.
#' @param nCluster the number of clusters; if \code{NULL} (default), it is determined automatically. User should not assign both n and nCluster.
#' 
#' @importFrom stats hclust
#' @import NbClust
#' @export hclust
#' @export
#'
#' @return an updated scData class object with the following attributes updates
#'
#' \item{nCluster}{the number of clusters}
#' \item{nPC}{the number of PCs used in clustering}
#' \item{clusters}{a vector assigning each cell to a cluster}
#'
#' @examples
#' example(cidr)
setMethod("scCluster", "scData", function(object, n, nCluster, nPC){
    object@nPC <- nPC
    exp_clustering <- object@PC[, c(1:nPC)]
    if (!is.null(n) & !is.null(nCluster)) {
        stop("Invalid input: user should not assign both n and nCluster.")
    } else if (!is.null(n)){
        y <- NbClust(exp_clustering, method=object@cMethod, index="ch", min.nc=1, max.nc = n)
        object@nCluster <- y$Best.nc[1]
        object@clusters <- y$Best.partition
    } else if (!is.null(nCluster)) {
        object@nCluster <- nCluster
        y <- hclust(dist(exp_clustering), method=object@cMethod)
        clusters <- cutree(y, k=object@nCluster)
        object@clusters <- clusters
    } else {
        y <- NbClust(exp_clustering, method=object@cMethod, index="ch", min.nc=1, max.nc=nPC*2+2)
        object@nCluster <- y$Best.nc[1]
        object@clusters <- y$Best.partition
    }
    return(object)
})
