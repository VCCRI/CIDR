#' @title Clustering through Imputation and Dimensionality Reduction
#'
#' @description Ultrafast and accurate clustering through imputation and dimensionality
#' reduction for single-cell RNA-Seq (scRNA-Seq) data.
#'
#' @author Peijie Lin <p.lin@victorchang.edu.au>, Michael Troup
#'
#' @docType package
#' @name cidr-package
#' @aliases cidr
#' @useDynLib cidr
#'
#' @examples
#' par(ask=FALSE)
#' ## Generate simulated single-cell RNA-Seq tags.
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
#' sData <- nPC(sData)
#' nCluster(sData)
#' sData <- scCluster(sData)
#' ## Two dimensional visualization: different colors denote different cell types,
#' ## while different plotting symbols denote the clusters output by cidr.
#' plot(sData@PC[,c(1,2)], col=cols,
#'      pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")
#' ## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
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

## class scData - single-cell RNA-Seq data object with attributes relevant to
## clustering through imputation and dimensionality reduction
setClass("scData", representation(tags="matrix",
                                  tagType="character",
                                  sampleSize="numeric",
                                  librarySizes="vector",
                                  nData="matrix",
                                  priorTPM="numeric",
                                  dThreshold="vector",
                                  wThreshold="numeric",
                                  pDropoutCoefA="numeric",
                                  pDropoutCoefB="numeric",
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
         prototype(nCluster=0, correction="none", priorTPM=1))

#' @title scData Constructor
#'
#' @rdname scDataConstructor
#'
#' @description
#' \code{scDataConstructor} creates a new scData class object from a tag table.
#'
#' @details
#' Creates an object in scData (single-cell RNA-Seq dataset) class.
#' Attributes of the class include scalar, vector and matrix
#' data types necessary for the \emph{CIDR} analysis - such as tag table, 
#' library  sizes, dropout candidates, imputation weighting threshold. The 
#' tags can be  raw counts (default) or counts per million (cpm).  Raw counts 
#' are preferrable as the individual library sizes, as determined by the raw 
#' counts, are used to determine dropout candidates.  
#'
#' @param tags a matrix of tags where the rows correspond to features (genes, transcripts, etc) and the columns correspond to cells.
#' @param tagType - \code{"raw"} for when tags are raw counts ; \code{"cpm"} when tags are counts per million ; default is \code{raw}.
#' @export
#' @return an scData class object.
#' @examples
#' ## Generate simulated single-cell RNA-Seq tags.
#' N=3 ## 3 cell types
#' k=50 ## 50 cells per cell type
#' sData <- scSimulator(N=N, k=k)
#' ## The input for cidr should be a tag matrix.
#' ## The default tagType is "raw" - meaning raw counts.
#' tags <- as.matrix(sData$tags)
#' ## create a new scData object
#' sData <- scDataConstructor(tags)
#' ## print the first 5 library sizes
#' sData@librarySizes[1:5]
#' ## print a portion of the data matrix of the class - contains raw tags
#' sData@tags[1:5, 30:34]
#' ## print part of the data matrix of the class - log tag per million
#' sData@nData[1:5, 30:34]
#' 
#' ## Example on using tags that are counts per million (cpm)
#' ## Note that we would only use cpm if we didn't have the raw counts.
#' tags_cpm <- t(t(tags)/colSums(tags))*1000000
#' ## create a new scData object, specifying the tagType parameter
#' sData <- scDataConstructor(tags_cpm, tagType="cpm")
#' ## print the first 5 library sizes
#' ## Note that if only the cpm data is available, we do not know the 
#' ## library sizes.  In this case CIDR sets all the library sizes to
#' ## 1 million.
#' sData@librarySizes[1:5]
#' ## print a portion of the data matrix of the class - contains raw tags
#' sData@tags[1:5, 30:34]
#' ## print part of the data matrix of the class - log tag per million
#' sData@nData[1:5, 30:34]
scDataConstructor <- function(tags, tagType="raw"){
    validTagTypes <- c("raw", "cpm")
    if (!(tagType %in% validTagTypes)) {
        stop("Invalid tagType parameter supplied: ", tagType, ".  Valid Tags: ", 
             paste(validTagTypes, collapse=", ")) 
    }        
    tags <- tags[rowSums(tags)>0,]
    object <- new("scData", tags=tags, tagType=tagType)
    object@sampleSize <- ncol(tags)
    if (tagType=="cpm") {
      object@librarySizes <- rep(1000000, object@sampleSize)
      object@nData <- log2(tags+object@priorTPM)
    } else {
      object@librarySizes <- colSums(tags)
      object@nData <- log2(t(t(tags)/object@librarySizes)*1000000+object@priorTPM)
    }
    return(object)
}

setGeneric("determineDropoutCandidates", function(object, min1=3, min2=8, N=2000, alpha=0.1, fast=TRUE, zerosOnly=FALSE, bw_adjust=1){
    standardGeneric("determineDropoutCandidates")
})

#' @title Determine Dropout Candidates
#'
#' @rdname determineDropoutCandidates
#' @name determineDropoutCandidates
#'
#' @description
#' determines which entries in a single-cell RNA-Seq dataset are dropout candidates.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param object the scData class object.
#' @param min1,min2 technical parameters used in estimating the minimum point between the first two modes of the density curve of logTPM for each cell.
#' @param alpha a cutoff quantile in the range [0,1]. Thresholds outside this will be adjusted to the quantile boundary.
#' @param N number of cells to consider when determining the threshold value for dropout candidates; used in conjunction with the \code{fast} parameter.
#' @param fast Boolean; if \code{TRUE} (default - unless \code{tagType} is \code{cpm}), implements a fast version for datasets with a sample size greater than N. NOTE: set to \code{FALSE} if \code{tagType} is \code{cpm}.
#' @param zerosOnly Boolean; if \code{TRUE}, only zeros are considered as dropout candidates; by default \code{FALSE}.
#' @param bw_adjust bandwidth adjustment factor; \emph{CIDR} uses the default bandwidth selection method ‘nrd0’ in the kernel density estimation; see \code{stats::density} help page for more details. 
#' @export
#' @return an updated scData class object with the following attributes updated
#'
#' \item{dThreshold}{a vector of library dependent dropout candidate thresholds.}
#' \item{dropoutCandidates}{a matrix with the same dimension as nData. The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.}
#' @examples
#' example(cidr)
setMethod("determineDropoutCandidates", "scData", function(object, min1, min2, N, alpha, fast, zerosOnly, bw_adjust){
    if(zerosOnly){
        object@dThreshold <- log2(rep(1, object@sampleSize)/object@librarySizes*1000000+object@priorTPM)
        object@dropoutCandidates <- (object@tags==0)
    } else {
        topLibraries <- 1:object@sampleSize
        ## only use "fast" method if input is raw counts (not cpm)
        if(fast & (object@sampleSize>N) & object@tagType=="raw"){
            topLibraries <- order(object@librarySizes,decreasing = TRUE)[1:N]
        } else {
            N <- object@sampleSize
        }
        dTs <- rep(0,N)
        LT1 <- log2(rep(min1, N)/object@librarySizes[topLibraries]*1000000+object@priorTPM)
        LT2 <- log2(rep(min2, N)/object@librarySizes[topLibraries]*1000000+object@priorTPM)
        object@dropoutCandidates <- array(NA, dim=dim(object@nData))
        for(i in 1:N){
            dfn_m <- density(object@nData[, topLibraries[i]], kernel="epanechnikov", n=1024, from=LT2[i], adjust=bw_adjust)
            dfn_max <- dfn_m$x[which.max(dfn_m$y)]
            dfn <- density(object@nData[, topLibraries[i]], kernel="epanechnikov", n=1024, from=LT1[i], to=dfn_max, adjust=bw_adjust)
            dTs[i] <- dfn$x[which.min(dfn$y)]
        }
        if(fast & (object@sampleSize>N) & object@tagType=="raw"){
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

setGeneric("wThreshold", function(object, cutoff=0.5, plotTornado=FALSE){
    standardGeneric("wThreshold")
})

#' @title Imputation Weighting Threshold
#'
#' @description
#' Determines the imputation weighting threshold.
#'
#' @rdname wThreshold
#' @name wThreshold
#'
#' @param object an scData class object.
#' @param cutoff parameter in the range (0,1), used in the calculation of imputation weighting threshold. Default is 0.5.
#' @param plotTornado Boolean; if \code{TRUE}, the \emph{Tornado Plot} is produced.
#'
#' @details
#' This method finds a function P(u) that maps the average expression level of a gene to the
#' probability of a dropout occurring. It does this by fitting a negative logistic function
#' to the empirical dropouts vs average expression data. The imputation weighting threshold
#' is calculated as the value of u at which P(u) = \code{cutoff}.
#'
#' @importFrom minpack.lm nlsLM
#' @export
#'
#' @return an updated scData class object with the following attribute updated
#'
#' \item{wThreshold}{imputation weighting threshold.}
#' \item{pDropoutCoefA}{the steepness parameter of the negative logistic function that fits the data.}
#' \item{pDropoutCoefB}{the midpoint parameter of the negative logistic function that fits the data.}
#'
#' @examples
#' example(cidr)
setMethod("wThreshold", "scData", function(object, cutoff, plotTornado){
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
    
    if (plotTornado){
        x <- 1:10000/10000*max(averLcpm)
        y <- 1/(1+exp(a*(x-b)))
        smoothScatter(averLcpm, dropoutRates, nrpoints=0,
                      xlab="Average of Expressed Entries (logTPM)",
                      ylab="Empirical Dropout Rate",
                      main="Tornado Plot",bty="l")
        points(x, y, col="RED", type="l")
        points(threshold, cutoff, col="RED", pch=16)
    }
    
    object@wThreshold <- threshold
    object@pDropoutCoefA <- a
    object@pDropoutCoefB <- b
    return(object)
})

setGeneric("scDissim", function(object, correction=FALSE, threads=0, useStepFunction=TRUE) {
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
#' @param object an scData class object.
#' @param correction Boolean; if \code{TRUE}, Cailliez correction is applied; by default \code{FALSE}.
#' @param threads integer; number of threads to be used; by default \code{0}, which uses all available threads.
#' @param useStepFunction Boolean; if \code{TRUE} (default), a step function is used as the imputation weighting function;
#'                                 if \code{FALSE}, a logistic function fitted by \code{\link{wThreshold}} is 
#'                                 used as the imputation weighting function.
#'                                 The logistic function implementation was written by Willem Van Der Byl.
#'                                 
#'
#' @importFrom Rcpp evalCpp
#' @importFrom ade4 cailliez
#' @import RcppParallel
#' @export
#'
#' @return an updated scData class object with the following attribute updated
#'
#' \item{dissim}{the \emph{CIDR} dissimilarity matrix.}
#'
#' @examples
#' example(cidr)
setMethod("scDissim", "scData", function(object, correction, threads, useStepFunction){
      ## the user can choose the number of threads
      threads_int <- as.integer(threads)
      if(!is.na(threads_int) && (threads_int > 0) && (threads_int < defaultNumThreads())) {
          ## user chooses valid thread number - set it
          numThreads <- threads_int
          RcppParallel::setThreadOptions(numThreads=threads_int)
      } else {
          ## reset to default
          numThreads <- defaultNumThreads()
          RcppParallel::setThreadOptions(numThreads=RcppParallel::defaultNumThreads())
      }
      N <- ncol(object@nData)
      Dist <- array(0, dim=c(N, N))
      if (useStepFunction) {
        D <- cpp_dist(Dist, object@dropoutCandidates, object@nData, N, object@wThreshold)
      } else {
        D <- cpp_dist_weighted(Dist, object@dropoutCandidates, object@nData, N, object@pDropoutCoefA, object@pDropoutCoefB)
      }
      D <- sqrt(D)
      D <- D+t(D)
      if(correction){
        D <- as.matrix(cailliez(as.dist(D)))
        object@correction <- "Cailliez"
      }
      object@dissim <- D
      return(object)
})

setGeneric("scPCA", function(object, plotPC=TRUE) {
    standardGeneric("scPCA")
})

#' @title Single-cell Principal Coordinates Analysis
#'
#' @description
#' performs principal coordinates analysis on the \emph{CIDR} dissimilarity matrix.
#'
#' @rdname scPCA
#' @name scPCA
#'
#' @param object an scData class object.
#' @param plotPC Boolean; if \code{TRUE}, a plot of PC variance explained is produced.
#'
#' @export
#'
#' @return an updated scData class object with the following attributes updated
#'
#' \item{eigenvalues}{all eigenvalues (positive and negative) output by the principal coordinates analysis.}
#' \item{PC}{principal coordinates.}
#' \item{variation}{proportion of variation explained by each of the principal coordinates.}
#'
#' @examples
#' example(cidr)
setMethod("scPCA", "scData", function(object,plotPC){
    y <- cidrPcoa(object@dissim)
    variation <- y$values
    ## store all eigenvalues - neg, 0, & pos
    object@eigenvalues <- variation
    ## for variation, only deal with positive eigenvalues
    variation <- variation[variation>0]
    object@PC <- y$vectors[, 1:length(variation)]
    variation <- variation/sum(variation)
    object@variation <- variation

    if(plotPC) plot(object@variation, xlab="PC", ylab="Proportion", main="Proportion of Variation")
    return(object)
})


setGeneric("nPC", function(object) {
    standardGeneric("nPC")
})

#' @title Determine nPC
#'
#' @description
#' determines the optimal number of principal coordinates (nPC) to be used in clustering.
#'
#' @rdname nPC
#' @name nPC
#'
#' @param object an scData class object.
#'
#' @export
#'
#' @return an updated scData class object with the following attribute updated
#'
#' \item{nPC}{the number of principal coordinates to be used in clustering.}
#'
#' @examples
#' example(cidr)
setMethod("nPC", "scData", function(object){
    object@nPC <- calc_npc(object@variation)
    return(object)
})


setGeneric("nCluster", function(object, n=NULL, nPC=NULL, cMethod="ward.D2") {
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
#' outputs the Calinski-Harabasz Index versus number of clusters plot,
#' useful for deciding the number of clusters.
#' Since \emph{CIDR} automatically decides the number of clusters,
#' this method is optional and only used if the user wants to alter the default number of clusters.
#'
#' @param object an scData class object.
#' @param nPC the number of PCs used in clustering; by default 4.
#' @param n maximum number of clusters; if \code{NULL} (default), it is set to be nPC*2+2.
#' @param cMethod hierarchical clustering method; by default "ward.D2".
#'
#' @importFrom clusterCrit intCriteria
#' @export
#'
#' @examples
#' example(cidr)
setMethod("nCluster", "scData", function(object, n, nPC, cMethod){
    if(!is.null(nPC)){
        object@nPC <- nPC
    } else {
        nPC <- object@nPC
    }
    if(is.null(n)){
        n <- nPC*3+3
    }

    exp_clustering <- object@PC[, c(1:nPC)]
    CH <- NbClust(exp_clustering, method=cMethod, index="ch", min.nc=1, max.nc=n)$All.index
    
    plot(1:n, CH[1:n], type="b",
         xlab="Number of Clusters", ylab="Calinski-Harabasz Index",
         bty="l")
})

setGeneric("scCluster", function(object, n=NULL, nCluster=NULL, nPC=NULL, cMethod="ward.D2") {
    standardGeneric("scCluster")
})

#' @title Single-cell Clustering
#'
#' @description
#' performs heirarchical clustering on \emph{CIDR} principal coordinates.
#'
#' @rdname scCluster
#' @name scCluster
#'
#' @param object an scData class object.
#' @param nPC the number of PCs used in clustering; by default 4.
#' @param n Calinski-Harabasz Index is used to decide which number between 2 and n is optimal as the number of clusters; if \code{NULL} (default), it is set to be nPC*2+2. User should not assign both n and nCluster.
#' @param nCluster the number of clusters; if \code{NULL} (default), it is determined automatically. User should not assign both n and nCluster.
#' @param cMethod hierarchical clustering method; by default "ward.D2".
#' 
#' @importFrom stats hclust
#' @import NbClust
#' @export hclust
#' @export
#'
#' @return an updated scData class object with the following attributes updated
#'
#' \item{nCluster}{the number of clusters.}
#' \item{nPC}{the number of PCs used in clustering.}
#' \item{clusters}{a vector assigning each cell to a cluster.}
#' \item{cMethod}{hierarchical clustering method.}
#'
#' @examples
#' example(cidr)
setMethod("scCluster", "scData", function(object, n, nCluster, nPC, cMethod){
    if(!is.null(nPC)){
        object@nPC <- nPC
    } else {
        nPC <- object@nPC
    }
    object@cMethod <- cMethod

    exp_clustering <- object@PC[, c(1:nPC)]
    if (!is.null(n) & !is.null(nCluster)) {
        stop("Invalid input: user should not assign both n and nCluster.")
    } else if (!is.null(nCluster)) {
        object@nCluster <- nCluster
        object@clusters <- NbClust(exp_clustering, method=object@cMethod, index="ch", min.nc=object@nCluster, max.nc=object@nCluster)$Best.partition
    } else {
        if (is.null(n)) {
            n <- nPC*3+3
        }
        CH <- NbClust(exp_clustering, method=object@cMethod, index="ch", min.nc=1, max.nc=n)$All.index
        l <- length(CH)
        a <- as.vector(CH[-c(1,l-1,l)]+CH[-c(1:3)] - 2*CH[-c(1,2,l)])
        b <- which.min(a)
        c <- which.min(a[-c(1:b)])
        if ((3*a[b+c])<a[b]){
            object@nCluster <- b+c+2
        } else {
            object@nCluster <- b+2
        }
        
        object@clusters <- NbClust(exp_clustering, method=object@cMethod, index="ch", min.nc=object@nCluster, max.nc=object@nCluster)$Best.partition
    } 
    return(object)
})
