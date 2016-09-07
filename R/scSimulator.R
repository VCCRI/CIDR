## library size 1 million

## function for simulating RNA-Seq Data

#' @title Single-cell RNA-Seq Tags Simulator
#'
#' @rdname scSimulator
#' @name scSimulator
#'
#' @description
#' simulates a table of tags generated from a single-cell RNA-Seq experiment.
#'
##@details
#'
#'
#' @param nDG the number of differentially expressed genes.
#' @param nMK the number of markers (expressed in one cell type only).
#' @param nNDG the number of non-differentially expressed genes.
#' @param k the number of cells in each cell type.
#' @param seed for reproducibility.
#' @param logmean mean of distrubution for each library.
#' @param logsd standard deviation of distribution for each library.
#' @param v dropout level parameter; higher v means a higher level of dropouts.
#' @export
#' @return a list object is returned containing the following components
#'
#' \item{annotation}{0 means a real zero; 1 means a non-zero entry; 2 means a dropout.  Note that a dropout may still have a small value.}
#' \item{expectedValues}{underlying real expression of the entry.}
#' \item{dropoutsSimulated}{expression values after simulation of dropouts before simulation of noise.}
#' \item{tags}{the simulated tags after simulation of both dropouts and noise.}
#'
#'
#' @examples
#' ## Generate simulated tags with default parameters: 3 cell types and 50 cells in each cell type.
#' sData <- scSimulator()
#' sData$tags[1:5, 1:5]
#'
scSimulator <- function(N=3, nDG=150, nMK=10, nNDG=20000, k=50,
                         seed=17,logmean=5.25,logsd=1, v=9.2){
    set.seed(seed)

    scmdSimulator <- function(N,logmean,logsd,nDG){
    trD <- array(0,dim=c(nDG,N))
    for (i in 1:N){
      trD[,i] <- rlnorm(nDG,meanlog = logmean,sdlog= logsd) -1
    }
    return(trD)
    }

    ## Simulating Differentially Expressed Genes

    DG <- scmdSimulator(N,logmean,logsd,nDG)

    ## Simulating Marker Genes
    mkSimulator <- function(N,logmean,logsd,nMK){
    y <- array(0,dim=c(nMK*N,N))
    for (i in 1:N){
      y[(i*nMK-nMK+1):(i*nMK) ,i] <- rlnorm(nMK,meanlog = logmean,sdlog = logsd) -1
    }
    return(y)
    }

    MK <- mkSimulator(N,logmean,logsd,nMK)

    DG <- rbind(DG,MK)

    ## Simulating the non-differentially expressed genes
    NDG0 <-  rlnorm(nNDG,meanlog = logmean,sdlog=logsd) -1

    NDG <- array(NA,dim=c(nNDG,N))
    for(i in 1:N){
    NDG[,i] <- NDG0
    }

    ## Expected Values Matrix
    EM <- rbind(DG,NDG)

    ## Simulating noise and drop-outs

    design <- rep(k,N)

    u <-0.7

    Simulator <- function(EM,design){
    a<-list() ## annotation
    b<-list() ## expression matrix - expected values
    c<-list() ## expression matrix - dropouts simulated
    d<-list() ## tags matrix
    TZ <- nrow(EM)
    for (i in 1:ncol(EM)){
      a[[i]]<- array(0,dim=c(TZ,design[i]))
      b[[i]]<- array(NA,dim=c(TZ,design[i]))
      c[[i]]<- array(NA,dim=c(TZ,design[i]))
      d[[i]]<- array(NA,dim=c(TZ,design[i]))

      a[[i]][EM[,i]>0,] <- 1

      ## probability function
      pi <- 1/(1 + exp(u * (log2(EM[,i]+1) - v)))
      for (j in 1:design[i]){
        b[[i]][,j] <- EM[,i]
        c[[i]][,j] <- EM[,i]
        for (k in 1:TZ){
          if(rbinom(1,1,pi[k])){
            c[[i]][k,j] <- max(rpois(1,1)-1,0)
            a[[i]][k,j] <- 2 * a[[i]][k,j]
          }
        }

        ## simulating noise
        nj<-rpois(TZ,lambda=round(c[[i]][,j]))*(a[[i]][,j]==1)
        nj <- nj+c[[i]][,j]*(a[[i]][,j]==2)
        nj[nj<0]<-0
        d[[i]][,j]<-nj
      }
    }

    A <- do.call(cbind,a)
    B <- do.call(cbind,b)
    C <- do.call(cbind,c)
    D <- do.call(cbind,d)

    G <-list(A,B,C,D)
    names(G) <-c("annotation","expectedValues","dropoutsSimulated","tags")
    return(G)
    }

    sData <- Simulator(EM,design)

    return(sData)
}
