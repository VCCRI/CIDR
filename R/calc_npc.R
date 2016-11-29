#' @title nPC Calculation
#'
#' @rdname calc_npc
#' @name calc_npc
#'
#' @description calculates nPC.
#'
#' @details
#' This function is used internally in the method \code{nPC};
#' see the help page for \code{nPC} for more details.
#'
#'
#' @param var a vector; proportion of variation explainded by each of the principal coordinates.
#'
#' @return This function is used internally in the method \code{nPC};
#' see the help page for \code{nPC} for more details.
#'
#'
#' @export calc_npc
#' 
#'

calc_npc <- function(var, N=1, cutoff_divisor=10) {
    NPC_DEFAULT <- 4
    d <- var[-length(var)] - var[-1]
    descending_d <- sort(d, decreasing=T)
    max_d <- descending_d[1]
    mean_N_max <- mean(descending_d[1:N])
    # get measure of spread of data as a proportion of the largest 3 differences
    spread <- mean(d)/mean(descending_d[1:3])*100
    #print(paste0("spread = ", spread))
    if (spread > 15) {
        ## data too evenly spread ?
        return(NPC_DEFAULT)
    } else if (spread > 10) {
        cutoff_divisor <- 5
    }
    cutoff <- mean_N_max/cutoff_divisor
    #print(paste0("cutoff: ", cutoff))
    groups <- list()
    groups[[1]] <- c(1)
    group_index <- 1
    for (i in 2:length(var)) {
        if (d[i-1] < cutoff) {
            #print(d[i-1])
            # include in current group
            groups[[group_index]] <- c(groups[[group_index]], i)
            ## check stopping criteria
            len_curr_group <- length(groups[[group_index]])
            if (len_curr_group > 7 || (i > 10 && len_curr_group > 3)) {
                #if (len_curr_group > 10) {
                ## return the last index in the previous group
                prev_group <- groups[[group_index - 1]]
                nPC <- prev_group[length(prev_group)]
                ## not helpful to return nPC=1, so use default in this case
                if (nPC == 1) {
                    nPC <- NPC_DEFAULT
                }
                return(nPC)
            }
        } else {
            # start new group
            group_index <- group_index + 1
            groups[[group_index]] <- i
        } 
    }
    ## must not have converged - return default
    return(NPC_DEFAULT) 
}
