
#' Compute degrees of freedom of one regular expression
#' 
#' @param x, list of vectors, regular expression
#' @noRd
#' @keywords internal
#' @return numeric, degrees of freedom
reg_exp_df <- function(x){
  prod(sapply(x, length))
}


#' Compute degrees of freedom for a list of regular expressions
#' 
#' @param regExp, list of list of vectors, list of regular expressions
#' @return numeric, degrees of freedom 
#' @noRd
#' @keywords internal
#' @examples
#' x1 <- list(c(1,2),3)
#' x2 <- list(1,c(2,4)) 
#' biogram:::count_df(list(x1, x2))
count_df <- function(reg_exp){
  sum(sapply(reg_exp, reg_exp_df))
}


#' Join regular expressions
#' 
#' @param e1 list of vectors, first regular expression
#' @param e2 list of vectors, second regular expression
#' @noRd
#' @keywords internal
#' @return list of vectors, smallest regular expression that contains e1 and e2 
#' @examples 
#' x1 <- list(c(1,2),3)
#' x2 <- list(1,c(2,4)) 
#' biogram:::join_reg_exp(x1, x2)
join_reg_exp <- function(exp1, exp2){
  mapply(FUN = function(x,y) unique(c(x,y)), exp1, exp2, SIMPLIFY=FALSE)
}

#' Clustering of sequences based on regular expression
#' 
#' Clusters sequences hierarchically with regular expressions.
#' At each step we minimize number of degrees of freedom for all regular
#' expressions needed to describe the data
#' 
#' @param ngrams list of elements
#' @export
#' @details Regular expression is a list of the length equal to the length
#' of the input sequences. Each element of the list represents a position in 
#' the sequence and contains amino acid, that are likely to occure on this
#' position.
#' @return List of four
##' \itemize{
##'  \item{"regExps"}{regular expression in best clustering}
##'  \item{"seqClustering"}{clustering of sequences in best clustering}
##'  \item{"allRegExps"}{all regular expressions.}
##'  \item{"allIndices"}{all clusterings}
##' }
##' @examples
##' data(human_cleave)
##' #cluster_reg_exp is computationally expensive
##' \donttest{
##' results <- cluster_reg_exp(human_cleave[1L:10, 1L:4])
##' }
cluster_reg_exp <- function(ngrams){
  n <- length(ngrams)
  regs <- list(NULL)
  regs[[n]] <- ngrams
  indices <- list(NULL)
  indices[[n]] <- lapply(1:n, identity)
  pb <- txtProgressBar(min = 1, max = n-1, style = 3)
  for(i in n:2){
    m <- matrix(ncol=i, nrow=i)
    dfC <- count_df(regs[[i]])
    for(j in 1:(i-1)){
      for(k in (j+1):i){
        m[j,k] <- dfC - count_df(regs[[i]][j]) - count_df(regs[[i]][k]) + count_df(list(join_reg_exp(regs[[i]][[j]], regs[[i]][[k]])))
      }
    }
    ind <- which(m == min(m, na.rm=T), arr.ind = TRUE)[1,]
    regs[[i-1]] <- append(regs[[i]][-ind], list(join_reg_exp(regs[[i]][[ind[1]]], regs[[i]][[ind[2]]])))
    indices[[i-1]] <- append(indices[[i]][-ind], list(unlist(indices[[i]][c(ind[1], ind[2])])))
    names(regs[[i-1]]) <- paste0("cluster_", 1L:(i-1))
    names(indices[[i-1]]) <- paste0("cluster_", 1L:(i-1))
    setTxtProgressBar(pb, n-i+1)
  }
  close(pb)
  plot(sapply(regs, count_df), ylab="Degrees of freedom", xlab="Number of clusters")
  points(which.min(sapply(regs, count_df)), min(sapply(regs, count_df)), col="red", lwd=5)
  clustering <- list(regExps=regs[[which.min(sapply(regs, count_df))]],
                     seqClustering=indices[[which.min(sapply(regs, count_df))]],
                     allRegExps=regs,
                     allIndices=indices)
  class(clustering) = "reg_exp_cluster"
  clustering
}
