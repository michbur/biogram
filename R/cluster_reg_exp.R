
#' Compute degrees of freedom of one regular expression
#' 
#' @param x, list of vectors, regular expression
#' @keywords internal
#' @return numeric, degrees of freedom
reg_exp_df <- function(x){
  prod(sapply(x,  length))
}


#' Compute degrees of freedom for a list of regular expressions
#' 
#' @param regExp, list of list of vectors, lsit of regular expressions
#' @keywords internal
#' @return numeric, degrees of freedom 
#' @examples
#' x1 <- list(c(1,2),3)
#' x2 <- list(1,c(2,4)) 
#' count_df(list(x1, x2))
count_df <- function(reg_exp){
  sum(sapply(reg_exp, reg_exp_df))
}


#' Join regular expressions
#' 
#' @param e1 list of vectors, first regular expression
#' @param e2 list of vectors, second regular expression
#' @keywords internal
#' @return list of vectors, samllest regular expression that contains e1 and e2 
#' @examples 
#' x1 <- list(c(1,2),3)
#' x2 <- list(1,c(2,4)) 
#' join_reg_exp(x1, x2)
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
#' @return List of four
##' \itemize{
##'  \item{"regExps"}{regular expression in best clustering}
##'  \item{"seqClustering"}{clustering of sequences in best clustering}
##'  \item{"allRegExps"=regs"}{all regular expressions}
##'  \item{"allIndices"}{all clusterings}
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
    regs[[i-1]] <- append(regs[[i]][-ind], list(joinRegExp(regs[[i]][[ind[1]]], regs[[i]][[ind[2]]])))
    indices[[i-1]] <- append(indices[[i]][-ind], list(unlist(indices[[i]][c(ind[1], ind[2])])))
    names(regs[[i-1]]) <- paste0("cluster_", 1:(i-1))
    names(indices[[i-1]]) <- paste0("cluster_", 1:(i-1))
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
