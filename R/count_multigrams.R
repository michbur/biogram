#' Detect and count multiple n-grams in sequences
#'
#' A convinient wrapper around \code{\link{count_ngrams}} for counting multiple
#' values of \code{n} and \code{d}.
#'
#' @inheritParams count_ngrams
#' @param n_d \code{list} list of n-grams' sizes and distances between elements of n-gram.
#' See Details.
#' @return a \code{integer} matrix with named columns. The naming conventions are the same
#' as in \code{\link{count_ngrams}}. 
#' @details Each element of \code{n_d} is a list consisting of two vectors. 
#' First element is a single \code{integer} value which determines the number of 
#' words in n-gram (equivalent of \code{n} from \code{\link{count_ngrams}}). Second 
#' element must be an \code{integer} vector describing distances between words in 
#' n-gram (equivalent of \code{d} from \code{\link{count_ngrams}}).
#' @export
#' @examples 
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_multigrams(list(list(2, 1), list(2, 3)), seqs, 1L:4, pos = TRUE)

count_multigrams <- function(n_d, seq, u, d = 0, pos = FALSE, scale = FALSE, threshold = 0) {
  n_loops <- length(n_d)
  do.call(cbind, lapply(1L:n_loops, function(current_loop) {
      count_ngrams(seq, n_d[[current_loop]][[1]], u, n_d[[current_loop]][[2]], 
                   pos, scale, threshold)
  }))
}


