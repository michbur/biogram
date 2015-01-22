#' Detect and count multiple n-grams in sequences
#'
#' A convinient wrapper around \code{\link{count_ngrams}} for counting multiple
#' values of \code{n} and \code{d}.
#'
#' @inheritParams count_ngrams
#' @param ns \code{numeric} vector of n-grams' sizes. See Details.
#' @param ds \code{list} of distances between elements of n-grams. Each element of the list
#' is a vector used as distance for the respective n-gram size given by the \code{ns}
#' parameter.
#' @return a \code{integer} matrix with named columns. The naming conventions are the same
#' as in \code{\link{count_ngrams}}.
#' @export
#' @examples 
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_multigrams(c(3, 1), list(c(1, 0), 0), seqs, 1L:4, pos = TRUE)

count_multigrams <- function(ns, ds, seq, u, pos = FALSE, scale = FALSE, threshold = 0) {
  n_loops <- length(ns)
  do.call(cbind, lapply(1L:n_loops, function(current_loop) {
      count_ngrams(seq, ns[current_loop], u, ds[[current_loop]], 
                   pos, scale, threshold)
  }))
}


