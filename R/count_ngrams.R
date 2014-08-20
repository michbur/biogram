#' Detect and count n-grams in single sequence
#'
#' Counts all n-grams present in given sequence.
#'
#' @param seq \code{integer} vector or matrix describing sequence(s). 
#' @param n size of n-grams.
#' @param u unigrams (\code{integer}, \code{numeric} or \code{character} vector).
#' @param d distance between elements of n-gram (0 means consecutive elements). See
#' Details.
#' @param pos \code{logical}, if \code{TRUE} n_grams contains position information.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @param threshold \code{integer}, if not equal to 0, data is binarized into
#' two groups (larger or equal to threshold, smaller than threshold).
#' @return a named \code{integer} vector. Elements of n-gram are separated by dot.
#' @note List of possible n-grams must be calculated outside of the function.
#' @details The length of \code{distance} vector should be always \code{n} - 1. For example 
#' when \code{n} = 3, \code{d} = c(1, 2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2, 0, 1) means A__AA_A. If vector \code{d} has length 1, it is recycled to
#' length \code{n} - 1.
#' @note does not work: n = 1 and dist != 0
#' @export
#' @seealso 
#' Create vector of possible n-grams: \code{\link{create_ngrams}}.
#' Get indices of n-grams: \code{\link{get_ngrams_ind}}.
#' @examples 
#' #trigrams for nucleotides
#' count_ngrams(sample(1L:4, 50, replace = TRUE), 3, 1L:4, pos = TRUE)
#' #multigrams
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_ngrams(seqs, c(1, 2, 3), 1L:4, pos = TRUE)

count_ngrams <- function(seq, n, u, d = 0, pos = FALSE, scale = FALSE, threshold = 0) {
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  len_seq <- ncol(seq)
  
  do.call(cbind, lapply(n, function(current_n) {
    #create list of n-grams for given n
    possib_ngrams <- create_ngrams(current_n, u)
    
    do.call(cbind, lapply(d, function(current_d) {
      #look for n-gram indices for given d
      ngram_ind <- get_ngrams_ind(len_seq, current_n, current_d)
      
      #use attr(ngram_ind, "d") instead of current_d because of distance recycling
      max_grams <- len_seq - current_n - sum(attr(ngram_ind, "d")) + 1
      grams <- apply(seq, 1, function(i)
        seq2ngrams(i, ind = ngram_ind, max_grams))
      ncol_grams <- ncol(grams)
      if (pos) {
        n_grams_number <- length(create_ngrams(current_n, u, max_grams))
        res <- t(vapply(1L:ncol(grams), function(ngram_column)
          unlist(lapply(possib_ngrams, function(current_ngram)
            grams[, ngram_column] == current_ngram)), rep(0, n_grams_number)))
        
        colnames(res) <- create_ngrams(current_n, u, max_grams)
        res
      } else {
        res <- vapply(possib_ngrams, function(current_ngram)
          vapply(1L:ncol_grams, function(ngram_column)
            sum(grams[, ngram_column] == current_ngram), 0), rep(0, ncol_grams))
      }
      
      if (threshold > 0) {
        ind_sums <- rowSums(res)
        res <- res[ind_sums >= threshold, ]
      }
      
      if (scale)
        res <- res/max_grams
      
      res
    }))
  }))
}

count_ngrams_helper <- function(seq, feature_list, n, ind, pos) {
  #TO DO implement position
  #feature list(list of possible n-grams) is outside, because count_ngrams is meant to
  #be inside the loop
  #same for indices
  if (n > 1) {
    element_matrix <- do.call(cbind, lapply(ind, function(i) seq[i]))
    grams <- apply(element_matrix, 1, function(x) 
      paste(x, collapse="."))
  } else {
    grams <- seq
  }
  
  if (pos) {
    grams <- paste(1L:length(grams), grams, sep = "_")
  } 
  
  sapply(feature_list, function(i)
    sum(i == grams))
}

seq2ngrams <- function(seq, ind, max_grams) {
  #TO DO implement position
  #feature list(list of possible n-grams) is outside, because count_ngrams is meant to
  #be inside the loop
  #same for indices
  if (length(ind) > 1) {
    element_matrix <- vapply(ind, function(i) seq[i], rep(0, max_grams))
    grams <- apply(element_matrix, 1, function(x) 
      paste(x, collapse="."))
  } else {
    grams <- seq
  }
  
  grams
}
