#' Count specified n-grams
#'
#' Counts specified n-grams present in the input sequence(s).
#'
#' @param seq \code{integer} vector or matrix describing sequence(s). 
#' @param ngrams a vector of n-grams. Must have the same \code{n}.
#' @return a \code{\link[slam]{simple_triplet_matrix}} where columns represent
#' n-grams and rows sequences.
#' @export
#' @examples
#' seqs <- matrix(c(1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 3, 4, 1, 2, 2, 4), nrow = 2)
#' count_specified(seqs, ngrams = c("2_1.1.1_0.0", "2_2.2.2_0.0", "3_1.1.2_0.0"))

count_specified <- function(seq, ngrams) {
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
    
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  df <- ngrams2df(ngrams)
  
  #splitted ngrams
  sn_grams <- strsplit(df[, "ngram"], ".", fixed = TRUE)
  
  #n in ngram
  n <- unique(vapply(sn_grams, length, 1))
  if(length(n) > 1)
    stop("n-grams must have the same n.")
  
  distances <- lapply(strsplit(df[, "distance"], ".", fixed = TRUE), as.numeric)
  
  res <- vapply(1L:nrow(df), function(ngram_id)
    vapply(1L:n_seqs, function(single_seq) {
      #all possible n-gram positions
      all_ngram_pos <- get_ngrams_ind(len_seq, n, distances[[ngram_id]])
      #positions of the n-gram of interest
      single_ngram_pos <- sapply(all_ngram_pos, function(single_pos) 
        single_pos[df[ngram_id, "position"]])
      as.numeric(all(as.character(seq[single_seq, single_ngram_pos]) == sn_grams[[ngram_id]]))
    }, 0), rep(0, n_seqs))
  colnames(res) <- ngrams
  as.simple_triplet_matrix(res)
}