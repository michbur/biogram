#' Detect and count n-grams in sequences
#'
#' Counts all n-grams present in sequences.
#'
#' @param seq \code{integer} vector or matrix describing sequence(s). 
#' @param n \code{integer} size of n-gram.
#' @param u unigrams (\code{integer}, \code{numeric} or \code{character} vector).
#' @param d \code{integer} vector of distances between elements of n-gram (0 means 
#' consecutive elements). See Details.
#' @param pos \code{logical}, if \code{TRUE} n_grams contains position information.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @param threshold \code{integer}, if not equal to 0, data is binarized into
#' two groups (larger or equal to threshold, smaller than threshold).
#' @return a \code{integer} matrix with named columns. Elements of n-gram are separated 
#' by dot. If \code{pos}, the left side of name means actual position of the n-gram 
#' (separated by \code{_}). the Right side of name is vector of distance(s) used separated by
#' \code{_}. See \code{Note} for examples.
#' @note List of possible n-grams must be calculated outside of the function.
#' @details The length of \code{distance} vector should be always \code{n} - 1. For example 
#' when \code{n} = 3, \code{d} = c(1, 2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2, 0, 1) means A__AA_A. If vector \code{d} has length 1, it is recycled to
#' length \code{n} - 1.
#' @note 46_4.4.4_0_1 means n-gram 44_4 on position 46.
#' 12_2.1_2 means n-gram 2__1 on position 12.
#' @export
#' @seealso 
#' Create vector of possible n-grams: \code{\link{create_ngrams}}.
#' Get n-grams from analyzed sequence: \code{\link{seq2ngrams}}.
#' Get indices of n-grams: \code{\link{get_ngrams_ind}}.
#' @examples 
#' #trigrams for nucleotides
#' count_ngrams(sample(1L:4, 50, replace = TRUE), 3, 1L:4, pos = TRUE)
#' #trigrams from multiple sequences
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_ngrams(seqs, 3, 1L:4, pos = TRUE)

count_ngrams <- function(seq, n, u, d = 0, pos = FALSE, scale = FALSE, threshold = 0) {
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  #create list of n-grams for n
  possib_ngrams <- create_ngrams(n, u)
  
  #look for n-gram indices for d
  ngram_ind <- get_ngrams_ind(len_seq, n, d)
  
  #use attr(ngram_ind, "d") instead of d because of distance recycling
  max_grams <- len_seq - n - sum(attr(ngram_ind, "d")) + 1
  
  #extract n-grams from sequene
  grams <- vapply(1L:n_seqs, function(i)
    seq2ngrams_helper(seq[i, ], ind = ngram_ind, max_grams), rep("a", max_grams))
  
  if (pos) {
    #get positioned possible n-grams
    pos_possib_ngrams <- create_ngrams(n, u, max_grams)
    n_pos_possib_ngrams <- length(pos_possib_ngrams)
    
    #it turned out to be faster than matrix substitution
    res <- t(vapply(1L:n_seqs, function(current_sequence)
      vapply(possib_ngrams, function(current_ngram)
        grams[, current_sequence] == current_ngram, rep(0, max_grams)), 
      rep(0, n_pos_possib_ngrams)))
    colnames(res) <- pos_possib_ngrams
    
    res
  } else {
    res <- t(vapply(possib_ngrams, function(current_ngram)
      vapply(1L:n_seqs, function(current_sequence)
        sum(grams[, current_sequence] == current_ngram), 0), rep(0, n_seqs)))
  }
  
  if (threshold > 0) {
    ind_sums <- rowSums(res)
    res <- res[ind_sums >= threshold, ]
  }
  
  if (scale)
    res <- res/max_grams
  colnames(res) <- paste(colnames(res), paste0(attr(ngram_ind, "d"), collapse = "_"), 
                         sep = "_")
  res
}


count_ngrams_helper <- function(seq, feature_list, n, ind, pos) {
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

#ultrafast function for n-gram extraction
seq2ngrams_helper <- function(seq, ind, max_grams) {
  #get all consecutive n-grams from sequence
  #length(ind) > 1 - equivalent of n > 1 
  if (length(ind) > 1) {
    #element_matrix contains elements of n-gram in matrix structure
    element_matrix <- vapply(ind, function(i) seq[i], rep(0, max_grams))
    grams <- vapply(1L:nrow(element_matrix), function(ith_row)
      paste(element_matrix[ith_row, ], collapse="."), "a")
  } else {
    #in case of unigrams take sequence - output must be character
    grams <- as.character(seq)
  }
  
  #all n-grams from seuqence - character vector
  grams
}

