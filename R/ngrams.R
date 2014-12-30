#' Get All Possible N-Grams
#'
#' Creates vector of all posible n_grams (for given \code{n}).
#'
#' @inheritParams count_ngrams
#' @param possible_grams number of possible n-grams. If not \code{NULL} n-grams do not
#' contain information about position
#' @return a character vector. Elements of n-gram are separated by dot.
#' @note Input data must be a matrix or data frame of numeric elements.
#' @details See Details section of \code{\link{count_ngrams}} for more 
#' information about n-grams naming convention. The possible information about distance 
#' must be added by hand (see examples).
#' @export
#' @examples 
#' #bigrams for standard aminoacids
#' create_ngrams(2, 1L:20)
#' #bigrams for standard aminoacids with positions, 10 amino acid long sequence, so 
#' #only 9 bigrams can be located in sequence
#' create_ngrams(2, 1L:20, 9)
#' #bigrams for DNA with positions, 10 nucleotide long sequence, distance 1, so only 8 bigrams
#' #in sequence
#' #paste0 adds information about distance at the end of n-gram
#' paste0(create_ngrams(2, 1L:4, 8), "_0")


create_ngrams <- function(n, u, possible_grams = NULL) {
  grid_list <- lapply(1L:n, function(i) u)
  res <- apply(expand.grid(grid_list), 1, function(x)
    paste(x, collapse = "."))
  if (!is.null(possible_grams))
    res <- as.vector(sapply(res, function(i) paste(1L:possible_grams, i, sep = "_")))
  res
}


#' Get Indices of N-Grams
#'
#' Computes list of n-gram elements positions in sequence.
#'
#' @param len_seq \code{integer} value describing sequence's length.
#' @inheritParams count_ngrams
#' @return A list with number of elements equal to \code{n}. Every element is a 
#' vector containing locations of given n-gram letter. For example, first element of
#' list contain indices of first letter of all n-grams. The attribute \code{d} of output
#' contains distances between letter used to compute locations (see Details).
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_ngrams}}.
#' @export
#' @examples 
#' #positions trigrams in sequence of length 10
#' get_ngrams_ind(10, 9, 0)

get_ngrams_ind <- function(len_seq, n, d) {
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive letter (a vector of distances)
  
  #calculate indices of n-grams elements
  ind <- lapply(1L:n, function(i) 
    (1 + i - 1):(len_seq - n + i))
  
  if(length(d) != 1 && length(d) != n - 1)
    stop("Length of d must be 1 or n - 1")
  
  
  if(n > 1) {
    #if distance vector is too short, recycle it
    if(length(d) == 1 && n > 2)
      d <- rep(d, n - 1)
    
    if(sum(d) > 0) {
      ind[-1] <- lapply(1L:length(d), function(i)
        ind[[i + 1]] + sum(d[1L:i]))
      not_taken <- ind[[1]][(length(ind[[1]]) - sum(d) + 1):length(ind[[1]])]
      ind <- lapply(ind, function(i) i[-not_taken])
    }
    
    attr(ind, "d") <- d
  } else {
    #distance is a nonsense for unigrams
    attr(ind, "d") <- 0
  }
  
  ind
}


#' Extract N-Grams From Sequence
#'
#' Extracts vector of n-grams present in sequence(s).
#'
#' @inheritParams count_ngrams
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_ngrams}}.
#' @return A \code{character} matrix of n-grams, where every row corresponds to a
#' different sequence.
#' @export
#' @examples 
#' #trigrams from multiple sequences
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' seq2ngrams(seqs, 3, 1L:4)

seq2ngrams <- function(seq, n, u, d = 0) {
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  #look for n-gram indices for d
  ngram_ind <- get_ngrams_ind(len_seq, n, d)
  
  max_grams <- calc_max_grams(len_seq, n, ngram_ind)
  
  #extract n-grams from sequene
  res <- t(vapply(1L:n_seqs, function(i) {
    grams <- seq2ngrams_helper(seq[i, ], ind = ngram_ind, max_grams)
    paste(grams, paste0(attr(ngram_ind, "d"), collapse = "_"), 
          sep = "_")
  }, rep("a", max_grams)))
  if (max_grams == 1)
    res <- t(res)
  res
}

