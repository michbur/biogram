#' Detect and count n-grams in single sequence
#'
#' Counts all n-grams present in given sequence.
#'
#' @param seq \code{integer} vector describing sequence.
#' @param n size of n-grams.
#' @param u unigrams (\code{integer}, \code{numeric} or \code{character} vector).
#' @param d distance between elements of n-gram (0 means consecutive elements). See
#' Details.
#' @param pos \code{logical}, if \code{TRUE} n_grams contains position information.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @return a named \code{integer} vector. Elements of n-gram are separated by dot.
#' @note List of possible n-grams must be calculated outside of the function.
#' @details The length of \code{distance} vector should be always n - 1. For example 
#' when \code{n} = 3, \code{d} = c(1, 2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2, 0, 1) means A__AA_A.
#' @note does not work: n = 1 and dist != 0
#' @export
#' @seealso For convenient wrapper for multidistance calculations 
#' see \code{\link{gramize_data}}.
#' @examples 
#' #trigrams for nucleotides
#' count_ngrams(sample(1L:4, 50, replace = TRUE), 3, 1L:4, pos = TRUE)
#' #multigrams
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_ngrams(seqs, c(1, 2, 3), 1L:4, pos = TRUE)

count_ngrams <- function(seq, n, u, d = 0, pos = FALSE, scale = FALSE) {
  
  if (class(seq) == "matrix") {
    len_seq <- ncol(seq)
  } else {
    len_seq <- length(seq)
  }
  
  do.call(cbind, lapply(n, function(current_n) {
    possib_ngrams <- create_ngrams(current_n, u)
    t(do.call(cbind, lapply(d, function(current_d) {
      max_grams <- len_seq - current_n - sum(current_d) + 1
      ngram_ind <- get_ngrams_ind(len_seq, current_n, current_d)
      
      if (pos)
        possib_ngrams <- create_ngrams(current_n, u, max_grams)
      
      if (class(seq) == "matrix") {
        res <- apply(seq, 1, function(i)
          count_ngrams_helper(i, feature_list = possib_ngrams, current_n, 
                              ind = ngram_ind, pos))
      } else {
        res <- count_ngrams_helper(seq, feature_list = possib_ngrams, current_n, 
                            ind = ngram_ind, pos)
      }
      
      names(res) <- paste0(names(res), "_", paste(attr(ngram_ind, "d"), collapse = "_"))
      
      if (scale)
        res <- res/max_grams
      
      res
    })))
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


#' Get all possible n-grams
#'
#' Creates vector of all posible n_grams.
#'
#' @inheritParams count_ngrams
#' @param possible_grams number of possible n-grams. If not \code{NULL} n-grams do not
#' contain information about position
#' @return a character vector. Elements of n-gram are separated by dot.
#' @note Input data must be a matrix or data frame of numeric elements.
#' @details N-gram coding: \code{p_1.2.1} means position \code{p}, \code{1.2.1} means 
#' unigrams constituing n-gram.  
#' @export
#' @examples 
#' #bigrams for standard aminoacids
#' create_ngrams(2, 1L:20)
#' #bigrams for standard aminoacids with positions, 10 nucleotide long sequence
#' create_ngrams(2, 1L:20, 9)

create_ngrams <- function(n, u, possible_grams = NULL) {
  grid_list <- lapply(1L:n, function(i) u)
  res <- apply(expand.grid(grid_list), 1, function(x)
    paste(x, collapse = "."))
  if (!is.null(possible_grams))
    res <- as.vector(sapply(res, function(i) paste(1L:possible_grams, i, sep = "_")))
  res
}


#' Get indices of n-grams
#'
#' Computes list of n-gram elements positions in sequence.
#'
#' @param len_seq \code{integer} value describing sequence's length.
#' @inheritParams count_ngrams
#' @return A list with number of elements equal to \code{n}. Every element is a 
#' vector containing locations of given n-gram letter. For example, first element of
#' list contain indices of first letter of all n-grams. The attribute \code{d} of output
#' contains distances between letter used to compute locations (see Details).
#' @details The length of \code{d} vector should be always n - 1. For example 
#' when \code{n} = 3, \code{d} = c(1, 2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2, 0, 1) means A__AA_A. Shorter \code{d} vectors will be recycled.
#' @export
#' @examples 
#' #trigrams in sequence of length 10
#' get_ngrams_ind(10, 9, 0)

get_ngrams_ind <- function(len_seq, n, d) {
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive letter (a vector of distances)
  
  #calculate indices of n-grams elements
  ind <- lapply(1L:n, function(i) 
    (1 + i - 1):(len_seq - n + i))
  
  if(length(d) != 1 && length(d) != n - 1)
    stop("Length of d must be 1 or n - 1")
  
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
  
  ind
}





