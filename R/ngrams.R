#' Get all possible n-grams
#'
#' Creates vector of all posible n_grams (for given \code{n}).
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





