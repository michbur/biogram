#' Build n-gram table
#'
#' Builds table containing counts n-grams from numeric sequence table.
#'
#' @param dat sequence table.
#' @param dists list of distances, see Examples
#' @inheritParams count_ngrams
#' @return a dataframe of all n-grams counts.
#' @details Input data must be a matrix or data frame of numeric elements.
#' @export
#' @examples
#' random_seqs <- matrix(sample(1L:4, 800, replace = TRUE), nrow = 40)
#' #distance for 1-grams, 2-grams and 3-grams
#' dist_list <- list(matrix(0L:7), expand.grid(0L:3, 0L:3))
#' gramize_data(random_seqs, dist_list)
#' 

gramize_data <- function(dat, dists, pos = FALSE, scale = FALSE) {
  res <- lapply(dists, function(i)
    mdcount_ngrams(dat, length(length(i[,1]) + 1), i, pos = pos, scale = scale))
  do.call(cbind, res)
}

#helper function, calculate n-grams for multiple distances
mdcount_ngrams <- function(seqs, n_gram, dists, pos, scale) {
  #possible number of n-grams
  possib_grams <- NULL
  if(pos)
    possib_grams <- ncol(seqs) - n_gram
  cur_features <- create_ngrams(n_gram, 1:4, pos)
  do.call(cbind, lapply(1:nrow(dists), function(cur_dist) {
    gram_table <- apply(seqs, 1, function(single_seq) 
      count_ngrams(seq = single_seq, 
                   feature_list = cur_features, 
                   n = n_gram, 
                   d = dists[cur_dist, ], 
                   scale = scale))
    rownames(gram_table) <- paste0(rownames(gram_table), "_", paste(dists[cur_dist, ], 
                                                                    collapse = "_"))
    t(gram_table)
  }))
}


#' Get all possible n-grams
#'
#' Creates vector of all posible n_grams.
#'
#' @param n size of gram (i.e. 2-grams: "11", "12", ...).
#' @param u unigrams.
#' @param possible_grams number of possible n-grams. If not \code{NULL} n-grams do not
#' contain information about position
#' @return a character vector. Elements of n-gram are separated by dot.
#' @details Input data must be a matrix or data frame of numeric elements.
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
  if (!is.null(NULL))
    res <- as.vector(sapply(res, function(i) paste(1L:possible_grams, i, sep = "_")))
  res
}

#' Detect and count n-grams in single sequence
#'
#' Counts all n-grams present in given sequence.
#'
#' @param seq \code{integer} vector describing sequence.
#' @param feature_list list of features (see Note).
#' @param n size of n-grams.
#' @param d distance between elements of n-gram (0 means consecutive elements). See
#' Details.
#' @param pos \code{logical}, if \code{TRUE} n_grams contains position information.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @return a named \code{integer} vector. Elements of n-gram are separated by dot.
#' @note List of possible n-grams must be calculated outside of the function.
#' @details The length of \code{distance} vector should be always n - 1. For example 
#' when \code{n} = 3, \code{d} = c(1, 2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2, 0, 1) means A__AA_A.
#' @export
#' @seealso For convenient wrapper for multidistance calculations 
#' see \code{\link{gramize_data}}.
#' @examples 
#' #trigrams for nucleotides
#' count_ngrams(sample(1L:4, 30, replace = TRUE), create_ngrams(3, 1L:4), 3)

count_ngrams <- function(seq, feature_list, n, d = 0, pos = FALSE, scale = FALSE) {
  #feature list(list of possible n-grams) is outside, because count_ngrams is meant to
  #be inside the loop
  if (n > 1) {
    grams <- count_ngrams_helper(seq, n, d)
  } else {
    grams <- seq
  }
  
  if (pos) {
    paste(1L:length(grams), grams, sep = "_")
  } 
  
  res <- sapply(feature_list, function(i)
      sum(i == grams))
  
  if (scale)
    res <- res/(length(seq) - n - sum(d) + 1)
  res
}

#helper function for counting n-grams when n > 1
count_ngrams_helper <- function(seq, n, d) {
  #creates grams from given function
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive letter (a vector of distances)
  ind <- lapply(1:n, function(i) 
    (1 + i - 1):(length(seq) - n + i))
  
  if(length(d) != 1 && length(d) != n - 1)
    stop("Length of d must be 1 or n - 1")
  if(length(d) == 1 && n > 2)
    d <- rep(d, n - 1)
  
  if(sum(d) > 0) {
    ind[-1] <- lapply(1L:length(d), function(i)
      ind[[i + 1]] + sum(d[1:i]))
    not_taken <- ind[[1]][(length(ind[[1]]) - sum(d) + 1):length(ind[[1]])]
    ind <- lapply(ind, function(i) i[-not_taken])
  }
  
  pair_matrix <- do.call(cbind, lapply(ind, function(i) seq[i]))
  apply(pair_matrix, 1, function(x) 
    paste(x, collapse="."))
}




