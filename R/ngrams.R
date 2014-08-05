#' Build n-gram table
#'
#' Builds table containing counts of 1-, 2- and 3- grams from numeric sequence table.
#'
#' @param dat sequence table.
#' @param d2 distances for 2-grams.
#' @param d3 distances for 3-grams.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @return a dataframe of all n-grams counts.
#' @details Input data must be a matrix or data frame of numeric elements.
#' @export

gramize_data <- function(dat, d2, d3, scale) {
  a_grams1 <- nuc_grams(dat, 1, matrix(0), scale = scale)
  a_grams2 <- nuc_grams(dat, 2, d2, scale = scale)
  a_grams3 <- nuc_grams(dat, 3, d3, scale = scale)
  do.call(cbind, list(a_grams1, a_grams2, a_grams3))
}

#' Get all possible n-grams
#'
#' Creates vector of all posible n_grams.
#'
#' @param n size of gram (i.e. 2-grams: "11", "12", ...).
#' @param u unigrams.
#' @return a character vector. Elements of n-gram are separated by dot.
#' @details Input data must be a matrix or data frame of numeric elements.
#' @export
#' @examples 
#' #bigrams for nucleotides
#' create_ngrams(2, 1L:4)

create_ngrams <- function(n, u) {
  grid_list <- lapply(1L:n, function(i) u)
  apply(expand.grid(grid_list), 1, function(x)
    paste(x, collapse = "."))
}




#' Count existing n-grams
#'
#' Counts all n-grams present in given sequence.
#'
#' @param seq \code{integer} vector describing sequence.
#' @param feature_list list of features (see Details).
#' @param n size of n-grams.
#' @param d distance between elements of n-gram (0 means consecutive elements).
#' @param scale \code{logical}, if \code{TRUE} output data is normalized.
#' @return a named \code{integer} vector. Elements of n-gram are separated by dot.
#' @details List of possible n-grams must be calculated outside of the function
#' @export
#' @examples 
#' #trigrams for nucleotides
#' count_ngrams(sample(1L:4, 30, replace = TRUE), create_ngrams(3, 1L:4), 3)

count_ngrams <- function(seq, feature_list, n, d = 0, scale = FALSE) {
  #feature list(list of possible n-grams) is outside, because count_ngrams is meant to
  #be inside the loop
  if (n > 1) {
    grams <- create_grams(seq, n, d)
  } else {
    grams <- seq
  }
  res <- sapply(feature_list, function(i)
    sum(i == grams))
  if (scale)
    res <- res/(length(seq) - n - sum(d) + 1)
  res
}


create_grams <- function(seq, n, d) {
  #creates grams from given function
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive letter (a vector of distances)
  #for example d = c(1, 2) means A_A__A and d = c(2, 0, 1) means A__AA_A
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







nuc_grams <- function(seqs, n_gram, dists, scale) {
  cur_features <- create_ngrams(n_gram, 1:4)
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

