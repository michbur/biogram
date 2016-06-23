#' Construct and filter n-grams
#'
#' Builds and selects important n-grams stepwise. 
#'
#' @param target \code{integer} vector with target information (e.g. class labels).
#' @inheritParams count_ngrams
#' @param n_max size of constructed n-grams.
#' @param conf_level confidence level.
#' @param gap \code{logical}, if \code{TRUE} gaps are used. See Details.
#' @param use_heuristics, if \code{FALSE} then all n-grams are tested. This may
#' slow down computations significantly 
#' @return a vector of n-grams.
#' @details 
#' 
#' \code{construct_ngrams} starts by 
#' extracting unigrams from the sequences, pasting them together in all combination and 
#' choosing from them significant features (with p-value below \code{conf_level}). The 
#' chosen n-grams are further extended to the specified by \code{n_max} size by pasting 
#' unigrams at both ends.
#' 
#' The \code{gap} parameter determines if \code{construct_ngrams} performs the
#' feature selection on exact n-grams (\code{gap} equal to FALSE) or on all features in the 
#' Hamming distance 1 from the n-gram (\code{gap} equal to TRUE).
#' 
#' @export
#' @seealso
#' Feature filtering method: \code{\link{test_features}}.
#' @examples
#' # to make the example faster, we run construct_ngrams() on the 
#' # subset of data
#' deg_seqs <- degenerate(human_cleave[c(1L:100, 801L:900), 1L:9],
#' list(`1` = c(1, 6, 8, 10, 11, 18),
#'      `2` = c(2, 13, 14, 16, 17),
#'      `3` = c(5, 19, 20),
#'      `4` = c(7, 9, 12, 15),
#'      '5' = c(3, 4)))
#' bigrams <- construct_ngrams(human_cleave[c(1L:100, 801L:900), "tar"], deg_seqs, 1L:5, 2)

construct_ngrams <- function(target, seq, u, n_max, conf_level = 0.95, gap = TRUE,
                             use_heuristics = TRUE) {
  # build unigrams
  unigrams <- count_ngrams(seq, 1, u, pos = TRUE)
  # filter unigrams
  test_unigrams <- test_features(target, unigrams)
  signif_ngrams <- cut(test_unigrams, breaks = c(0, 1 - conf_level, 1))[[1]]
  
  res <- list()
  
  if(use_heuristics){
    res[[1]] <- signif_ngrams
  } else {
    res[[1]] <- signif_ngrams
    signif_ngrams <- paste0(create_ngrams(1, 1L:4, possible_grams = ncol(seq)), "_0")
  }
  
  
  for(i in 2L:n_max) {
    if(length(signif_ngrams) != 0) {
      #code from build_and_test because we need nplusgrams for non-heuristic version
      seq_length <- ncol(seq)
      nplusgrams <- unique(unlist(lapply(signif_ngrams, function(single_ngram)
        add_1grams(single_ngram, u, seq_length))))
      new_counts <- count_specified(seq, nplusgrams)
      new_test <- test_features(target, new_counts)
      
      if(use_heuristics){
        signif_ngrams <- cut(new_test, breaks = c(0, conf_level, 1))[[1]]
        res[[i]] <- signif_ngrams
      } else {
        res[[i]] <- cut(new_test, breaks = c(0, conf_level, 1))[[1]]
        signif_ngrams <- nplusgrams
      }
    } else {
      NULL
    }
  }
  res
}

build_and_test <- function(signif_ngrams, seq, u, target, conf_level) {
  seq_length <- ncol(seq)
  nplusgrams <- unique(unlist(lapply(signif_ngrams, function(single_ngram)
    add_1grams(single_ngram, u, seq_length))))
  new_counts <- count_specified(seq, nplusgrams)
  new_test <- test_features(target, new_counts)
  cut(new_test, breaks = c(0, conf_level, 1))[[1]]
}
