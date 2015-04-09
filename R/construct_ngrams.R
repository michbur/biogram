#' Construct and filter n-grams
#'
#' Step-by-step builds and selects important n-grams.
#'
#' @inheritParams test_features
#' @inheritParams count_ngrams
#' @param n_max size of constructed n-grams.
#' @return a vector of n-grams.
#' @details \code{construct_ngrams} extracts unigrams from the sequences and filters 
#' significant features. From filtered unigrams, the function builds all possible 
#' bigrams using \code{\link{add_1grams}} and so on, till the required size of 
#' the n-gram is achieved. On each step, obtained n-grams are filtered using 
#' \code{\link{test_features}}.
#' @export
#' @examples
#' \donttest{
#' deg_seqs <- degenerate(human_cleave[, 1L:9],
#' list(`1` = c(1, 6, 8, 10, 11, 18),
#'      `2` = c(2, 13, 14, 16, 17),
#'      `3` = c(5, 19, 20),
#'      `4` = c(7, 9, 12, 15),
#'      '5' = c(3, 4)))
#' bigrams <- construct_ngrams(human_cleave[, "tar"], deg_seqs, 1L:5, 2)
#' }

construct_ngrams <- function(target, seq, u, n_max) {
  #build unigrams
  unigrams <- count_ngrams(seq, 1, u, pos = TRUE)
  #filter unigrams
  test_unigrams <- test_features(target, unigrams)
  signif_ngrams <- cut(test_unigrams, breaks = c(0, 0.05, 1))[[1]]
  for(i in 2L:n_max) {
    signif_ngrams <- build_and_test(signif_ngrams, seq, i, 1L:5, 
                                    target)
  }
  signif_ngrams
}

build_and_test <- function(signif_ngrams, seqs, n, u, target) {
  nplusgrams <- unique(add_1grams(signif_ngrams))
  new_counts <- count_specified(seqs, nplusgrams)
  new_test <- test_features(target, new_counts)
  cut(new_test, breaks = c(0, 0.05, 1))[[1]]
}