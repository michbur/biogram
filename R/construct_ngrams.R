#' Construct and filter n-grams
#'
#' Step-by-step builds and selects important n-grams.
#'
#' @inheritParams position_ngrams
#' @return a vector of n-grams (where \code{n} is equal to the \code{n} of the input plus one) 
#' with position information.
#' @details n-grams are built by pasting existing n-grams with unigrams extracted 
#' from them.
#' @export
#' @seealso
#' Function used by \code{build_ngrams} to extract unigrams: \code{\link{position_ngrams}}.
#' @note All n-grams must have the same length (\code{n}).
#' @examples
#' deg_seqs <- degenerate(human_cleave[, 1L:9],
#' list(`1` = c(1, 6, 8, 10, 11, 18),
#'      `2` = c(2, 13, 14, 16, 17),
#'      `3` = c(5, 19, 20),
#'      `4` = c(7, 9, 12, 15),
#'      '5' = c(3, 4)))

construct_ngrams <- function(target, seq, u, n_max) {
  unigrams <- count_ngrams(deg_seqs, n, 1L:5, pos = TRUE)
  test_unigrams <- test_features(human_cleave[, "tar"], unigrams)
  signif_ngrams <- cut(test_unigrams, breaks = c(0, 0.05, 1))[[1]]
  for(i in 1L:n_max) {
    signif_ngrams <- build_and_test(signif_ngrams, deg_seqs, i, 1L:5, 
                                    human_cleave[, "tar"])
  }
  signif_ngrams
}

build_and_test <- function(signif_ngrams, seqs, n, u, target) {
  nplusgrams <- unique(build_ngrams(signif_ngrams))
  nplusgrams_table <- do.call(rbind, strsplit(nplusgrams, "_"))
  colnames(nplusgrams_table) <- c("position", "ngram", "distance")
  unique_dists <- unique(nplusgrams_table[, "distance"])
  new_counts <- do.call(cbind, lapply(unique_dists, function(unique_dist) {
    nplusgrams_counts <- count_ngrams(deg_seqs, n, u, 
                                      d = as.numeric(strsplit(unique_dist, ".", fixed = TRUE)[[1]]),
                                      pos = TRUE)
    nplusgrams_counts[, nplusgrams[nplusgrams_table[, "distance"] == unique_dist]]
  }))
  new_test <- test_features(target, new_counts)
  cut(new_test, breaks = c(0, 0.05, 1))[[1]]
}