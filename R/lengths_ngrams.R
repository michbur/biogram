#' Get lengths of the n-grams
#'
#' Computes the length of n-grams.
#'
#' @inheritParams decode_ngrams
#' @return A \code{numeric} vector of n-gram lengths.
#' @export
#' @examples 
#' lengths_ngrams(c("2_1.1.2_0.1", "3_1.1.2_2.0", "3_2.2.2_0.0"))

lengths_ngrams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: paste(names(which(!validated_ngram)), collapse = ", ").")
  
  nchar(decode_ngrams(ngrams))
}
