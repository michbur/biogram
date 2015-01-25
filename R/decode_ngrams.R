#' Decode n-grams
#'
#' Tranforms a vector of n-grams into a human-friendly form.
#'
#' @param ngrams a \code{character} vector of n-grams.
#' @return a \code{character} vector of length equal to the number of n-grams. 
#' @export
#' @seealso
#' Assign n-grams to positions: \code{\link{position_ngrams}}.
#' 
#' Validate n-gram structure: \code{\link{is_ngram}}.
#' @examples
#' decode_ngrams(c("2_1.1.2_0.0", "3_1.1.2_0.0", "3_2.2.2_0.0"))
decode_ngrams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: paste(names(which(!validated_ngram)), collapse = ", ").")
  
  sngrams <- strsplit(ngrams, "_")
  vapply(sngrams, decode_single_ngrams, "a")
}

decode_single_ngrams <- function(splitted_ngram) {
  pos_inf <- ifelse(length(splitted_ngram) == 3, TRUE, FALSE)
  seq <- strsplit(splitted_ngram[1 + pos_inf], ".", fixed = TRUE)[[1]]
  dists <- strsplit(splitted_ngram[2 + pos_inf], ".", fixed = TRUE)[[1]]
  #distances in bar form
  bar_dists <- vapply(dists, function(i) 
    paste(rep("_", i), collapse = ""), "a")
  paste(c(vapply(1L:(length(seq) - 1), function(i)
    c(seq[i], bar_dists[i]), c("a", "a")), seq[length(seq)]), collapse = "")
}