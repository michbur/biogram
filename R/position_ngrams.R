#' Position n-grams
#'
#' Tranforms a vector of positioned n-grams into a list of positions filled with unigrams 
#' appearing on them.
#'
#' @param ngrams a vector of positioned n-grams (as created by \code{\link{count_ngrams}}).
#' @return a list of length equal to the number of unique positions present in n-grams. Each 
#' element of the list contains unigrams that are present on this position.
#' @export
#' @seealso
#' Change n-gram name to human-friendly form: \code{\link{decode_ngrams}}.
#' 
#' Validate n-gram structure: \code{\link{is_ngram}}.
#' @examples
#' position_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0"))

position_ngrams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: paste(names(which(!validated_ngram)), collapse = ", ").")
  
  sngrams <- strsplit(ngrams, "_")
  #check if there is information about position
  if(length(sngrams[[1]]) != 3)
    stop("n-grams do not have position information.")
  
  #table of positions
  pos_table <- do.call(rbind, lapply(sngrams, function(single_ngram) {
    unigrams <- strsplit(single_ngram[[2]], ".", fixed = TRUE)[[1]]
    dists <- strsplit(single_ngram[[3]], ".", fixed = TRUE)[[1]]
    #positions of unigrams
    uni_positions <- as.numeric(single_ngram[[1]])
    for (next_unigram in as.numeric(dists))
      uni_positions[length(uni_positions) + 1] <- next_unigram + uni_positions[length(uni_positions)] + 1
    data.frame(unigrams = as.numeric(unigrams), pos = uni_positions)
  }))
  
  res <- lapply(sort(unique(pos_table[["pos"]])), function(unique_pos) {
    sort(pos_table[pos_table[["pos"]] == unique_pos, "unigrams"])
  })
  
  names(res) <- sort(unique(pos_table[["pos"]]))
  res
}

