#' Decode n-grams
#'
#' Tranforms a vector of positioned n-grams into a list of positions filled with unigrams 
#' appearing on them.
#'
#' @param ngrams a vector of n-grams.
#' @return a list of length equal to the number of positions present in n-grams. Each element of 
#' the list contains unigrams that are present on this position.
#' @export

decode_ngrams <- function(ngrams) {
  
  sngrams <- strsplit(ngrams, "_")
  #check if there is information about position
  if(length(sngrams[[1]]) != 3)
    stop("Only n-grams with position can be decoded.")
  
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