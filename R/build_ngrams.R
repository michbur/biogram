#' Build n-grams
#'
#' Builds longer n-grams from shorter n-grams.
#'
#' @inheritParams position_ngrams
#' @return a vector of n-grams (where \code{n} is equal to the \code{n} of the input plus one) 
#' with position information.
#' @details n-grams are build iteratively by pasting existing n-grams with unigrams extracted 
#' from them.
#' @export
#' @seealso
#' Function used by \code{build_ngrams} to extract unigrams: \code{\link{position_ngrams}}.
#' @note All n-grams must have the same length (\code{n}).
#' @examples
#' build_ngrams(c("1_1_0", "2_1_0", "5_1_0", "7_1_0", "4_2_0", 
#' "5_2_0", "7_2_0", "8_5_0"))

build_ngrams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: paste(names(which(!validated_ngram)), collapse = ", ").")
  
  splitted_ngrams <- strsplit(ngrams, "_")
  if(unique(sapply(splitted_ngrams, length)) != 3)
    stop("Use only n-grams with position information.")
  n <- sapply(strsplit(sapply(splitted_ngrams, function(ngram) 
    ngram[2]), ".", fixed = TRUE), length)
  if(length(unique(n)) != 1)
    stop("Unequal n-gram size. Use n-grams with the same size (n).")
  if(unique(n) == 1)
    build_bigrams(ngrams)
}

#build unigrams
build_bigrams <- function(ngrams) {
  positioned_ngrams <- position_ngrams(ngrams, df = TRUE)
  positions <- unique(positioned_ngrams[["position"]])
  #remove last position, because there is nothing to paste it with
  positions <- positions[-which.max(positions)]
  
  res <- unlist(lapply(positions, function(single_position) {
    chosen_ugrams <- positioned_ngrams[positioned_ngrams[["position"]] == single_position, 
                                        "ngram"]
    other_ugrams <- positioned_ngrams[positioned_ngrams[["position"]] > single_position, ]
    #position in other_ugrams is now distance between single_position and their position
    other_ugrams[["position"]] <- other_ugrams[["position"]] - single_position - 1
    
    lapply(chosen_ugrams, function(single_ugram)
      apply(other_ugrams, 1, function(other_ugram)
        paste0(single_position, "_", #position 
               single_ugram, ".", other_ugram[1], #ngram 
               "_", other_ugram[2]))) #distance
  }))
  names(res) <- NULL
  res
}