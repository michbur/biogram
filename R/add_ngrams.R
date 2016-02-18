#' Add 1-grams
#'
#' Builds (n+1)-grams from n-grams.
#'
#' @param ngram a single n-gram.
#' @param u \code{integer}, \code{numeric} or \code{character} vector of all possible 
#' unigrams.
#' @return vector of n-grams (where \code{n} is equal to the \code{n} of the input plus one) 
#' with position information.
#' @details n-grams are built by pasting every possible unigram in the every possible free 
#' position. 
#' @export
#' @examples
#' add_1grams("1_2.3.4_3.0", 1L:5)

add_1grams <- function(ngram, u) {
  validated_ngram <- sapply(ngram, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: ", paste(names(which(!validated_ngram)), 
                                     collapse = ", "))
  
  # check if there is information about position
  sngrams <- strsplit(ngram, "_")
  pos_inf <- ifelse(length(sngrams[[1]]) == 3, TRUE, FALSE)
  
  
  splitted_ngrams <- c("_", strsplit(decode_ngrams(ngram), "")[[1]], "_")
  
  if(pos_inf) {
    positioned_ngrams <- position_ngrams(ngram, df = TRUE, 
                                         unigrams_output = FALSE)
    
    #remove adding unigrams in position 0
    if(positioned_ngrams[["position"]] == 1) {
      splitted_ngrams <- splitted_ngrams[-1]
    }
  }
  
  
  res <- lapply(which(splitted_ngrams == "_"), function(single_position) {
    sapply(u, function(single_u) {
      splitted_ngrams[single_position] <- single_u
      code_ngrams(paste0(splitted_ngrams, collapse = ""))
    })
  })
  
  if(pos_inf) {
    res <- lapply(res[-1], function(single_res) {
      paste0(positioned_ngrams[["position"]], "_", single_res)
    })
    
    if(positioned_ngrams[["position"]] == 1) {
      res[[1]] <- paste0(positioned_ngrams[["position"]], "_", res[[1]])
    } else {
      res[[1]] <- paste0(positioned_ngrams[["position"]] - 1, "_", res[[1]])
    }
    
  }
  
  unlist(res)
}

