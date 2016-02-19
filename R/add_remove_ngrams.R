#' Add 1-grams
#'
#' Builds (n+1)-grams from n-grams.
#'
#' @param ngram a single n-gram.
#' @param u \code{integer}, \code{numeric} or \code{character} vector of all possible 
#' unigrams.
#' @param seq_length length of an origin sequence.
#' @return vector of n-grams (where \code{n} is equal to the \code{n} of the input plus one). 
#' @details n-grams are built by pasting every possible unigram in the every possible free 
#' position. The total length of n-gram (n plus total distance between elements of the n-gram) 
#' is limited by the length of an origin sequence, because the n-gram cannot be longer than 
#' an origin sequence.
#' @seealso Reverse function: \code{\link{gap_ngrams}}.
#' @export
#' @examples
#' add_1grams("1_2.3.4_3.0", 1L:4, 8)
#' 
#' add_1grams("a.a_1", c("a", "b", "c"), 4)

add_1grams <- function(ngram, u, seq_length) {
  
  # no need to validate n-grams, decode does it for us
  decoded <- strsplit(decode_ngrams(ngram), "")[[1]]
  
  # check if there is information about position
  sngrams <- strsplit(ngram, "_")
  pos_inf <- ifelse(length(sngrams[[1]]) == 3, TRUE, FALSE)
  
  if(pos_inf) {
    start_position <- position_ngrams(ngram, df = TRUE, unigrams_output = FALSE)[["position"]]
    add_1grams_positioned(decoded, start_position, seq_length, u, add_position = TRUE)
  } else {
    possible_pos <- 1L:(seq_length - length(decoded) + 1)
    # added unique, because without position information some n-gram might be redundant
    unique(unlist(lapply(possible_pos, function(single_pos)
      add_1grams_positioned(decoded, start_position = single_pos, seq_length, u, add_position = FALSE)
    )))
  }
}

add_1grams_positioned <- function(decoded, start_position, seq_length, u, add_position) {
  splitted_ngrams <- c(rep("_", start_position - 1), 
                       decoded,
                       rep("_", seq_length - length(decoded) - start_position + 1))
  
  ngrams <- lapply(which(splitted_ngrams == "_"), function(single_position) {
    sapply(u, function(single_u) {
      splitted_ngrams[single_position] <- single_u
      code_ngrams(paste0(splitted_ngrams, collapse = ""))
    })
  })
  
  if(add_position) {
    # add position information
    positions <- 1L:seq_length
    positions[positions > start_position] <- start_position
    positions <- positions[splitted_ngrams == "_"]
    
    ngrams <- lapply(1L:length(ngrams), function(position_id) 
      paste0(positions[position_id], "_", ngrams[[position_id]])
    )
  }
  
  unlist(ngrams)
}


#' Gap n-grams
#'
#' Introduces gaps in the n-grams.
#'
#' @inheritParams position_ngrams
#' @return A \code{character} vector of (n-1)-grams with introduced gaps.
#' @details A single element of the input n-gram at a time will be replaced 
#' by a gap. For example, introducing gaps in n-gram \code{2_1.1.2_0.1} 
#' will results in three n-grams: \code{3_1.2_1} (where the \code{2_1_0} unigram 
#' was replaced by a gap), \code{2_1.2_2} and \code{2_1.1_0}.
#' @seealso Reverse function: \code{\link{add_1grams}}.
#' @export
#' @examples 
#' gap_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0"))
#' gap_ngrams(c("1.1.2_0.1", "1.1.2_0.0", "2.2.2_0.0"))

gap_ngrams <- function(ngrams) {
  # check if unigrams are there
  if(any(nchar(ngrams2df(ngrams)[["ngram"]]) == 1))
    stop("'ngrams' have n bigger than 1.")
  
  unique(unlist(lapply(ngrams, gap_single_ngram)))
}


gap_single_ngram <- function(ngram) {
  # check if there is information about position
  sngrams <- strsplit(ngram, "_")
  pos_inf <- ifelse(length(sngrams[[1]]) == 3, TRUE, FALSE)
  
  # no need to validate n-grams, decode does it for us
  decoded <- strsplit(decode_ngrams(ngram), "")[[1]]
  
  res <- sapply(which(decoded != "_"), function(single_position) {
    decoded[single_position] <- "_"
    code_ngrams(paste0(decoded, collapse = ""))
  })
  
  if(pos_inf) {
    start_position <- position_ngrams(ngram, df = TRUE, unigrams_output = FALSE)[["position"]]
    
    # if the first element of n-gram is removed, the position should be ajudsted
    pos_adj <- which(decoded[-1] != "_")[1]
    
    positions <- rep(start_position, length(res))
    positions[1] <- positions[1] + pos_adj
    
    res <- vapply(1L:length(res), function(res_id)
      paste0(positions[res_id], "_", res[res_id]), "a")
  }
  
  res
}