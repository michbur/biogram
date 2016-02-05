#' Count specified n-grams
#'
#' Counts specified n-grams in the input sequence(s).
#'
#' @param seq vector or matrix describing sequence(s). 
#' @param ngrams vector of n-grams.
#' @return A \code{\link[slam]{simple_triplet_matrix}} where columns represent
#' n-grams and rows sequences.
#' @export
#' @details \code{\link{count_specified}} counts only selected n-grams declared by
#' user in the \code{ngrams} parameter. Declared n-grams must be written using the
#' \code{biogram} notation.
#' @seealso Count all possible n-grams: \code{\link{count_ngrams}}.
#' @examples
#' seqs <- matrix(c(1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 3, 4, 1, 2, 2, 4), nrow = 2)
#' count_specified(seqs, ngrams = c("1.1.1_0.0", "2.2.2_0.0", "1.1.2_0.0"))
#' 
#' seqs <- matrix(sample(1L:5, 200, replace = TRUE), nrow = 20)
#' count_specified(seqs, ngrams = c("2_4.2_0", "2_1.4_0", "3_1.3_0",
#'                                  "2_4.2_1", "2_1.4_1", "3_1.3_1",
#'                                  "2_4.2_2", "2_1.4_2", "3_1.3_2"))

count_specified <- function(seq, ngrams) {
  # validate n-grams
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: ", paste(names(which(!validated_ngram)), collapse = ", "))
  
  # if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  # length of sequence
  len_seq <- ncol(seq)
  # number of sequences
  n_seqs <- nrow(seq)
  
  df <- ngrams2df(ngrams)
  
  # splitted ngrams
  sn_grams <- strsplit(df[, "ngram"], ".", fixed = TRUE)
  
  ngrams_list <- lapply(1L:nrow(df), function(single_ngram) 
    list(ngram = strsplit(df[single_ngram, "ngram"], ".", fixed = TRUE)[[1]],
         distance = as.numeric(strsplit(df[single_ngram, "distance"], ".", fixed = TRUE)[[1]])
    )
  )
  
  # when ncol(df) == 3 n-grams are positioned
  if(ncol(df) == 3)
    ngrams_list <- lapply(1L:nrow(df), function(single_ngram) 
      c(ngrams_list[[single_ngram]], position = as.numeric(df[single_ngram, "position"]))
    )
  
  names(ngrams_list) <- ngrams

  # when ncol(df) == 3 n-grams are positioned
  res <- if(ncol(df) == 3) {
    vapply(ngrams_list, function(single_ngram) 
      count_single_positioned_ngram(single_ngram, seq, len_seq), rep(0, nrow(seq)))
  } else {
    vapply(ngrams_list, function(single_ngram) 
      count_single_unpositioned_ngram(single_ngram, seq, len_seq), rep(0, nrow(seq)))
  }
  
  if(class(res) == "numeric") {
    res <- matrix(res, ncol = 1)
  }
  # name columns
  colnames(res) <- ngrams
  
  as.simple_triplet_matrix(res)
}


count_single_positioned_ngram <- function(single_ngram, seq, len_seq) {
  all_ngram_pos <- do.call(rbind, get_ngrams_ind(len_seq, length(single_ngram[["ngram"]]), single_ngram[["distance"]]))
  as.numeric(rowSums(seq[, all_ngram_pos[, single_ngram[["position"]]]] == single_ngram[["ngram"]]) == 2)
}

count_single_unpositioned_ngram <- function(single_ngram, seq, len_seq) {
  all_ngram_pos <- do.call(rbind, get_ngrams_ind(len_seq, length(single_ngram[["ngram"]]), single_ngram[["distance"]]))
  
  vapply(1L:nrow(seq), function(single_seq) {
    sum(apply(all_ngram_pos, 2, function(single_ngram_pos)
      as.numeric(all(as.character(seq[single_seq, single_ngram_pos]) == single_ngram[["ngram"]]))), na.rm = TRUE)
  }, 0)
}

