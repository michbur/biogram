#' Create feature according to given contingency matrix
#'
#' Creates a matrix of features and target based on the values from contingency matrix.
#' 
#' @param n11 number of elements for which both target and feature equal 1.
#' @param n01 number of elements for which target and feature equal 1,0 
#' respectively.
#' @param n10 number of elements for which target and feature equal 0,1 
#' respectively.
#' @param n00 number of elements for which both target and feature equal 0.
#' @export
#' @return a matrix of 2 columns and n11+n10+n01+n00 rows. Columns represent
#' target and feature vectors, respectively.
#' @examples
#' # equivalent of 
#' #         target
#' # feature 10 375
#' #        15 600
#' target_feature <- create_feature_target(10, 375, 15, 600)
create_feature_target <- function(n11, n01, n10, n00){
  tar <- c(rep(1, n11), rep(1, n01), rep(0, n10), rep(0, n00))
  feat <- c(rep(1, n11), rep(0, n01), rep(1, n10), rep(0, n00))
  cbind(tar, feat)
}

#' Very fast 2d cross-tabulation
#'
#' Quickly cross-tabulates two binary vectors.
#'
#' @inheritParams calc_ig
#' @return a vector of length four: 
#' \enumerate{
#' \item target +, feature+
#' \item target +, feature-
#' \item target -, feature+
#' \item target -, feature-
#' }
#' @details Input looks odd, but the function was build to be as fast 
#' as possible subroutine of \code{\link{calc_ig}}, which works on
#' many features but only one target.
#' @note Binary vector means a numeric vector with 0 or 1.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' library(bit) # used to code vector as bit
#' fast_crosstable(as.bit(tar), length(tar), sum(tar),  feat)

fast_crosstable <- function(target_b, len_target, pos_target, feature) {
  # no input tests - every if clause slows a little bit
  feature_b = as.bit(feature) # from bit library, faster than any other type
  
  # target positive and feature positive
  n_tar_f <- sum(feature_b & target_b) # simple boolean algebra to speed it more
  # feature positive
  pos_f <- sum(feature_b)
  
  c(n_tar_f, # tar +, feature +
    pos_target - n_tar_f, # tar +, feature -
    pos_f - n_tar_f, # tar -, feature +
    len_target - pos_target - pos_f + n_tar_f) # tar -, feature -
}

#' Convert sequence list to matrix
#'
#' Converts a list of sequences to a matrix.
#'
#' @param seq_list list of sequences (e.g. as returned by 
#' the \code{\link[seqinr]{read.fasta}} function).
#' @return a matrix of sequences.
#' @details a matrix with the number of rows equal to the number of 
#' sequences and the number of columns equal to the length of the 
#' longest sequence. If sequences do not have the same length, 
#' shorter sequences are filled with NAs.
#' @export
#' @examples 
#' # list of sequences
#' sample_list <- lapply(round(runif(5, 6, 15), 0), function(i)
#'   sample(c("A", "C", "G", "T"), i, replace = TRUE))
#' # convert to matrix
#' list2matrix(sample_list)


list2matrix <- function(seq_list) {
  max_len <- max(lengths(seq_list))
  t(sapply(seq_list, function(single_seq)
    c(single_seq, rep(NA, max_len - length(single_seq)))
  ))
}