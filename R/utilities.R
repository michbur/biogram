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

#' 2d cross-tabulation
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
#' @details Input looks odd, but the function was build to be fast
#' subroutine of \code{\link{calc_ig}}, which works on
#' many features but only one target.
#' @note Binary vector means a numeric vector with 0 or 1.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' fast_crosstable(tar, length(tar), sum(tar),  feat)

fast_crosstable <- function(target, len_target, pos_target, feature) {
  # target positive and feature positive
  n_tar_f <- sum(feature == 1 & target == 1) 
  # feature positive
  pos_f <- sum(feature)
  
  c(n_tar_f, # tar +, feature +
    pos_target - n_tar_f, # tar +, feature -
    pos_f - n_tar_f, # tar -, feature +
    len_target - pos_target - pos_f + n_tar_f) # tar -, feature -
}


#' Binarize
#'
#' Binarizes a matrix.
#' 
#' @param x \code{matrix} or \code{\link[slam]{simple_triplet_matrix}}.
#' @export
#' @return a \code{matrix} or \code{\link[slam]{simple_triplet_matrix}} 
#' (depending on the input).

binarize <- function(x) {
  nonbin_matrix <- if(is.simple_triplet_matrix(x)) {
    as.matrix(x)
  } else {
    x
  }
    
  nonbin_matrix <- nonbin_matrix > 0
  storage.mode(nonbin_matrix) <- "integer"
  
  if(is.simple_triplet_matrix(x)) {
    as.simple_triplet_matrix(nonbin_matrix)
  } else {
    nonbin_matrix
  }
}
