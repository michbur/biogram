#' Calculate IG for single feature
#'
#' Computes information gain of feature and target vectors.
#'
#' @param feature feature vector.
#' @param target_b target in bits (as per \code{\link[bit]{as.bit}}).
#' @param len_target length of target vector.
#' @param pos_target number of positive cases in target vector.
#' @param ES numeric value of target entropy.
#' @return a single numeric value - information gain in nats.
#' @details Input looks strange, but the function was build to be as fast 
#' as possible subroutine of \code{\link{calc_ig}}, which works on
#' many features but only one target.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' prop <- c(100 - sum(tar), sum(tar))/100
#' entr <- - sum(prop*log(prop))
#' library(bit) #used to code vector as bit
#' calc_ig_single(feat, as.bit(tar), 100, sum(tar), entr)

calc_ig_single <- function(feature, target_b, len_target, pos_target, ES) {
  crosstable <- fast_crosstable(target_b, len_target, pos_target, feature)
  counts_feature <- c(crosstable[2] + crosstable[4], crosstable[1] + crosstable[3])
  
  #entropy - conditional entrophy
  ES + (crosstable[1] * log(crosstable[1]/counts_feature[2]) +
          crosstable[3] * log(crosstable[3]/counts_feature[2]) +
          crosstable[2] * log(crosstable[2]/counts_feature[1]) + 
          crosstable[4] * log(crosstable[4]/counts_feature[1]))/len_target
}


#' Very fast 2d cross-tabulation
#'
#' Quickly cross-tabulates two binary vectors.
#'
#' @inheritParams calc_ig_single
#' @return a vector of length four: 
#' \enumerate{
#' \item target +, feature+
#' \item target +, feature-
#' \item target -, feature+
#' \item target -, feature-
#' }
#' @details Input looks strange, but the function was build to be as fast 
#' as possible subroutine of \code{\link{calc_ig}}, which works on
#' many features but only one target.
#' @note Binary vector means numeric vector with 0 or 1.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' library(bit) #used to code vector as bit
#' fast_crosstable(as.bit(tar), 100, sum(tar),  feat)

fast_crosstable <- function(target_b, len_target, pos_target, feature) {
  feature_b = as.bit(feature) #from bit library, faster than any other type
  
  #target positive and feature positive
  n_tar_f <- sum(feature_b & target_b) #simple boolean algebra to speed it more
  #feature positive
  pos_f <- sum(feature_b)
  
  c(n_tar_f, # tar +, feature +
    pos_target - n_tar_f, # tar +, feature -
    pos_f - n_tar_f, # tar -, feature +
    len_target - pos_target - pos_f + n_tar_f) # tar -, feature -
}


#' Calculate IG of features
#'
#' Computes information gain between features and target vector.
#'
#' @param target target vector.
#' @param features \code{integer} matrix of features with number of rows equal 
#' to the length of target vector.
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{test_ig}} for Monte Carlo information gain test.
#' @export
#' @examples 
#' calc_ig(sample(0L:1, 100, replace = TRUE), matrix(sample(0L:1, 400, replace = TRUE), ncol = 4))

calc_ig <- function(target, features) {
  tar_bit <- as.bit(target)
  l_tar <- length(target)
  pos_tar <- sum(target)
  props_tar <- c(l_tar - pos_tar, pos_tar)/l_tar
  #entrophy
  ES <- - sum(props_tar*log(props_tar))
  apply(features, 2, function(single_feature) 
    calc_ig_single(single_feature, tar_bit, l_tar, pos_tar, ES))
}


#' Test information gains
#'
#' Computes p-value of features' information gain.
#'
#' @inheritParams calc_ig
#' @param times number of times procedure should be repetead
#' @return a numeric vector of lenth equal to the number of features containing computed
#' information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{calc_ig}}
#' @export
#' @examples calc_ig(sample(0L:1, 100, replace = TRUE), 
#' matrix(sample(0L:1, 400, replace = TRUE), ncol = 4))

test_ig <- function(target, features, times)
  rowMeans(calc_ig(target, features) <= 
             replicate(times, calc_ig(sample(target), features)))
