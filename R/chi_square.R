#' Calculate Chi-squared-based measure
#'
#' Computes Chi-squared-based measure between features and target vector.
#'
#' @inheritParams calc_ig
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed Chi-square values.
#' @seealso \code{\link{test_features}}.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' 
#' The function was designed to be as fast as possible subroutine of 
#' \code{\link{calc_criterion}} and might be cumbersome if directly called by a user.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' prop <- c(100 - sum(tar), sum(tar))/100
#' entr <- - sum(prop*log(prop))
#' library(bit) #used to code vector as bit
#' calc_chisq(feat, as.bit(tar), 100, sum(tar), entr)
calc_chisq <- function(feature, target_b, len_target, pos_target, ES) {
  crosstable <- fast_crosstable(target_b, len_target, pos_target, feature)
  counts_feature <- c(crosstable[2] + crosstable[4], crosstable[1] + crosstable[3])
  
  browser()
}