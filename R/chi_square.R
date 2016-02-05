#' Calculate Chi-squared-based measure
#'
#' Computes Chi-squared-based measure between features and target vector.
#'
#' @inheritParams calc_ig
#' @return A \code{numeric} vector of length 1 representing computed Chi-square values.
#' @seealso \code{\link{test_features}}.
#' 
#' \code{\link[stats]{chisq.test}} - Pearson's chi-squared test for count data.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' 
#' The function was designed to be as fast as possible subroutine of 
#' \code{\link{calc_criterion}} and might be cumbersome if directly called by a user.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' library(bit) # used to code vector as bit
#' calc_cs(feat, as.bit(tar), 100, sum(tar))
calc_cs <- function(feature, target_b, len_target, pos_target) {
  crosstable_m <- matrix(fast_crosstable(target_b, len_target, pos_target, feature), nrow = 2, byrow = TRUE)

  # copied from original chisq.test in stats package
  sr <- rowSums(crosstable_m)
  sc <- colSums(crosstable_m)
  n <- sum(crosstable_m)
  E <- outer(sr, sc, "*")/n
  sum((crosstable_m - E)^2/E)
}