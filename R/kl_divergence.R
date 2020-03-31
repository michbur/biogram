#' Calculate KL divergence of features
#'
#' Computes Kullback-Leibler divergence between features and target vector.
#'
#' @inheritParams calc_ig
#' @return A \code{numeric} vector of length 1 representing Kullback-Leibler divergence 
#' value.
#' @seealso \code{\link{test_features}}.
#' Kullback-Leibler divergence is calculated using \code{\link[entropy]{KL.plugin}}.
#' 
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' 
#' The function was designed to be as fast as possible subroutine of 
#' \code{\link{calc_criterion}} and might be cumbersome if directly called by a user.
#' @export
#' @references Kullback S, Leibler RA \emph{On information and sufficiency}. Annals
#' of Mathematical Statistics 22 (1):79-86, 1951.
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' calc_kl(feat, tar, 100, sum(tar))
calc_kl <- function(feature, target, len_target, pos_target) {
  crosstable <- fast_crosstable(target, len_target, pos_target, feature)
  counts_feature <- c(crosstable[2] + crosstable[4], crosstable[1] + crosstable[3])
  
  KL.plugin(crosstable[c(2, 4)]/counts_feature[1], crosstable[c(1, 3)]/counts_feature[2])
}