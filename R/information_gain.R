#' Calculate IG for single feature
#'
#' Computes information gain of single feature and target vector.
#'
#' @param feature feature vector.
#' @param target_b target in bits (as per \code{\link[bit]{as.bit}}).
#' @param len_target length of target vector.
#' @param pos_target number of positive cases in target vector.
#' @param ES numeric value of target entropy.
#' @return a single numeric value - information gain in nats.
#' @details Input looks strange, but the function was designed to be as fast 
#' as possible subroutine of \code{\link{calc_ig}} and generally should not be directly
#' called by user.
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
  
  log_crosstable <- c(entlog(crosstable[1] %/e% counts_feature[2]),
                      entlog(crosstable[3] %/e% counts_feature[2]),
                      entlog(crosstable[2] %/e% counts_feature[1]),
                      entlog(crosstable[4] %/e% counts_feature[1]))
  
  #entropy - conditional entrophy
  ES + (crosstable[1] * log_crosstable[1] +
          crosstable[3] * log_crosstable[2] +
          crosstable[2] * log_crosstable[3] + 
          crosstable[4] * log_crosstable[4])/len_target
}

#' Calculate IG of features
#'
#' Computes information gain between features and target vector.
#'
#' @inheritParams test_features
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{test_features}}.
#' @export
#' @note During calculations \eqn{0 \log 0  = 0}. For justification see References.
#' @references Cover TM, Thomas JA \emph{Elements of Information Theory, 2nd Edition}
#' Wiley, 2006.
#' @examples 
#' calc_ig(sample(0L:1, 100, replace = TRUE), matrix(sample(0L:1, 400, replace = TRUE), ncol = 4))
#' @seealso
#' Other implemented \code{\link{criterions}}.
#' 
#' Subroutine used in \code{\link[base]{apply}} statement in \code{calc_ig}: 
#' \code{\link{calc_ig_single}}

calc_ig <- function(target, features) {
  tar_bit <- as.bit(target)
  l_tar <- length(target)
  pos_tar <- sum(target)
  props_tar <- c(l_tar - pos_tar, pos_tar)/l_tar
  #entrophy
  ES <- - sum(props_tar * entlog(props_tar))
  apply(features, 2, function(single_feature) 
    calc_ig_single(single_feature, tar_bit, l_tar, pos_tar, ES))
}

#logarithm safe for entropy calculation
entlog <- function(x, ...)
  ifelse(x == 0, 0, log(x, ...))

#division safe for entropy calculation
"%/e%" <- function(x, y)
  ifelse(x == 0 && y == 0, 0, x/y)