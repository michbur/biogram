#' Permutation test for feature selection
#'
#' Performs a feature selection on positioned N-gram data using a Fisher's 
#' permutation test.
#'
#' @inheritParams calc_ig
#' @param times number of times procedure should be repetead
#' @param the criterion used in permutation test.
#' @details Currently implemented criterions:
#' \itemize{
#' \item{"\code{ig}" - information gain}
#' }
#' @return a numeric vector of lenth equal to the number of features containing computed
#' information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{calc_ig}}
#' @export
#' @examples calc_ig(sample(0L:1, 100, replace = TRUE), 
#' matrix(sample(0L:1, 400, replace = TRUE), ncol = 4))

test_features <- function(target, features, times, criterion = "ig") {
  if(criterion == "ig")
    chosen_test <- calc_ig
  rowMeans(calc_ig(target, features) <= 
             replicate(times, chosen_test(sample(target), features)))
}