#' criterion_distribution class
#'
#' A result of \code{\link{distr_crit}} function.
#'
#' @details An object of class \code{criterion_distribution} is a numeric matrix. 
#' @section Data:
#' \describe{
#'   \item{1st column:}{possible values of criterion.}
#'   \item{2nd column:}{probability density function.}
#'   \item{3rd column:}{cumulative distribution function.}
#' }
#' @section Attributes:
#' \describe{
#'   \item{XX}{fill me.}
#' }
#' @name criterion_distribution
#' @docType class
NULL

create_criterion_distribution <- function(criterion, pdf, range, unsort_criterion) {
  dist <- cbind(criterion, 
                pdf, 
                1 - rev(cumsum(rev(pdf))),
                range,
                unsort_criterion)
  colnames(dist) <- c("criterion", "pdf", "cdf", "range", "unsort_criterion")
  dist
}