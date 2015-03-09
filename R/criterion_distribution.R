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
#'   \item{plot_data}{A matrix with values of the criterion and their probabilities.}
#'   \item{nice_name}{'Nice' name of the criterion.}
#' }
#' @name criterion_distribution
#' @docType class
NULL

create_criterion_distribution <- function(criterion, pdf, range, unsort_criterion,
                                          unsort_prob, nice_name) {
  dist <- cbind(criterion, 
                pdf, 
                1 - rev(cumsum(rev(pdf))))
  colnames(dist) <- c("criterion", "pdf", "cdf")

  attr(dist, "plot_data") <- matrix(c(unsort_criterion, unsort_prob), ncol = 2,
                                    dimnames = list(range, c("unsort_criterion",
                                                             "unsort_prob")))
  attr(dist, "nice_name") <- nice_name
  
  dist
}