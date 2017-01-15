#' regional_param class
#'
#' List of rules defining the region.
#'
#' @details An object of the \code{regional_param} class is a list consisting of all rules 
#' necessary to properly build a region.
#' @section Attributes:
#' \describe{
#'   \item{reg_len}{the number of unigrams inside the region. Might be 0}
#'   \item{prop_ranges}{required intervals of properties of unigrams in the region}
#'   \item{exactness}{a \code{numeric} value between 0 and 1 defining how stricly 
#' unigrams are kept within \code{prop_ranges}. If 1, only unigrams within 
#' \code{prop_ranges} are inside the region. if 0.9, there is 10% chance that 
#' unigrams that are not in the \code{prop_ranges} will be inside the region.}
#' }
#' @name regional_param
#' @docType class
#' @seealso 
#' \code{\link{generate_sequence}}
NULL

create_region <- function(reg_len, prop_ranges, exactness) {
  if(reg_len < 0) {
    stop("reg_len cannot be smaller than 0.")
  }
  
  if(exactness < 0 | exactness > 1) {
    stop("exactness must be between 0 and 1.")
  }
  
  res <- list(reg_len = reg_len,
              prop_ranges = prop_ranges,
              exactness = exactness)
  class(res) <- "regional_param"
  res
}