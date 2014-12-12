#' Choosen criterions
#'
#' Checks if the criterion is viable or matches it to the list of implemented criterions.
#' 
#' @param input_criterion \code{character} string, criterion from input.
#' @param criterion_names list of implemented criterions, always in lowercase.
#' @export
#' @return a chosen criterion.

check_criterion <- function(input_criterion, criterion_names = c("ig")) {
  #think twice about grep
  res <- method.names[grepl(tolower(input_criterion), method.names)]
  
  if (length(res) == 0)
    stop("Name ", input_criterion, " cannot be associated with any available criterion.")
  
  if (length(res) > 1)
    stop("Name ", input_criterion, " is too ambiguous. Rerun with more precise name.")
  
  #TO DO - should also return the full name of criterion for purpose of summaries/plots
  res
}
