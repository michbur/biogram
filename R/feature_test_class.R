#' feature_test class
#'
#' A result of \code{\link{test_features}} function.
#'
#' @details sth
#' @name feature_test
NULL

#constructor
create_feature_test <- function(p_value, criterion, adjust, times) {
  if (!is.numeric(p_value)) 
    stop("p_values must be numeric")
  #other tests should be here
  structure(list(p_value = p_value, 
                 criterion = criterion,
                 adjust = adjust,
                 times = times), class = "feature_test")
}


#' Summarize tested features
#'
#' Summarizes results of \code{\link{test_features}} function.
#'
#' @param object of class \code{\link{feature_test}}.
#' @param significance a threshold level. A feature with p-value equal or below it is considered 
#' significant.
#' @param ... ignored
#' @return nothing.
#' @export
#' @keywords print methods manip

summary.feature_test <- function(object, significance = 0.05, ...) {
  cat(paste0("Total number of features: ", 
             length(object[["p_value"]]), "\n"))
  cat(paste0("Number of significant features: ", 
             sum(object[["p_value"]] <= significance), "\n"))
  cat(paste0("Criterion used: ", 
             object[["criterion"]], "\n"))
  cat(paste0("p-values adjustment method: ", 
             object[["adjust"]], "\n"))
}