#' feature_test class
#'
#' A result of \code{\link{test_features}} function.
#'
#' @details sth
#' @name feature_test
#' @docType class
NULL

#constructor
create_feature_test <- function(p_value, criterion, adjust, times) {
  if (!is.numeric(p_value)) 
    stop("p_values must be numeric")
  
  #add names if they are missing
  if(is.null(names(p_value)))
    names(p_value) <- paste0("feature", 1L:length(p_value))
  
  res <- p_value
  attributes(res) <- list(names = names(p_value),
                          criterion = criterion,
                          adjust = adjust,
                          times = times,
                          class = "feature_test")
  res
}

#' Summarize tested features
#'
#' Summarizes results of \code{\link{test_features}} function.
#'
#' @param object of class \code{\link{feature_test}}.
#' @param conf_level confidence level. A feature with p-value equal or smaller is considered 
#' significant.
#' @param ... ignored
#' @return nothing.
#' @export
#' @keywords print methods manip
summary.feature_test <- function(object, conf_level = 0.95, ...) {
  cat("Total number of features:", 
             length(object), "\n")
  cat("Number of significant features:", 
             sum(object <= 1 - conf_level), "\n")
  cat("Criterion used:", 
             attr(object, "criterion"), "\n")
  cat("Feature test:", 
             ifelse(is.na(attr(object, "times")), "QuiPT",
                    paste0("Fisher's permutation test (",  attr(object, "times"),
                           " permutations)")), "\n")
  cat("p-values adjustment method:", 
             attr(object, "adjust"), "\n")
  
}

#' Aggregate tested features
#'
#' Aggregates results of \code{\link{test_features}} function into groups based on their 
#' significance.
#'
#' @param x an object of class \code{\link{feature_test}}.
#' @param significances a vector of significances along which p-values are classified. See description of
#' \code{\link[base]{cut}} function. 
#' @param ... ignored
#' @return a named list of length equal to the length of \code{significances} minus one. Each elements of 
#' the list contains names of the n-grams belonging to the given significance group.
#' @export
#' @keywords print methods manip
aggregate.feature_test <- function(x, significances = c(0, 0.0001, 0.01, 0.05, 1), ...) {
  cutted_pvals <- cut(x, breaks = c(0, 0.0001, 0.01, 0.05, 1), include.lowest = TRUE)
  #aggregate does not cut here, because it does not return standard list output
  #dat <- aggregate(ngrams ~ cutted_pvals, data = data.frame(ngrams = names(x), cutted_pvals), 
  #                    function(i)
  #                      as.character(i))
  dat <- data.frame(ngrams = names(x), cutted_pvals)
  res <- lapply(levels(cutted_pvals), function(i)
    as.character(dat[dat[["cutted_pvals"]] == i, "ngrams"]))

  names(res) <- levels(cutted_pvals)
  res
}