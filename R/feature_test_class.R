#' feature_test class
#'
#' A result of \code{\link{test_features}} function.
#'
#' @details An object of class \code{test_features} is a numeric vector of p-values. Additional 
#' attributes characterizes futher the details of test which returned these p-values.
#' @section Attributes:
#' \describe{
#'   \item{criterion}{the criterion used in permutation test.}
#'   \item{adjust}{the name of p-value adjusting method}
#'   \item{times}{the number of permutations or \code{NA} if QuiPT was chosen.}
#' }
#' @name feature_test
#' @docType class
NULL

create_feature_test <- function(p_value, criterion, adjust, times, occ) {
  if (!(is.numeric(p_value) && is.vector(p_value)))
    stop("p_values must be numeric")
  
  #add names if they are missing
  if(is.null(names(p_value)))
    names(p_value) <- paste0("feature", 1L:length(p_value))
  
  res <- p_value
  attributes(res) <- list(names = names(p_value),
                          criterion = criterion,
                          adjust = adjust,
                          times = times,
                          class = "feature_test",
                          occ = occ)
  res
}

#' Summarize tested features
#'
#' Summarizes results of \code{\link{test_features}} function.
#'
#' @param object of class \code{\link{feature_test}}.
#' @param conf_level confidence level. A feature with p-value equal to or smaller than the 
#' confidence is considered significant.
#' @param ... ignored
#' @return nothing.
#' @export
#' @keywords manip
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
#' @param frequencies a vector of frequencies along which occurences of n-grams belonging
#' to positive or negative group are classified.
#' @param split if not \code{NULL}, splits the output of the function. Must have value:
#' \code{NULL}, \code{"significances"}, \code{"positives"} and \code{"negatives"}. For
#' further information see Value
#' @param ... ignored
#' @return the value of function depends on the \code{split} parameter. 
#' If it is not \code{NULL}, the function returns a named list of length equal to the length 
#' of \code{significances} (when \code{split} equals \code{"significances"}) or 
#' \code{frequencies} (when \code{split} equals \code{"positives"} or \code{"negatives"})
#' minus one. Each elements of the list contains names of the n-grams belonging to the given 
#' significance or frequency group.
#' 
#' If the \code{split} parameter is \code{NULL}, the function returns data frame with four
#' columns containing respectively: names of n-grams, p-values, frequencies in positive
#' and negative group aggregated to groups along values of \code{significances} and 
#' \code{frequencies} vectors.
#' 
#' @export
#' @keywords manip
aggregate.feature_test <- function(x, significances = c(0, 0.0001, 0.01, 0.05, 1), 
                                   frequencies = c(0, 0.05, 0.1, 0.2, 1), 
                                   split = "significances", ...) {

  cutted_pvals <- cut(x, breaks = significances, include.lowest = TRUE)
  #aggregate does not cut here, because it does not return standard list output
  #dat <- aggregate(ngrams ~ cutted_pvals, data = data.frame(ngrams = names(x), cutted_pvals), 
  #                    function(i)
  #                      as.character(i))
  occ_pos <- cut(attr(x, "occ")["pos", ], 
                 breaks = frequencies, include.lowest = TRUE)
  occ_neg <- cut(attr(x, "occ")["neg", ], 
                 breaks = frequencies, include.lowest = TRUE)
  dat <- data.frame(ngram = names(x), 
                    p_value = cutted_pvals,
                    occ_pos = occ_pos,
                    occ_neg = occ_neg)
  
  if(!is.null(split)) {
    if(!(split %in% c("significances", "positives", "negatives")))
      stop("'split' must have one of following values: 'significances', 'positives', 
           'negatives' or NULL")
    
    split_factor <- switch(split,
                           significances = list(levels(cutted_pvals),
                                                "p_value"),
                           positives = list(levels(occ_pos),
                                            "occ_pos"),
                           negatives = list(levels(occ_neg),
                                            "occ_neg"))
    
    dat <- lapply(split_factor[[1]], function(i)
      as.character(dat[dat[[split_factor[[2]]]] == i, "ngram"]))
    
    names(dat) <- split_factor[[1]]
  }
  
  dat
}