# Check chosen criterion
# 
# Checks if the criterion is viable or matches it to the list of implemented 
# criterions.
# 
# @param input_criterion \code{character} string, criterion from input.
# @param criterion_names list of implemented criterions, always in lowercase.
# @export
# @return a list of three: 
# \itemize{
# \item{criterion name,}
# \item{its function,}
# \item{nice name for outputs.}
# }
# @seealso
# Calculate the value of criterion: \code{\link{calc_criterion}}.
check_criterion <- function(input_criterion, criterion_names = c("ig", "kl", "cs")) {
  # think twice about grep
  valid_name <- criterion_names[grepl(tolower(input_criterion), criterion_names)]
  
  if (length(valid_name) == 0)
    stop("Name ", input_criterion, " cannot be associated with any available criterion.")
  
  if (length(valid_name) > 1)
    stop("Name ", input_criterion, " is too ambiguous. Rerun with more precise name.")
  
  
  criterion_data <- switch(valid_name,
                           ig = list(crit_function = calc_ig, nice_name = "Information Gain"),
                           kl = list(crit_function = calc_kl, nice_name = "Kullback-Leibler divergence"),
                           cs = list(crit_function = calc_cs, nice_name = "Chi-squared-based measure"))
  # TO DO - should also return the full name of criterion for purpose of summaries/plots
  c(crit_name = valid_name, criterion_data)
}

#' Calculate value of criterion
#'
#' Computes a chosen statistical criterion for each feature versus target vector.
#'
#' @details The permutation test implemented in \code{biogram} uses several criterions to filter 
#' important features. Each can be used by \code{\link{test_features}} by specifying the 
#' \code{criterion} parameter.
#' 
#' @inheritParams test_features
#' @param criterion_function a function calculating criterion. For a full list, see 
#' \code{\link{test_features}}.
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{test_features}}.
#' @export
#' @examples 
#' tar <- sample(0L:1, 100, replace = TRUE)
#' feats <- matrix(sample(0L:1, 400, replace = TRUE), ncol = 4)
#' 
#' # Information Gain
#' calc_criterion(tar, feats, calc_ig)
#' 
#' # hi-squared-based measure
#' calc_criterion(tar, feats, calc_cs)
#' 
#' # Kullback-Leibler divergence
#' calc_criterion(tar, feats, calc_kl)
calc_criterion <- function(target, features, criterion_function) {
  tar_bit <- as.bit(target)
  l_tar <- length(target)
  pos_tar <- sum(target)
  
  apply(features, 2, function(single_feature) 
    criterion_function(single_feature, tar_bit, l_tar, pos_tar))
}