#' Compute criterion distribution 
#' 
#' Computes criterion distribution under null hypothesis for all contingency 
#' tables possible for a feature and a target.
#' 
#' @param target \{0,1\}-valued target vector. See Details.
#' @param feature \{0,1\}-valued feature vector. See Details.
#' @param criterion criterion used for calculations of distribution. 
#' See \code{\link{calc_criterion}} for the list of avaible criteria.
#' @param iter_limit limit the number of calculated contingence matrices. If
#' \code{NULL}, computes all possible contingence matrices.
#' @export
#' @details both \code{target} and \code{feature} vectors may contain only 0 
#' and 1.
#' @return An object of class \code{\link{criterion_distribution}}.
#' @seealso \code{\link{calc_criterion}}.
#' @keywords distribution
#' @examples
#' target_feature <- create_feature_target(10, 375, 15, 600) 
#' distr_crit(target = target_feature[,1], feature = target_feature[,2])

distr_crit <- function(target, feature, criterion = "ig", iter_limit = 200) {
  n <- length(target)
  if (length(feature) != n) {
    stop("Target and feature have different lengths.")
  }
  if (!all(target %in% c(0, 1))) {
    stop("Target is not {0,1}-valued vector.")
  }
  if (!all(feature %in% c(0,1)) ) {
    stop("Feature is not {0,1}-valued vector.")
  }
  
  valid_criterion <- check_criterion(criterion)
  
  crit_function <- function(target, features)
    calc_criterion(target, features, valid_criterion[["crit_function"]])
  
  non_zero_target <- sum(target)
  non_zero_feat <- sum(feature)
  p <- non_zero_target/n
  q <- non_zero_feat/n
  
  # min_iter and max_iter limit the function to possible contingence matrices
  max_iter <- min(non_zero_target, non_zero_feat)
  min_iter <- max(0, non_zero_target + non_zero_feat - n)
  cross_tab <- fast_crosstable(target, length(target), sum(target), feature)
  # if(cross_tab[3L] == 0)
  #   max_iter <- sort(cross_tab)[2]

  if (is.null(iter_limit))
    iter_limit <- max_iter
  
  # values of criterion for different contingency tables
  crit_range <- max_iter - min_iter
  possible_crit_values <- if(crit_range > iter_limit) {
    round(seq(from = min_iter, to = max_iter, length.out = iter_limit), 0)
  } else {
    min_iter:max_iter
  }
  
  diff_conts <- sapply(possible_crit_values, function(i) {
    # to do - check if other criterions also follow this distribution
    
    k <- c(i, non_zero_feat - i, non_zero_target - i, 
           n - non_zero_target - non_zero_feat + i)
    
    prob_log <- dmultinom(x = k,
                          size = n,
                          prob = c(p*q, (1-p)*q, p*(1-q), (1-p)*(1-q)),
                          log = TRUE)
    # feature-target data - different contingency tables
    ft_data <- do.call(create_feature_target, as.list(k))
    # values of criterion
    vals <- unname(crit_function(ft_data[,1], ft_data[, 2, drop = FALSE]))
    c(prob_log = prob_log, vals = vals)
  })
  
  dist_temp <- exp(diff_conts["prob_log", ])/sum(exp(diff_conts["prob_log", ]))
  
  # We get the same IG values for different contingency tables
  # therefore we need to combine them for distribution
  dist_temp <- dist_temp[order(diff_conts["vals", ])]
  val_temp <- diff_conts["vals", ][order(diff_conts["vals", ])]
  
  j <- 1
  criterion_distribution <- dist_temp[1]
  criterion_values <- val_temp[1]
  
  if(length(val_temp) > 1)
    for(i in 2L:length(val_temp)) {
      if (abs(val_temp[i - 1] - val_temp[i]) < 1e-10) {
        criterion_values[j] <- criterion_values[j]
        criterion_distribution[j] <- criterion_distribution[j] + dist_temp[i]
      }
      else {
        j <- j + 1
        criterion_values[j] <- val_temp[i]
        criterion_distribution[j] <- dist_temp[i]
      }
    }
  
  create_criterion_distribution(criterion_values, 
                                criterion_distribution, 
                                possible_crit_values, 
                                diff_conts["vals", ],
                                exp(diff_conts["prob_log", ])/sum(exp(diff_conts["prob_log", ])),
                                valid_criterion[["nice_name"]])
}
