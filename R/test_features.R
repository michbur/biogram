#'Calculate log of partial binomial probabilities 
#'
#'@param p probability of success
#'@param m number of trials
#'@param n number of successes
#'
#'@return temp numeric, log(m!/[(m-n)! n!] p^n
#'@keywords internal
calc_partial_binomial_log_probability <- function(p, m, n) { 
  temp = lgamma(m + 1.0);
  temp = temp - lgamma(n + 1.0) - lgamma(m - n + 1.0);
  temp = temp + n*log(p);
  temp
}

#' Calculate the probability of a certain contingency table
#' 
#' Contingency table for two Bernoulli random variables. Feature and
#' target. 
#' 
#' @param p numeric, a probability of feature equaling 1
#' @param q numeric, a probability of feature equaling 1
#' @param n11 number of elements for which both target and feature equal 1
#' @param n01 number of elements for which both target and feature equal 1,0
#' @param n10 number of elements for which both target and feature equal 0,1
#' @param n00 number of elements for which both target and feature equal 0
#' 
#' @keywords internal
#' @return log_probability probability of getting given contingency table
#'         under p and q
calc_cont_table_log_probability <- function(p, q, n11, n01, n10, n00) {
  n <- n11 + n01 + n10 + n00
  log_probability <- calc_partial_binomial_log_probability(p*q, n, n11)
  log_probability <- log_probability + 
    calc_partial_binomial_log_probability((1-p)*q, n-n11, n01)
  log_probability <- log_probability + 
    calc_partial_binomial_log_probability(p*(1-q), n-n11-n01, n10)
  log_probability <- log_probability + 
    calc_partial_binomial_log_probability((1-p)*(1-q), n-n11-n01-n10, n00)
  return(log_probability)
}

#' Create feature with given contingency matrix
#'
#' @param n11 number of elements for which both target and feature equal 1
#' @param n01 number of elements for which both target and feature equal 1,0
#' @param n10 number of elements for which both target and feature equal 0,1
#' @param n00 number of elements for which both target and feature equal 0
#' @export
#' @return a matrix of 2 columns and n11+n10+n01+n00 rows. First column is
#' target, second is feature
#' @examples
#' temp <- create_feature_target(10, 375, 15, 600)
create_feature_target <- function(n11, n01, n10, n00){
  tar <- rep(1, n11)
  feat <- rep(1, n11)
  tar <- c(tar, rep(1, n01))
  feat <- c(feat, rep(0, n01))
  tar <- c(tar, rep(0, n10))
  feat <- c(feat, rep(1, n10))
  tar <- c(tar, rep(0, n00))
  feat <- c(feat, rep(0, n00))
  return(cbind(tar, feat))
}

#' Computes Information Gain distribution for feature, target under null hypothesis
#' 
#' @param target {0,1}-valued target vector
#' @param feature {0,1}-valued feature vector
#' @param graphical.output default value is FALSE
#' @export
#' @return a matrix of 3 rows. First row are possible IG values, 
#'  second in probability density function, third is cimulative 
#'  distribution function
#'  @examples
#' temp <- create_feature_target(10, 375, 15, 600) 
#' ig_distribution(temp[,1], temp[,2], graphical.output = TRUE)
ig_distribution <- function(target, feature, graphical.output = FALSE) {
  n <- length(target)
  if (length(feature) != n) {
    stop("target and feature have different lengths")
  }
  if (length(target[target %in% c(0,1)]) != n ) {
    stop("target is not {0,1}-valued vector")
  }
  if (length(feature[feature %in% c(0,1)]) != n ) {
    stop("feature is not {0,1}-valued vector")
  }
  
  prob_log <- NULL
  ig <- NULL
  non_zero_target <- sum(target)
  non_zero_feat <- sum(feature)
  p <- non_zero_target/n
  q <- non_zero_feat/n
  for (i in 0:min(non_zero_target,non_zero_feat)){
    prob_log[i+1] <- 
      calc_cont_table_log_probability(p, q, i, non_zero_feat-i, 
                                      non_zero_target-i, 
                                      n-non_zero_target-non_zero_feat+i)
    dane <- create_feature_target(i, non_zero_feat-i, non_zero_target-i,
                                  n-non_zero_target-non_zero_feat+i)
    ig[i+1] <- calc_ig(dane[,1], dane[,2, drop=F])
  }
  
  ig_dist_temp <- exp(prob_log)/sum(exp(prob_log))
  if (graphical.output){
    par(mar=c(5,4,4,5)+0.1)
    plot(0:min(non_zero_target,non_zero_feat), ig, col="red", 
         xlab="Number of cases with feature=1 and target=1",
         ylab="Information gain")
    par(new=TRUE)
    plot(0:min(non_zero_target,non_zero_feat), ig_dist_temp, type="l", 
         col="green", xaxt="n", yaxt="n",xlab="",ylab="")
    axis(4)
    mtext("density",side=4,line=3)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
        new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("top", legend = c("Information gain", "Probability"), xpd = TRUE, 
           horiz = TRUE, fill=c("red", "green"),  bty = "n", cex = 1)
  }
  
  ig_dist_temp <- ig_dist_temp[order(ig)]
  ig_val_temp <- ig[order(ig)]
  j <- 1
  ig_distribution <- ig_dist_temp[1]
  ig_values <- ig_val_temp[1]
  for(i in 2:length(ig_val_temp)) {
    if (abs(ig_val_temp[i-1]-ig_val_temp[i])<1e-10) {
      ig_values[j] <- ig_values[j]
      ig_distribution[j] <- ig_distribution[j] + ig_dist_temp[i]
    }
    else {
      j <- j + 1
      ig_values[j] <- ig_val_temp[i]
      ig_distribution[j] <- ig_dist_temp[i]
    }
  }

  dist <- rbind(ig_values, ig_distribution, cumsum(ig_distribution))
  rownames(dist) <- c("IG", "pdf", "cdf")
  return(dist)
}


#' Permutation test for feature selection
#'
#' Performs a feature selection on positioned N-gram data using a Fisher's 
#' permutation test.
#'
#' @inheritParams calc_ig
#' @param criterion the criterion used in permutation test.
#' @details Currently implemented criterions:
#' \itemize{
#' \item{"\code{ig}" - information gain}
#' }
#' @return a numeric vector of lenth equal to the number of features containing computed
#' information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{calc_ig}}, \code{\link{ig_distribution}}
#' @export
#' @examples
#' temp <- create_feature_target(10, 390, 0, 600) 
#' temp2 <- create_feature_target(9, 391, 1, 599) 
#' test_features_fast(temp[,1], cbind(temp[,2], temp[,2]))
test_features_fast <- function(target, features, criterion = "ig") {
  apply(features, 2, function(feature) {
    feature <- as.matrix(feature, ncol=1)
    n <- length(target)
    if (length(feature) != n) {
      stop("target and feature have different lengths")
    }
    if (length(target[target %in% c(0,1)]) != n ) {
      stop("target is not {0,1}-valued vector")
    }
    if (length(feature[feature %in% c(0,1)]) != n ) {
      stop("feature is not {0,1}-valued vector")
    }
    if(criterion != "ig")
      stop("Only Information Gain criterion is avaialble")
    result <- NULL
    result$estimate <- calc_ig(target = target, features = feature)
    dist <- ig_distribution(target, feature, graphical.output = FALSE)
    result$p.value <- 1-dist[3, which.max(dist[1,]>=result$estimate-1e-15)]
    class(result) <- "htest"
    result$data.name <- "Target variable, feature variable"
    result$method <- "Information gain permutation test"
    names(result$estimate) <- "IG for feature"
    return(result)
  })
}


#' Permutation test for feature selection
#'
#' Performs a feature selection on positioned N-gram data using a Fisher's 
#' permutation test.
#'
#' @inheritParams calc_ig
#' @param times number of times procedure should be repetead
#' @param criterion the criterion used in permutation test.
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