#' Create feature according to given contingency matrix
#'
#' Creates a matrix of features and target based on the values from contingency matrix.
#' 
#' @param n11 number of elements for which both target and feature equal 1.
#' @param n01 number of elements for which target and feature equal 1,0 
#' respectively.
#' @param n10 number of elements for which target and feature equal 0,1 
#' respectively.
#' @param n00 number of elements for which both target and feature equal 0.
#' @export
#' @return a matrix of 2 columns and n11+n10+n01+n00 rows. Columns represent
#' target and feature vectors, respectively.
#' @examples
#' #equivalent of 
#' #         target
#' #feature 10 375
#' #        15 600
#' target_feature <- create_feature_target(10, 375, 15, 600)
create_feature_target <- function(n11, n01, n10, n00){
  tar <- c(rep(1, n11), rep(1, n01), rep(0, n10), rep(0, n00))
  feat <- c(rep(1, n11), rep(0, n01), rep(1, n10), rep(0, n00))
  cbind(tar, feat)
}

#' Computes Information Gain distribution for feature, target under null hypothesis
#' 
#' @param target \{0,1\}-valued target vector. See Details.
#' @param feature \{0,1\}-valued feature vector. See Details.
#' @param graphical.output default value is \code{FALSE}, if \code{TRUE}
#'        probability density function is plotted
#' @param criterion the criterion used for calculations of distribution.
#' @export
#' @details both \code{target} and \code{feature} vectors may contain only 0 and 1.
#' @return A matrix of 3 rows:
#' \describe{
#'   \item{1st row:}{possible IG values.}
#'   \item{2nd row:}{probability density function.}
#'   \item{3rd row:}{cumulative distribution function.}
#' }
#' @seealso \code{\link{calc_ig}}
#' @examples
#' target_feature <- create_feature_target(10, 375, 15, 600) 
#' ig_distribution(target = target_feature[,1], feature = target_feature[,2], 
#' graphical.output = TRUE)
ig_distribution <- function(target, feature, graphical.output = FALSE, criterion = "ig") {
  n <- length(target)
  if (length(feature) != n) {
    stop("target and feature have different lengths")
  }
  if (!all(target %in% c(0, 1))) {
    stop("target is not {0,1}-valued vector")
  }
  if (!all(feature %in% c(0,1)) ) {
    stop("feature is not {0,1}-valued vector")
  }
  
  valid_criterion <- check_criterion(criterion)
  #here implement switch as external function
  if(valid_criterion == "ig")
    chosen_test <- calc_ig
  
  prob_log <- NULL
  ig <- NULL
  non_zero_target <- sum(target)
  non_zero_feat <- sum(feature)
  p <- non_zero_target/n
  q <- non_zero_feat/n
  for (i in 0:min(non_zero_target,non_zero_feat)){
    prob_log[i+1] <- 
      dmultinom(x = c(i, non_zero_feat-i, non_zero_target-i, 
                  n-non_zero_target-non_zero_feat+i),
                size = n,
                prob = c(p*q, (1-p)*q, p*(1-q), (1-p)*(1-q)),
                log = TRUE)
    dane <- create_feature_target(i, non_zero_feat-i, non_zero_target-i,
                                  n-non_zero_target-non_zero_feat+i)
    ig[i+1] <- chosen_test(dane[,1], dane[,2, drop=F])
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
  
  # We get the same IG values for different contingency tables
  # therefore we need to combine them for distribution
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

  dist <- rbind(ig_values, 
                ig_distribution, 
                1-rev(cumsum(rev(ig_distribution))))
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
#' @return a vector of objects of htest class that relate to each feature tested
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{calc_ig}}, \code{\link{ig_distribution}}
#' @export
#' @examples
#' tar_feat1 <- create_feature_target(10, 390, 0, 600) 
#' tar_feat2 <- create_feature_target(9, 391, 1, 599)
#' tar_feat3 <- create_feature_target(8, 392, 0, 600)
#' test_features_fast(tar_feat1[,1], cbind(tar_feat1[,2], tar_feat2[,2], 
#' tar_feat3[,2]))
test_features_fast <- function(target, features, criterion = "ig") {
  #TO DO - here we will go with switch
  if (criterion != "ig") {
    stop("Only Information Gain criterion is currently implemented")
  }

  apply(features, 2, function(feature) {
    if (length(feature) != length(target)) {
      stop("target and feature have different lengths")
    }
    if (!all(target %in% c(0, 1))) {
      stop("target is not {0,1}-valued vector")
    }
    if (!all(feature %in% c(0,1)) ) {
      stop("feature is not {0,1}-valued vector")
    }
  })
  # compute distribution once
  feature_size <- unique(colSums(features))
  dists <- lapply(feature_size, function(i){
    t <- create_feature_target(i, sum(target)-i, 0, length(target)-sum(target)) 
    return(i=ig_distribution(t[,1], t[,2], graphical.output = FALSE, criterion = criterion))
  })
  names(dists) <- feature_size
  
  apply(features, 2, function(feature) {
    feature <- as.matrix(feature, ncol=1)
    n <- length(target)
    
    result <- NULL
    result[["estimate"]] <- calc_ig(target = target, features = feature)
    dist <- dists[[paste(sum(feature))]]
    result[["p.value"]] <- 1 - dist[3, which.max(dist[1, ] >= result[["estimate"]] - 1e-15)]
    class(result) <- "htest"
    result[["estimate"]] <- "Target variable, feature variable"
    result[["method"]] <- "Information gain permutation test"
    names(result[["estimate"]]) <- "IG for feature"
    result
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
  #TO DO implement more criterions, do stuff like switch
  #add variant of https://github.com/michbur/chipPCR/blob/master/R/check.method.R
  if(criterion == "ig")
    chosen_test <- calc_ig
  rowMeans(chosen_test(target, features) <= 
             replicate(times, chosen_test(sample(target), features)))
}