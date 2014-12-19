#' Permutation test for feature selection
#'
#' Performs a feature selection on positioned N-gram data using a Fisher's 
#' permutation test.
#'
#' @param target target vector.
#' @param features \code{integer} matrix of features with number of rows equal 
#' to the length of target vector.
#' @param criterion the criterion used in permutation test. See \code{\link{criterions}} for the
#' list of possible criterions.
#' @param quick logical, if \code{TRUE} Quick Permutation Test (QuiPT) is used.
#' @param times number of times procedure should be repetead. Ignored if \code{quick} is 
#' \code{TRUE}.
#' @details Currently implemented criterions:
#' \itemize{
#' \item{"\code{ig}" - information gain}
#' }
#' @return a vector of objects of htest class that relate to each feature tested
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso See \code{\link{criterion_distribution}} for insight on QuiPT.
#' @export
#' @keywords nonparametric
#' @references 
#' Radivojac P, Obradovic Z, Dunker AK, Vucetic S, 
#' \emph{Feature selection filters based on the permutation test} in 
#' Machine Learning: ECML 2004, 15th European 
#' Conference on Machine Learning, Springer, 2004.
#' @examples
#' tar_feat1 <- create_feature_target(10, 390, 0, 600) 
#' tar_feat2 <- create_feature_target(9, 391, 1, 599)
#' tar_feat3 <- create_feature_target(8, 392, 0, 600)
#' test_features(tar_feat1[,1], cbind(tar_feat1[,2], tar_feat2[,2], 
#' tar_feat3[,2]))
test_features <- function(target, features, criterion = "ig", quick = TRUE, times = 1e5) {
  
  valid_criterion <- check_criterion(criterion)
  
  #few tests for data consistency
  if (!all(target %in% c(0, 1))) {
    stop("target is not {0,1}-valued vector")
  }
  if (nrow(features) != length(target)) {
    stop("target and feature have different lengths")
  }
  
  apply(features, 2, function(feature) {
    if (!all(feature %in% c(0,1)) ) {
      stop("feature is not {0,1}-valued vector")
    }
  })
  
  if(quick) {
    
    # compute distribution once
    feature_size <- unique(if (class(features) == "simple_triplet_matrix") {
      col_sums(features)
    } else {
      colSums(features)
    })
    
    dists <- lapply(feature_size, function(i){
      t <- create_feature_target(i, sum(target) - i, 0, length(target) - sum(target)) 
      criterion_distribution(t[, 1], t[, 2], graphical.output = FALSE, criterion = criterion)
    })
    names(dists) <- feature_size
    
    apply(features, 2, function(feature) {
      feature <- as.matrix(feature, ncol = 1)
      n <- length(target)
      
      result <- NULL
      result[["estimate"]] <- valid_criterion[["crit_function"]](target = target, features = feature)
      dist <- dists[[paste(sum(feature))]]
      result[["p.value"]] <- 1 - dist[3, which.max(dist[1, ] >= result[["estimate"]] - 1e-15)]
      class(result) <- "htest"
      result[["estimate"]] <- "Target variable, feature variable"
      result[["method"]] <- "Information gain permutation test"
      names(result[["estimate"]]) <- "Criterion value for feature"
      result
    })
  } else {
    #slow version
    rowMeans(valid_criterion[["crit_function"]](target, features) <= 
               replicate(times, valid_criterion[["crit_function"]](sample(target), features)))
  }
}


