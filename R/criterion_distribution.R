#' Compute criterion distribution 
#' 
#' Computes criterion distribution for feature, target under null hypothesis
#' 
#' @param target \{0,1\}-valued target vector. See Details.
#' @param feature \{0,1\}-valued feature vector. See Details.
#' @param graphical.output default value is \code{FALSE}, if \code{TRUE}
#'        probability density function is plotted
#' @param criterion the criterion used for calculations of distribution. 
#' See list of possible \code{\link{criterions}}.
#' @export
#' @details both \code{target} and \code{feature} vectors may contain only 0 and 1.
#' @return A matrix of 3 rows:
#' \describe{
#'   \item{1st row:}{possible values of criterion.}
#'   \item{2nd row:}{probability density function.}
#'   \item{3rd row:}{cumulative distribution function.}
#' }
#' @seealso \code{\link{criterions}}.
#' @keywords distribution
#' @examples
#' target_feature <- create_feature_target(10, 375, 15, 600) 
#' criterion_distribution(target = target_feature[,1], feature = target_feature[,2], 
#' graphical.output = TRUE)
criterion_distribution <- function(target, feature, graphical.output = FALSE, criterion = "ig") {
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
  
  non_zero_target <- sum(target)
  non_zero_feat <- sum(feature)
  p <- non_zero_target/n
  q <- non_zero_feat/n
  
  #values of criterion for different contingency tables
  diff_conts <- sapply(0L:min(non_zero_target, non_zero_feat), function(i) {
    #to do - check if other criterions also follow this distribution
    prob_log <- dmultinom(x = c(i, non_zero_feat - i, non_zero_target - i, 
                                n-non_zero_target - non_zero_feat + i),
                          size = n,
                          prob = c(p*q, (1-p)*q, p*(1-q), (1-p)*(1-q)),
                          log = TRUE)
    #feature-target data - different contingency tables
    ft_data <- create_feature_target(i, non_zero_feat - i, non_zero_target - i,
                                     n-non_zero_target - non_zero_feat + i)
    #values of criterion
    vals <- unname(valid_criterion[["crit_function"]](ft_data[,1], ft_data[, 2, drop = FALSE]))
    c(prob_log = prob_log, vals = vals)
  })
  
  dist_temp <- exp(diff_conts["prob_log", ])/sum(exp(diff_conts["prob_log", ]))
  if (graphical.output){
    #TO DO - remember that par manipulations changes pars for all plots in future. Revert old parameters
    #after plotting. Very clunky solution below.
    old_par <- par(c("mar", "fig", "oma"))
    par(mar = c(5,4,4,5) + 0.1)
    plot(0L:min(non_zero_target, non_zero_feat), diff_conts["vals", ], col="red", 
         xlab = "Number of cases with feature=1 and target=1",
         ylab = valid_criterion[["nice_name"]])
    par(new = TRUE)
    plot(0L:min(non_zero_target,non_zero_feat), dist_temp, type = "l", 
         col = "green", xaxt = "n", yaxt = "n",xlab = "",ylab = "")
    axis(4)
    mtext("density",side = 4,line = 3)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
        new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("top", legend = c(valid_criterion[["nice_name"]], "Probability"), xpd = TRUE, 
           horiz = TRUE, fill = c("red", "green"), bty = "n", cex = 1)
    par(mar = old_par[["mar"]], fig = old_par[["fig"]], oma = old_par[["oma"]])
  }
  
  # We get the same IG values for different contingency tables
  # therefore we need to combine them for distribution
  dist_temp <- dist_temp[order(diff_conts["vals", ])]
  val_temp <- diff_conts["vals", ][order(diff_conts["vals", ])]
  
  j <- 1
  criterion_distribution <- dist_temp[1]
  criterion_values <- val_temp[1]
  
  for(i in 2L:length(val_temp)) {
    #     if(is.na(abs(val_temp[i - 1] - val_temp[i]) < 1e-10))
    #       browser()
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
  
  dist <- rbind(criterion_values, 
                criterion_distribution, 
                1 - rev(cumsum(rev(criterion_distribution))))
  rownames(dist) <- c("criterion", "pdf", "cdf")
  dist
}