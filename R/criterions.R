#' Possible critetions
#'
#' Permutation test implemented in \code{biogram} uses several criterions to filter 
#' important features. Each can be used by \code{\link{test_features}} by specifying
#' \code{criterion} parameter.
#' 
#' Possible criterions are:
#' \describe{
#'   \item{ig}{Information Gain. Calculated using \code{\link{calc_ig}}.}
#' }
#' 
#' @section Information gain:
#' 
#' TO DO here some words about IG.
#' \deqn{IG = E(S) - E(S|F)}
#' \deqn{E(S) = -\frac{N_{pos}}{N} \log \frac{N_{pos}}{N} - \frac{N_{neg}}{N} \log \frac{N_{neg}}{N}}
#' \deqn{E(S|F) = -\frac{N_{f+}}{N} \left( \frac{N_{pos,f+}}{N_{f+}} \log \frac{N_{pos,f+}}{N_{f+}}+ \frac{N_{neg,f+}}{N_{f+}} \log \frac{N_{neg,f+}}{N_{f+}} \right) - \frac{N_{f-}}{N} \left( \frac{N_{pos,f-}}{N_{f-}} \log \frac{N_{pos,f-}}{N_{f-}}  + \frac{N_{neg,f-}}{N_{f-}} \log \frac{N_{neg,f-}}{N_{f-}} \right)}
#' @references Cover TM, Thomas JA \emph{Elements of Information Theory, 2nd Edition}
#' Wiley, 2006.
#' @name criterions
NULL