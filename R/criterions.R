#' Critetions
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
#' The information gain term is used here (improperly) as a synonym of mutual information. 
#' It is defined as:
#' \deqn{IG(X; Y) = \sum_{y \in Y} \sum_{x \in X} p(x, y) \log \left(\frac{p(x, y)}{p(x) p(y)}  \right)}
#' 
#' In biogram package information gain is calculated using following relationship: 
#' \eqn{IG = E(S) - E(S|F)}
#' @references 
#' Cover TM, Thomas JA \emph{Elements of Information Theory, 2nd Edition}
#' Wiley, 2006.
#' @name criterions
NULL