#' Compute similarity index
#'
#' Computes similarity index between two encodings.
#' @inheritParams calc_ed
#' @return the value of similarity index.
#' @details Briefly, the similarity index is a fraction of elements that have 
#' the same pairing in both encodings. Pairing is a binary variable, that has 
#' value 1 if two elements are in the same group and 0 if not. For more details, 
#' see references.
#' @export
#' @seealso 
#' \code{\link{calc_ed}}: calculate the encoding distance between two encodings.
#' @references Stephenson, J.D., and Freeland, S.J. (2013). Unearthing the Root 
#' of Amino Acid Similarity. J Mol Evol 77, 159-169.
#' @examples
#' # example from Stephenson & Freeland, 2013 (Fig. 6)
#' enc1 <- list(`1` = "A",
#'              `2` = c("F", "E"),
#'              `3` = c("C", "D", "G"))
#' 
#' enc2 <- list(`1` = c("A", "G"),
#'              `2` = c("C", "D", "E", "F"))
#' 
#' enc3 <- list(`1` = c("D", "G"),
#'              `2` = c("E", "F"),
#'              `3` = c("A", "C"))
#'              
#' calc_si(enc1, enc2)
#' calc_si(enc2, enc3)
#' calc_si(enc1, enc3)

calc_si <- function(a, b) {
  if(!validate_encoding(a, unlist(b)))
    stop("Encodings ('a' and 'b') must contain the same elements.")
  
  # similarity matrix
  comp <- encoding2matrix(a) == encoding2matrix(b)
  diag(comp) <- 0
  
  sum(comp)/(nrow(comp) * (nrow(comp) - 1))
}

# converts encoding to a matrix as introduced by Stephenson 2013
encoding2matrix <- function(x) {
  x_df <- encoding2df(x, sort = TRUE)
  
  res <- sapply(sort(levels(x_df[["element"]])), function(single_element) {
    x_gr <- x_df[x_df[["element"]] == single_element, "gr_id"]
    x_df[["gr_id"]] == x_gr
  })
  diag(res) <- FALSE
  
  res
}