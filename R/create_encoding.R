#' Create encoding
#'
#' Reduces an alphabet using physicochemical properties.
#' @param prop \code{matrix} of properties with number of column equal to the 
#' length of the alphabet. Column must be named after elements of the 
#' alphabet. Each row represents a different physicochemical property. 
#' @param len length of the resulting encoding. Must be larger than zero and 
#' smaller than number of elements in the alphabet.
#' @return An encoding.
#' @details The encoding is a list of groups to which elements of an alphabet 
#' should be reduced. All elements of the alphabet (all 
#' amino acids or all nucleotides) should appear in the encoding.
#' @export
#' @keywords hclust
#' @seealso 
#' \code{\link{calc_ed}}: calculate the encoding distance between two encodings.
#' \code{\link{encoding2df}}: converts an encoding to a data frame.
#' \code{\link{validate_encoding}}: validate a structure of an encoding.
#' @examples
#' enc1 = list(`1` = c("a", "t"), 
#'             `2` = c("g", "c"))
#' encoding2df(enc1)

create_encoding <- function(prop, len) {
  cl <- hclust(dist(t(prop)), method = "ward.D2")
  gr <- cutree(cl, k = len)
  
  agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
  #inside encodings, amino acids are ordered alphabetically
  agg_gr <- lapply(agg_gr, sort)
  #groups are sorted by their length
  agg_gr[order(lengths(agg_gr))]
}

#' Convert encoding to data frame
#'
#' Converts an encoding to a data frame.
#' @param x encoding.
#' @param sort if \code{TRUE} rows are sorted according to elements.
#' @return data frame with two columns. First column represents an index of a 
#' group in the supplied encoding and the second column contains all elements of 
#' the encoding.
#' @details The encoding is a list of groups to which elements of an alphabet 
#' should be reduced. All elements of the alphabet (all 
#' amino acids or all nucleotides) should appear in the encoding.
#' @export
#' @keywords manip
#' @seealso 
#' \code{\link{calc_ed}}: calculate the encoding distance between two encodings.
#' \code{\link{encoding2df}}: converts an encoding to a data frame.
#' \code{\link{validate_encoding}}: validate a structure of an encoding.
#' @examples
#' create_encoding(aaprop[1L:5, ], 5)

encoding2df <- function(x, sort = FALSE) {
  res <- do.call(rbind, lapply(1L:length(x), function(gr_id) {
    data.frame(gr_id = gr_id, element = x[[gr_id]])
  }))
  
  if(sort)
    res <- res[order(levels(res[["element"]])), ]
  
  res
}

#' Validate encoding
#'
#' Checks the structure of an encoding.
#' @param x encoding.
#' @param u \code{integer}, \code{numeric} or \code{character} vector of all
#' elements belonging to the encoding. See Details.
#' @return \code{TRUE} if the \code{x} is a correctly reduced \code{u}, 
#' \code{FALSE} in any other cases.
#' @details The encoding is a list of groups to which elements of an alphabet 
#' should be reduced. All elements of the alphabet (all 
#' amino acids or all nucleotides) should appear in the encoding.
#' @export
#' @keywords manip
#' @seealso 
#' \code{\link{calc_ed}}: calculate the encoding distance between two encodings.
#' \code{\link{encoding2df}}: converts an encoding to a data frame.
#' @importFrom stats cutree dist hclust
#' @examples
#' enc1 = list(`1` = c("a", "t"), 
#'             `2` = c("g", "c"))
#' # see if enc1 is the correctly reduced nucleotide (DNA) alphabet
#' validate_encoding(enc1, c("a", "c", "g", "t"))
#' 
#' # enc1 is not the RNA alphabet, so the results is FALSE
#' validate_encoding(enc1, c("a", "c", "g", "u"))
#' 
#' # validate_encoding works also on other notations
#' enc2 = list(a = c(1, 4),
#'             b = c(2, 3))
#' validate_encoding(enc2, 1L:4)

validate_encoding <- function(x, u) {
  if(!is.list(x))
    stop("'x' must have 'list' class.")
  all(sort(unlist(x)) == sort(u))
}