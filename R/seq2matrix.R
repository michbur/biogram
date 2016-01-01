#' Convert list of sequences to matrix
#'
#' Converts list of sequences to matrix.
#' @param seq_list \code{list} of sequences, e.g. \code{\link[seqinr]{SeqFastaAA}}.
#' @keywords manip
#' @return A matrix with the number of rows equal to the number of sequences and the 
#' number of columns equal to the length of the longest sequence.
#' @note  
#' Since matrix must have specified number of columns, shorter sequences end with 
#' NA.
#' @export
#' @keywords manip
#' @examples
#' list2matrix(list(s1 = c("c", "g", "g", "t"),
#'                  s2 = c("g", "t", "c", "t", "t", "g")))

list2matrix <- function(seq_list) {
  max_len <- max(lengths(seq_list))
  t(sapply(seq_list, function(i)
    c(i, rep(NA, max_len) - length(i))
  ))
}
