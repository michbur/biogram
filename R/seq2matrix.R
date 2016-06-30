#' Convert list of sequences to matrix
#'
#' Converts list of sequences to matrix.
#' @param seq_list list of sequences (e.g. as returned by 
#' the \code{\link[seqinr]{read.fasta}} function).
#' @keywords manip
#' @return A matrix with the number of rows equal to the number of sequences and the 
#' number of columns equal to the length of the longest sequence.
#' @note  
#' Since matrix must have specified number of columns, ends of shorter sequences are 
#' completed with NAs.
#' @export
#' @keywords manip
#' @examples
#' list2matrix(list(s1 = c("c", "g", "g", "t"),
#'                  s2 = c("g", "t", "c", "t", "t", "g"),
#'                  s3 = c("a", "a", "t")))

list2matrix <- function(seq_list) {
  max_len <- max(lengths(seq_list))
  do.call(rbind, lapply(seq_list, function(i)
    c(i, rep(NA, max_len - length(i)))
  ))
}

