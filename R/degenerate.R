#' Degenerate protein sequence
#'
#' 'Degenerates' protein sequence by aggregating aminoacids to bigger groups.
#' 
#' @param seq \code{character} vector representing single aminoacid sequence.
#' @param aa_group list of aminoacid groups to which sequence should be aggregated.
#' @keywords manip
#' @return a \code{character} vector.
#' @export
#' @keywords manip
#' @examples
#' sample_seq <- c(1, 3, 1, 3, 4, 4, 3, 1, 2)
#' table(sample_seq)
#' 
#' #aggregate sequence to purins and pyrimidines
#' deg_seq <- degenerate(sample_seq, list(W = c(1, 4), S = c(2, 3)))
#' table(deg_seq)

degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}