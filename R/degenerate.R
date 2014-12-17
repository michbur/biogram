#' Degenerate protein sequence
#'
#' 'Degenerates' amino acid or nucleic sequence by aggregating 
#' elements to bigger groups.
#' 
#' @param seq \code{character} vector representing single sequence.
#' @param aa_group list of groups to which elements of sequence should be aggregated.
#' @keywords manip
#' @return a \code{character} vector containing aggregated elements.
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
  if (!all(unique(seq) %in% unlist(aa_group)))
    stop("The sequence contains elements not present in any of groups.")
  
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}