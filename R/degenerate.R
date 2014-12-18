#' Degenerate protein sequence
#'
#' 'Degenerates' amino acid or nucleic sequence by aggregating 
#' elements to bigger groups.
#' 
#' @param seq \code{character} vector representing single sequence.
#' @param aa_group list of groups to which elements of sequence should be aggregated.
#' @keywords manip
#' @return a \code{character} vector containing aggregated elements.
#' @note Both sequence and \code{aa_group} should contain lower-case letters.
#' Upper-case will be automatically converted with a message.
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
  seq <- tolower(seq)
  if (!all(unique(seq) %in% unlist(aa_group)))
    stop("The sequence contains elements not present in any of groups.")
  
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}

#' Convert letters to numbers
#'
#' Converts biological sequence of letters to number notation.
#' @param seq \code{character} vector representing single sequence.
#' @param what the type of sequence. Can be \code{rna}, \code{dna} or \code{prot}.
#' @keywords manip
#' @return a \code{numeric} vector containing converted elements.
#' @export
#' @keywords manip
#' @examples
#' sample_seq <- c("a", "d", "d", "g", "a", "g", "n", "a", "l")
#' letters2numbers(sample_seq, "prot")

letters2numbers <- function(seq, what = "rna") {
  elements_list <- switch(what,
                          rna = c("a", "c", "g", "u"),
                          dna = c("a", "c", "g", "t"),
                          prot = c("a", "c", "d", "e", "f", 
                                   "g", "h",  "i", "k", "l", 
                                   "m", "n", "p", "q", "r", 
                                   "s", "t", "v", "w", "y"))
  names(elements_list) <- 1L:length(elements_list)
  as.numeric(degenerate(seq, elements_list))
}