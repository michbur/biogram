#' Degenerate protein sequence
#'
#' 'Degenerates' amino acid or nucleic sequence by aggregating 
#' elements to bigger groups.
#' 
#' @param seq \code{character} vector or matrix representing single sequence.
#' @param element_groups encoding of elements: list of groups to which elements 
#' of sequence should be aggregated. Must have unique names.
#' @keywords manip
#' @return A \code{character} vector or matrix (if input is a matrix) 
#' containing aggregated elements.
#' @note  
#' Characters not present in the \code{element_groups} will be converted to NA with a 
#' warning.
#' @export
#' @seealso \code{\link{l2n}} to easily convert information stored in biological sequences from 
#' letters to numbers.
#' \code{\link{calc_ed}} to calculate distance between encodings.
#' @keywords manip
#' @examples
#' sample_seq <- c(1, 3, 1, 3, 4, 4, 3, 1, 2)
#' table(sample_seq)
#' 
#' # aggregate sequence to purins and pyrimidines
#' deg_seq <- degenerate(sample_seq, list(w = c(1, 4), s = c(2, 3)))
#' table(deg_seq)

degenerate <- function(seq, element_groups) {
  tmp_seq <- seq
  if (!all(unique(tmp_seq) %in% unlist(element_groups))) {
    warning("'seq' contains elements not present in 'element_groups'. 
            Element(s): ", 
            paste0(setdiff(unique(tmp_seq), unlist(element_groups)), collapse = ", "), 
            " will be replaced by NA.")
    tmp_seq[!(tmp_seq %in% unlist(element_groups))] <- NA
  }

  if(is.null(names(element_groups))) {
    warning("'element_groups' is unnamed. Assumed names of groups are their ordinal numbers.")
    names(element_groups) <- 1L:length(element_groups)
  }
  
  if(length(unique(names(element_groups))) != length(names(element_groups))) {
    stop("'element_groups' must have unique names.")
  }
    
  for (i in 1L:length(element_groups)) {
    tmp_seq[tmp_seq %in% element_groups[[i]]] <- names(element_groups)[i]
  }
  
  if(is.matrix(seq))
    dim(tmp_seq) <- dim(seq)
  
  tmp_seq
}

#' Convert letters to numbers
#'
#' Converts biological sequence from letter to number notation.
#' @inheritParams degenerate
#' @param seq_type the type of sequence. Can be \code{rna}, \code{dna} or \code{prot}.
#' @keywords manip
#' @return a \code{numeric} vector or matrix containing converted elements.
#' @export
#' @keywords manip
#' @seealso 
#' \code{l2n} is a wrapper around \code{\link{degenerate}}.
#' 
#' Inverse function: \code{\link{n2l}}.
#' @examples
#' sample_seq <- c("a", "d", "d", "g", "a", "g", "n", "a", "l")
#' l2n(sample_seq, "prot")

l2n <- function(seq, seq_type) {
  elements_list <- return_elements(seq_type)
  names(elements_list) <- 1L:length(elements_list)
  seq <- tolower(seq)
  deg_seq <- as.numeric(degenerate(seq, elements_list))
  if(is.matrix(seq))
    deg_seq <- matrix(deg_seq, ncol = ncol(seq))
  deg_seq
}


#' Convert numbers to letters
#'
#' Converts biological sequence from number to letter notation.
#' @param seq \code{integer} vector or matrix representing single sequence.
#' @param seq_type the type of sequence. Can be \code{rna}, \code{dna} or \code{prot}.
#' @keywords manip
#' @return a \code{character} vector or matrix containing converted elements.
#' @export
#' @keywords manip
#' @seealso 
#' \code{n2l} is a wrapper around \code{\link{degenerate}}.
#' 
#' Inverse function: \code{\link{l2n}}.
#' @examples
#' sample_seq <- c(1, 3, 3, 6, 1, 6, 12, 1, 10)
#' n2l(sample_seq, "prot")

n2l <- function(seq, seq_type) {
  names_list <- return_elements(seq_type)
  elements_list <- 1L:length(names_list)
  names(elements_list) <- names_list
  degenerate(seq, elements_list)
}

#' Convert encoding from full to simple format
#'
#' Converts an encoding from the full format to the simple format.
#' @param x encoding.
#' @export
#' @examples 
#' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#'            `2` = c("k", "h"), 
#'            `3` = c("d", "e"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' full2simple(aa1)
#' 
full2simple <- function(x) {
  single_enc <- x
  element_df <- do.call(rbind, lapply(1L:length(single_enc), function(i) {
    data.frame(gr = rep(names(single_enc[i]), length(single_enc[[i]])),
               element = single_enc[[i]], stringsAsFactors = FALSE)
  }))
  
  element_df <- element_df[order(element_df[["element"]]), ]
  res <- element_df[["gr"]]
  names(res) <- element_df[["element"]]
  res
}


#' Convert encoding from simple to full format
#'
#' Converts an encoding from the simple format to the full format.
#' @param x encoding (see Details).
#' @details The encoding should be named. Each name should correspond to a different
#' amino acid or nucleotide.
#' @export
#' @examples 
#' aa1 = structure(c("1", "4", "3", "3", "4", "1", "2", "1", "2", "1", 
#'                   "1", "4", "1", "4", "4", "4", "4", "1", "4", "4"), 
#'                 .Names = c("a", "c", "d", "e", "f", "g", "h", "i", 
#'                            "k", "l", "m", "n", "p", "q", 
#'                            "r", "s", "t", "v", "w", "y"))
#' simple2full(aa1)
#' 
simple2full <- function(x) {
  if(is.null(names(x)))
    stop("'x' must be named.")
  single_enc <- x
  gr <- unique(sort(single_enc))
  res <- lapply(gr, function(i)
    names(x[x == i]))
  names(res) <- gr
  res
}

# an internal function returning elements for a specific sequence type: aa, dna, rna
return_elements <- function(seq_type) {
  if (!(seq_type %in% c("prot", "dna", "rna")))
    stop("The value of 'what' must be: 'dna', 'rna' or 'prot'.")
  switch(seq_type,
         rna = c("a", "c", "g", "u"),
         dna = c("a", "c", "g", "t"),
         prot = c("a", "c", "d", "e", "f", 
                  "g", "h",  "i", "k", "l", 
                  "m", "n", "p", "q", "r", 
                  "s", "t", "v", "w", "y"))
}