#' Regenerate n-grams
#'
#' 'Regenerates' amino acid or nucleic sequence written in a simplified 
#' alphabet by converting groups to regular expression.
#' 
#' @inheritParams is_ngram
#' @inheritParams degenerate
#' @keywords manip
#' @return A \code{character} string representing a POSIX regular expression.
#' @note  
#' Gaps (\code{_}) will be converted to any possible character from the alphabet
#' (nucleotides or amino acids).
#' @export
#' @seealso \code{\link{degenerate}} to easily convert information stored in biological sequences from 
#' letters to numbers.
#' \code{\link{calc_ed}} to calculate distance between simplified alphabets.
#' @keywords manip
#' @examples
#' regenerate("ssw", list(w = c(1, 4), s = c(2, 3)))

regenerate <- function(x, element_groups) {
  if (!all(unique(strsplit(x, "")[[1]]) %in% names(element_groups))) 
    stop("'x' contains groups not present in 'element_groups'")
  
  tmp_x <- x
  
  for(i in names(element_groups)) 
    tmp_x <- gsub(pattern = i, replacement = paste0(c("[", element_groups[[i]], "]"), collapse = ""), x = tmp_x)
  
  tmp_x
}
