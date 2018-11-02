#' Read FASTA files
#'
#' A lightweight tool to read nucleic or amino-acid sequences from a file 
#' in FASTA format.
#' 
#' @param file the name of the file which the data are to be read from.
#' @return a list of sequences.
#' @export
#' @seealso 
#' \code{\link[seqinr]{read.fasta}}: heavier function for processing FASTA files.
#' @examples 
#' \dontrun{
#'   read_fasta("https://www.uniprot.org/uniprot/P28307.fasta")
#'   }

read_fasta <- function(file) {
  all_lines <- readLines(file)
  prot_id <- cumsum(grepl("^>", all_lines))
  all_prots <- split(all_lines, prot_id)

  seq_list <- lapply(all_prots, function(ith_seq)
    unlist(strsplit(ith_seq[-1], split = "")))
  
  names(seq_list) <- sub(">", "", sapply(all_prots, function(ith_seq) ith_seq[1]), fixed = TRUE)
  
  seq_list
}
