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
  nonempty_lines <- all_lines[which(all_lines != "")]
  prot_id <- cumsum(grepl("^>", nonempty_lines))
  all_prots <- split(nonempty_lines, prot_id)

  seq_list <- lapply(all_prots, function(ith_seq)
    unlist(strsplit(ith_seq[-1], split = "")))
  
  names(seq_list) <- sub(">", "", sapply(all_prots, function(ith_seq) ith_seq[1]), fixed = TRUE)
  
  seq_list
}


#' Write FASTA files
#'
#' A lightweight tool to read nucleic or amino-acid sequences from a file 
#' in FASTA format.
#' 
#' @param seq a list of sequences.
#' @param file the name of the output file.
#' @param nchar the number of characters per line.
#' @export
#' @seealso 
#' \code{\link[seqinr]{write.fasta}}: heavier function for writing FASTA files.
write_fasta <- function(seq, file, nchar = 80) {
  char_vec <- unlist(lapply(1L:length(seq), function(ith_id) 
    c(paste0(">", names(seq[ith_id])), 
      lapply(split(seq[[ith_id]], floor(seq_along(seq[[ith_id]])/nchar)), paste0, collapse = ""))))
  writeLines(text = char_vec, con = file)
}


#' Read sequences from .txt file
#'
#' Read sequence data saved in text file.
#'
#' @param connection a \code{\link{connection}} to the text (.txt) file.
#' @return a list of sequences. 
#' @details The input file should contain one or more amino acid sequences separated by 
#' empty line(s).
#' @export
#' @keywords manip
#' @examples 
#' sequences <- read_txt(system.file("PlastoGram/sequences.txt", package = "PlastoGram"))

read_txt <- function(connection) {
  
  content <- readLines(connection)
  
  #test for empty content
  if(content[1] != "" || length(content) > 1) {
    if (sum(grepl(">", content, fixed = TRUE)) == 0) {
      if (content[1] != "")
        content <- c("", content)
      
      #number of empty lines
      nel <- 0
      #content without too many empty lines
      content2 <- c()
      for (i in 1L:length(content)) {
        if(content[i] == "") {
          nel <- nel + 1
        } else {
          nel <- 0
        }
        if (nel <= 1)
          content2 <- c(content2, content[i])
      }
      content <- content2
      content_end <- length(content)
      while(content[content_end] == "i")
        content_end <- content_end - 1
      prot_names <- sapply(1L:sum(content == ""), function(i)
        paste0(">sequence", i))
      content[content == ""] <- prot_names
    }
    read_fasta(textConnection(content))
  } else {
    warning("No text detected.")
    NULL
  } 
}
