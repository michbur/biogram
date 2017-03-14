#' Write encodings to a file
#'
#' Saves a list of encodings (or a single encoding to the file).
#' @param x encoding or list of encodings.
#' @param file ither a character string naming a file or a 
#' \code{\link[base]{connection}} 
#' open for writing. "" indicates output to the console.
#' @export
#' @examples 
#' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#'            `2` = c("k", "h"), 
#'            `3` = c("d", "e"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' write_encoding(aa1)
#' 
write_encoding <- function(x, file = "") {
  single_enc <- x
  element_df <- do.call(rbind, lapply(1L:length(single_enc), function(i) {
    data.frame(gr = rep(names(single_enc[i]), length(single_enc[[i]])),
               element = single_enc[[i]], stringsAsFactors = FALSE)
  }))
  
  element_df <- element_df[order(element_df[["element"]]), ]
  cat(c("ID:",
        paste0("SG:", paste0(element_df[["gr"]], collapse = ";")),
        paste0("EL:", paste0(element_df[["element"]], collapse = ";"))),
      sep = "\n",
      file = file)
}