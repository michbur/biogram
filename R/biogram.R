#' biogram - analysis of n-grams from biological sequences
#'
#' @description biogram specialises in analysis of n-grams (k-tuples) from nucleotide or
#' protein sequences.
#' @import slam
#' @importFrom bit as.bit
#' @author Michal Burdukiewicz, Piotr Sobczyk, Chris Lauber
#' @docType package
#' @name biogram-package
#' @aliases biogram
#' @examples
#' #use data set from package
#' data(human_cleave)
#' #first nine columns represent subsequent nine amino acids from cleavage sites
#' #degenerate the sequence to reduce the dimensionality of the problem
#' deg_seqs <- degenerate(human_cleave[, 1L:9], 
#'                       list(`1` = c(1, 6, 8, 10, 11, 18), 
#'                            `2` = c(2, 13, 14, 16, 17), 
#'                            `3` = c(5, 19, 20), 
#'                            `4` = c(7, 9, 12, 15), 
#'                            '5' = c(3, 4)))
#' #extract bigrams
#' bigrams <- count_ngrams(deg_seqs, 3, 1L:4, pos = TRUE)
NULL