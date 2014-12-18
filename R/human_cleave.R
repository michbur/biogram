#' Human signal peptides cleavage sites
#' 
#' uniprot query \code{keyword:signal AND taxonomy:"Homo sapiens (Human) [9606]" 
#' AND annotation:(type:signal AND evidence:ECO_0000269) AND reviewed:yes}
#' cleaves[i] + 5):(cleaves[i] + 13)]
#' 
#' @name human_cleave
#' @docType data
#' @format A data frame with 1296 observations on the following 10 variables.
#' \describe{ 
#' \item{P1}{amino acid in first position.} 
#' \item{P9}{amino acid in nineth position.} 
#' \item{tar}{target vector (1 if cleavage site, 0 if post-cleavage site).}}
#' @keywords datasets
#' @examples
#' 
#' data(human_cleave)
#' table(human_cleave[, 1])
#' 
NULL