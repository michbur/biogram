#' Tabulate n-grams
#' 
#' Tabulates counts of n-grams versus their class labels.
#' @param x a \code{character} string representing single n-gram.
#' @return \code{TRUE} if n-gram's notation is correct, \code{FALSE} if not.
#' @export
#' @examples
#' print(is_ngram("1_1.1.1_0.0"))
#' print(is_ngram("not_ngram"))



all_ngrams <- seq2ngrams(deg_seqs, 3, letters[1L:5], pos = TRUE)

tar <- human_cleave[, "tar"]

ngrams <- gr[[1]]

table_ngrams 

all_ngrams <- as.matrix(count_specified(deg_seqs, ngrams))
#counts of each value of target
n_tar <- tabulate(tar + 1)

#values of target
val_tar <- sort(unique(tar))

res <- t(vapply(ngrams, function(ngram_name) 
  vapply(val_tar, function(target_value) sum(all_ngrams[tar == target_value, ngram_name]), 0)/n_tar, c(0, 0)))

res_df <- data.frame(rownames(res), res)
rownames(res_df) <- NULL
colnames(res_df) <- c("ngram", paste0("target", val_tar))

as.data.frame(table(all_ngrams[, 1], human_cleave[, "tar"]))
