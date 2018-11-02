#' #' Simplify alphabets
#' #'
#' #' 'Degenerates' amino acid or nucleic sequence by aggregating 
#' #' elements to bigger groups.
#' #' 
#' #' @inheritParams count_ngrams
#' #' @inheritParams test_features
#' #' @param n_generations number of generations.
#' #' @param n_alph number of alphabets in a generation.
#' #' @param alph_len length of the alphabet.
#' #' @keywords manip
#' #' @return A list.
#' #' @export
#' 
#' 
#' simplify_alph <- function(target, seqs, n_generations = 100, n_alph = 200, alph_len = 6) {
#'   
#'   alphs <- lapply(1L:n_alph, function(dummy) {
#'     alph <- sample(1L:alph_len, size = 20, replace = TRUE)
#'     names(alph) <- biogram:::return_elements("prot")
#'     alph
#'   }) %>% lapply(fill_alph, alph_len)
#'   
#'   max_ftn <- c()
#'   mean_ftn <- c()
#'   best_alph <- list()
#'   
#'   
#'   for(i in 1L:n_generations) {
#'     alphs <- lapply(alphs, mutatate_alph, alph_len = alph_len)
#'     crossed <- cross_alph(alphs, alph_len = alph_len)
#'     new_gen <- c(crossed, sample(alphs, n_alph - length(crossed)))
#'     
#'     #target_vec <- rep(c(0, 1), sapply(feature_list, nrow))
#'     
#'     ftn <- 1 - sapply(new_gen, function(ith_alph) {
#'       
#'       #if(class(try(degenerate_ngrams(ngram_seqs, element_groups = simple2full(ith_alph)), silent = TRUE)) == "try-error")
#'       #  browser()
#'       
#'       deg_ngrams <- degenerate_ngrams(ngram_seqs, element_groups = simple2full(ith_alph))
#'       
#'       x <- colSums(deg_ngrams[target == 0, ])
#'       y <- colSums(deg_ngrams[target == 1, ])
#'       
#'       # cosine similarity
#'       as.vector(x %*% y / sqrt(x%*%x * y%*%y))
#'     })
#'     
#'     max_ftn <- c(max_ftn, max(ftn))
#'     mean_ftn <- c(mean_ftn, mean(ftn))
#'     best_alph <- c(best_alph, list(alphs[ftn == max(ftn)]))
#'     alphs <- new_gen[make_tournament(ftn)]
#'   }
#'   
#'   list(max_ftn = max_ftn,
#'        mean_ftn = mean_ftn, 
#'        best_alph = best_alph)
#' }
#' 
#' degenerate_ngrams <- function(x, element_groups) {
#'   
#'   deg_ngrams <- colnames(x) %>% 
#'     decode_ngrams %>% 
#'     strsplit("") %>% 
#'     lapply(degenerate, element_groups = element_groups) %>% 
#'     lapply(paste0, collapse = "") %>% 
#'     lapply(code_ngrams)
#'   
#'   res <- do.call(cbind, lapply(unique(deg_ngrams), function(ith_ngram) {
#'     row_sums(x[, ith_ngram == deg_ngrams, drop = FALSE])
#'   }))
#'   
#'   colnames(res) <- unique(deg_ngrams)
#'   
#'   res
#' }
#' 
#' # if the alphabet doesn't have defined number of groups
#' fill_alph <- function(x, alph_len = 4) {
#'   u_groups <- unique(x)
#'   while(length(u_groups) != alph_len) {
#'     missed_groups <- (1L:alph_len)[!(1L:alph_len %in% u_groups)]
#'     x[sample(1L:20, size = length(missed_groups), replace = FALSE)] <- missed_groups
#'     u_groups <- unique(x)
#'   } 
#'   x
#' }
#' 
#' # tmp <- rep(3, 20)
#' # fill_alph(tmp)
#' 
#' # mutation operator
#' mutatate_alph <- function(x, alph_len = 4, mutate_prob = 0.05) {
#'   changed <- runif(20) < mutate_prob
#'   x[changed] <- sample(1L:alph_len, size = sum(changed), replace = TRUE)
#'   fill_alph(x, alph_len)
#' }
#' 
#' # crossover two alphabets
#' cross_two_alph <- function(x, y) {
#'   chosen_group <- sample(x, size = 1)
#'   
#'   y[x == chosen_group] <- x[x == chosen_group]
#'   y
#' }
#' 
#' cross_alph <- function(alph_list, alph_len = 4, cross_prob = 0.5) {
#'   crossed <- runif(length(alph_list)) < cross_prob
#'   lapply(which(crossed), function(alph_id) {
#'     # second parent id
#'     sec_parent_id <- sample((1L:length(alph_list))[-alph_id], 1)
#'     fill_alph(cross_two_alph(alph_list[[alph_id]], alph_list[[sec_parent_id]]), alph_len)
#'   })
#' }
#' 
#' # set.seed(1)
#' # tmp_alphs <- lapply(1L:2, function(dummy)
#' #   sample(1L:4, size = 20, replace = TRUE)
#' # )
#' # 
#' # tmp_alphs
#' # cross_two_alph(tmp_alphs[[1]], tmp_alphs[[2]])
#' 
#' 
#' make_tournament <- function(ftn_vec) 
#'   sapply(1L:length(ftn_vec), function(ftn_id1) {
#'     # second parent id
#'     ftn_id2 <- sample((1L:length(ftn_vec))[-ftn_id1], 1)
#'     ifelse(ftn_vec[ftn_id1] > ftn_vec[ftn_id2], ftn_id1, ftn_id2)
#'   })
