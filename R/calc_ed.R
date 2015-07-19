#' Calculate encoding distance
#' 
#' Compares two encodings and calculates encoding distance between them.
#' @param a list of groups to which elements of sequence should be aggregated.
#' @param b list of groups to which \code{a} should be compared. Must be shorter than \code{a} 
#' or have equal length.
#' @return an encoding distance.
#' @export
#' @examples
#' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#'            `2` = c("k", "h"), 
#'            `3` = c("d", "e"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' 
#' aa2 = list(`1` = c("g", "a", "p", "v", "m", "l", "q"), 
#'            `2` = c("k", "h", "d", "e", "i"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n"))
#' calc_ed(aa1, aa2)     
#'     

calc_ed <- function(a, b) {
  #compare temporary a to temporary b
  if(length(a) < length(b)) {
    warning("'b' must be longer than 'a'. Reverting a and b.") 
    tb <- a
    ta <- b
  } else {
    ta <- a
    tb <- b
  }
  
  if(any(lengths(ta) == 0))
    ta <- ta[lengths(ta) != 0]
  
  if(any(lengths(tb) == 0))
    tb <- tb[lengths(tb) != 0]
  
  if(length(unlist(a)) != length(unlist(b)))
    stop("Encodings ('a' and 'b') must contain the same number of elements.")
  
  if(!all(sort(unlist(a)) == sort(unlist(b))))
    stop("Encodings ('a' and 'b') must contain the same elements.")
  
  #encoding distance - distance between encodings
  ed <- 0
  
  #exclude identical subgroups in a and b, because hey do not affect ed.
  ident_gr <- which(sapply(ta, function(single_subgroup_a)
    sapply(tb, function(single_subgroup_b)
      identical(sort(single_subgroup_a), sort(single_subgroup_b)))), arr.ind = TRUE)
  
  if(length(ident_gr) != 0) {
    ta <- ta[-ident_gr[, "col"]]
    tb <- tb[-ident_gr[, "row"]]
  }
  
  #if encodings are identical, ta and tb are 0
  if(length(ta) == 0 && length(tb) == 0)
    return(ed)
  
  len_a <- length(ta)
  len_b <- length(tb)
  
  ed <- if(len_a != len_b) {
    if(len_b == 1) {
      sum(lengths(ta)[lengths(ta) != max(lengths(ta))])
    } else {
      #indices of groups to merge
      merges <- combn(1L:len_a, 2)
      
      #list of merged
      reduced_size <- lapply(1L:ncol(merges), function(single_merge)
        list(enc = c(list(unname(unlist(ta[merges[, single_merge]]))), ta[-merges[, single_merge]]),
             ed = c(min(lengths(ta[merges[, single_merge]])))))
      
      min(sapply(reduced_size, function(single_ta) calc_ed_single(single_ta[["enc"]], tb,
                                                                  single_ta[["ed"]])))
    }
  } else {
    calc_ed_single(ta, tb)
  }
  
  unname(ed)
}

calc_ed_single <- function(ta, tb, ed = 0) {
  
  # comparision table, for case when I have enough time to write something 
  # smart, not all permutations 
  #   comp_tab <- cbind(bgr = sort(rep(1L:length(ta), length(tb))),
  #                     do.call(rbind, lapply(tb, function(single_subgroup_b) {
  #                       counts <- do.call(rbind, lapply(ta, function(single_subgroup_a)
  #                         data.frame(counts = sum(single_subgroup_a %in% single_subgroup_b))))
  #                       
  #                       #order matrix
  #                       #first row - sequence of groups
  #                       #second row - groups' order
  #                       ordmat <- matrix(c(1L:length(tb),order(counts[, "counts"], decreasing = TRUE)), 
  #                                        nrow = 2, byrow = TRUE)
  #                       cbind(agr = 1L:length(tb), position = ordmat[1, order(ordmat[2, ])],
  #                             counts)
  #                     })
  #                     ))
  #   
  #   agr <- comp_tab[comp_tab[, "position"] == 1, "agr"]
  #   bgr <- comp_tab[comp_tab[, "position"] == 1, "bgr"]
  
  
  len_b <- length(tb)
  perms <- expand.grid(lapply(1L:len_b, function(dummy) 1L:len_b))
  perms <- perms[apply(perms, 1, function(single_permutation) 
    length(unique(single_permutation))) == len_b, ]
  bgr <- 1L:len_b
  
  
  ed + min(apply(perms, 1, function(agr) {
    sum(lengths(lapply(1L:len_b, function(gr_id) {
      tb_group <- tb[[bgr[[gr_id]]]]
      #added amino acids from other groups
      unlist(lapply(agr[-gr_id], function(ida_from) {
        #id of the single amino acid
        ta[[ida_from]][which(ta[[ida_from]] %in% tb_group)]
      }))
    })))
  }))
}




