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
  
  len_a != len_b
  
  #indices of groups to merge
  merges <- combn(1L:len_a, 2)
  
  #list of merged
  reduced_size <- lapply(1L:ncol(merges), function(single_merge)
    list(enc = c(list(unname(unlist(ta[merges[, single_merge]]))), ta[-merges[, single_merge]]),
         ed = c(min(lengths(ta[merges[, single_merge]])))))
  

  
  #rows b
  #columns a
  comp_tab <- create_comp_tab(ta, tb)
  
  #loop for moving aa between groups
  #TRUE if groups are not similiar
  #diff <- comp_tab != 0 & comp_tab != 1
  if(any(!(comp_tab %in% c(0, 1)))) {
    #ids of groups in ta that are the closest to groups in tb
    ta_gr <- apply(comp_tab, 1, which.max)
    
    gr_id <- 1
    
    for(i in 1L:length(tb)) {
      tb_group <- tb[[gr_id]]
      ida_to <- ta_gr[gr_id]
      ta_ids <- (1L:length(ta))[-ida_to]
      for (ida_from in ta_ids) {
        #id of the single amino acid
        el_id <- which(ta[[ida_from]] %in% tb_group)
        #add amino acid to second group
        ta[[ida_to]] <- c(ta[[ida_to]], ta[[ida_from]][el_id])
        #remove amino acid from the first group
        ta[[ida_from]] <-  ta[[ida_from]][-el_id]
        ed <- ed + length(el_id)
      }
    }
    #remove empty sets
    if(any(lengths(ta) == 0)) {
      ta <- ta[!lengths(ta) == 0]
    }
    
  }
  
  if(any(rowSums(comp_tab) != 1)) {
    #concatenate smaller groups into one big
    #check which row in comp_tab has two 1s 
    conc_group <- unname(do.call(c, ta[as.logical(comp_tab[rowSums(comp_tab) != 1, ])]))
    #get lengths of groups
    gr_len <- lengths(ta[as.logical(comp_tab[rowSums(comp_tab) != 1, ])])
    #choose number of elements in smaller group(s) as ed
    ed <- ed + sum(gr_len[-which.max(gr_len)])
    ta <- c(ta[!as.logical(comp_tab[rowSums(comp_tab) != 1, ])],
            list(conc_group))
  }
  # very verbose output, good for debugging 
#   list(ed = ed,
#        a = sapply(ta, sort),
#        b = sapply(tb, sort))
  ed
}

create_comp_tab <- function(ta, tb) {
  comp_tab <- sapply(ta, function(single_subgroup_a)
    sapply(tb, function(single_subgroup_b)
      sum(single_subgroup_a %in% single_subgroup_b)))
  if(class(comp_tab) != "matrix") {
    comp_tab <- matrix(comp_tab, nrow = length(tb), ncol = length(ta))
  }
  comp_tab
}

calc_ed_single <- function(ta, tb, ed = 0) {
  comp_tab <- create_comp_tab(reduced_size[[1]][["enc"]], tb)
  ta_order <- apply(comp_tab, 1, order, decreasing = TRUE)
  ta_gr <- rep(NA, length(tb))
  gr_id <- 1
  
  #occurence matrix
  data.frame(b = rep(1L:length(b), length(b)),
             position = sort(rep(1L:length(b), length(b))))
  
  #reasign a to match position
  data.frame(a = sort(rep(1L:length(b), length(b))),
             count = as.vector(comp_tab))
  
  
  for(i in 1L:length(tb)) {
    tb_group <- tb[[gr_id]]
    ida_to <- ta_gr[gr_id]
    ta_ids <- (1L:length(ta))[-ida_to]
    for (ida_from in ta_ids) {
      #id of the single amino acid
      el_id <- which(ta[[ida_from]] %in% tb_group)
      #add amino acid to second group
      ta[[ida_to]] <- c(ta[[ida_to]], ta[[ida_from]][el_id])
      #remove amino acid from the first group
      ta[[ida_from]] <-  ta[[ida_from]][-el_id]
      ed <- ed + length(el_id)
    }
  }
}


