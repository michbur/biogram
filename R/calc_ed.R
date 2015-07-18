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
  
  min(sapply(reduced_size, function(single_ta) calc_ed_single(single_ta[["enc"]], tb,
                                                          single_ta[["ed"]])))
}

calc_ed_single <- function(ta, tb, ed = 0) {
  
  comp_tab <- cbind(bgr = sort(rep(1L:length(ta), length(tb))),
                    do.call(rbind, lapply(tb, function(single_subgroup_b) {
                      counts <- do.call(rbind, lapply(ta, function(single_subgroup_a)
                        data.frame(counts = sum(single_subgroup_a %in% single_subgroup_b))))
                      
                      #order matrix
                      #first row - sequence of groups
                      #second row - groups' order
                      ordmat <- matrix(c(1L:length(tb),order(counts[, "counts"], decreasing = TRUE)), 
                                       nrow = 2, byrow = TRUE)
                      cbind(agr = 1L:length(tb), position = ordmat[1, order(ordmat[2, ])],
                            counts)
                    })
                    ))
  
  agr <- comp_tab[comp_tab[, "position"] == 1, "agr"]
  bgr <- comp_tab[comp_tab[, "position"] == 1, "bgr"]
  
  for(gr_id in 1L:length(tb)) {
    tb_group <- tb[[bgr[[gr_id]]]]
    ida_to <- agr[gr_id]
    ta_ids <- sort(agr[-gr_id])
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
  
  #list(ed, ta, tb)
  ed
}




