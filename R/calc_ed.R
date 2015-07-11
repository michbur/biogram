#' Calculate encoding distance
#' 
#' Compares two encodings and calculates encoding distance between them.
#' @param a list of groups to which elements of sequence should be aggregated.
#' @param b list of groups to which \code{a} should be compared.
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
  #I'm comparing a to b
  #temporary b and temporary a
  ta <- a
  tb <- b
  
  if(any(lengths(ta) == 0))
    ta <- ta[lengths(ta) != 0]
  
  if(any(lengths(tb) == 0))
    tb <- tb[lengths(tb) != 0]
  
  if(length(unlist(a)) != length(unlist(b)))
    stop("'a' and 'b' must contain the same number of elements.")
  
  if(!all(sort(unlist(a)) == sort(unlist(b))))
    stop("'a' and 'b' must contain the same elements.")
  
  if(length(a) < length(b)) {
    warning("a cannot be shorter than b. Reverting a and b.") 
    tb <- a
    ta <- b
  }
  
  names(ta) <- NULL
  names(tb) <- NULL
  
  #encoding distance - distance between encodings
  ed <- 0
  
  #rows b
  #columns a
  comp_tab <- create_comp_tab(ta, tb)
  
  #loop for moving aa between groups
  #TRUE if groups are not similiar
  #diff <- comp_tab != 0 & comp_tab != 1
  
  while(any(!(comp_tab %in% c(0, 1)))) {
    #change 1 to -1 to use which
    comp_tab[comp_tab == 1] <- -1
    
    #cannot use which.max, it doesn't use arr.ind
    #choose groups to switch aa
    arrb <- which(comp_tab == max(comp_tab), arr.ind = TRUE)
    #case when both groups are 0.5 - use only the first row
    if(nrow(arrb) > 1) {
      arrb <- arrb[1, , drop = FALSE]
      idb <- arrb[, "col"]
      ida <- order(comp_tab[, arrb[, "col"]], decreasing = TRUE)[1]
    } else{
      #the biggest in column idb, second largest is ida
      idb <- arrb[, "col"]
      ida <- order(comp_tab[, arrb[, "col"]], decreasing = TRUE)[2]
    }
    
    #id of the single amino acid
    el_id <- which(ta[[idb]] %in% tb[[ida]])
    #add amino acid to second group
    ta[[ida]] <- c(ta[[ida]], ta[[idb]][el_id])
    #remove amino acid from the first group
    ta[[idb]] <-  ta[[idb]][-el_id]
    ed <- ed + 1

    #remove empty sets
    if(any(lengths(ta) == 0)) {
      ta <- ta[!lengths(ta) == 0]
    }
    
    comp_tab <- create_comp_tab(ta, tb)
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
  list(ed = ed,
       a = sapply(ta, sort),
       b = sapply(tb, sort))
  #ed
}

create_comp_tab <- function(ta, tb) {
  comp_tab <- sapply(ta, function(single_subgroup_a)
    sapply(tb, function(single_subgroup_b)
      sum(single_subgroup_a %in% single_subgroup_b)/length(single_subgroup_a)))
  if(class(comp_tab) != "matrix") {
    comp_tab <- matrix(comp_tab, nrow = length(tb), ncol = length(ta))
  }
  comp_tab
}
