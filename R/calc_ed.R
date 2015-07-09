#' Calculate encoding distance
#' 
#' Compares two encodings and calculates encoding distance between them.
#' @param a list of groups to which elements of sequence should be aggregated.
#' @param b list of groups to which \code{a} should be compared.
#' @return a list of three:
#' #' \itemize{
#' \item{ed: encoding distance.}
#' \item{a: \code{a} parameter after comparing to \code{b}.}
#' \item{b: \code{b} parameter.}
#' }
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

calc_ed <- function(a, b) {
  #I'm comparing a to b
  #temporary b and temporary a
  ta <- a
  tb <- b
  
  if(!all(sort(unlist(a)) == sort(unlist(b))))
    stop("Encodings must contain the same elements.")
  
  if(length(a) < length(b)) {
    warning("a cannot be shorter than b. Reverting a and b.") 
    tb <- a
    ta <- b
  }
  
  names(ta) <- NULL
  names(tb) <- NULL
  uneven <- length(ta) - length(tb)
  
  #encoding distance - distance between encodings
  ed <- 0
  
  #rows b
  #columns a
  comp_tab <- sapply(ta, function(single_subgroup_a)
    sapply(tb, function(single_subgroup_b)
      sum(single_subgroup_a %in% single_subgroup_b)/length(single_subgroup_a)))
  #loop for moving aa between groups
  #TRUE if groups are not similiar
  #diff <- comp_tab != 0 & comp_tab != 1

  while(any(!(comp_tab %in% c(0, 1)))) {
    #change 1 to -1 to use which
    comp_tab[comp_tab == 1] <- -1
    
    #cannot use which.max, it doesn't use arr.ind
    #choose groups to switch aa
    gr_id <- which(comp_tab == max(comp_tab), arr.ind = TRUE)
    gr_id <- order(comp_tab[, which(comp_tab == max(comp_tab), arr.ind = TRUE)[, "col"]], decreasing = TRUE)[c(1, 2)]
    arrb <- which(comp_tab == max(comp_tab), arr.ind = TRUE)
    arra <- which(comp_tab == max(comp_tab), arr.ind = TRUE)
    #the biggest in column idb, second largest is ida
    arra[, "row"] <- order(comp_tab[, arrb[, "col"]], decreasing = TRUE)[2]
    idb <- arrb[, "col"]
    ida <- arra[, "row"]
    
    #id of the single amino acid
    el_id <- which(ta[[idb]] %in% tb[[ida]])
    #add amino acid to second group
    ta[[ida]] <- c(ta[[ida]], ta[[idb]][el_id])
    #remove amino acid from the first group
    ta[[idb]] <-  ta[[idb]][-el_id]
    ed <- ed + 1
    
    #remove empty sets
    if(any(lengths(ta) == 0)) {
      uneven <- uneven - 1
      ta <- ta[!lengths(ta) == 0]
    }
    
    comp_tab <- sapply(ta, function(single_subgroup_a)
      sapply(tb, function(single_subgroup_b)
        sum(single_subgroup_a %in% single_subgroup_b)/length(single_subgroup_a)))
  }
  
  if(any(rowSums(comp_tab) != 1)) {
    #concatenate smaller groups into one big
    #check which row in comp_tab has two 1s 
    conc_group <- unname(do.call(c, ta[as.logical(comp_tab[rowSums(comp_tab) != 1, ])]))
    #choose smaller group as a group that is moved
    ed <- ed + min(lengths(ta[as.logical(comp_tab[rowSums(comp_tab) != 1, ])]))
    ta <- c(ta[!as.logical(comp_tab[rowSums(comp_tab) != 1, ])],
            list(conc_group))
  }
  
  list(ed = ed,
       a = sapply(ta, sort),
       b = sapply(tb, sort))
}

