#' Calculate encoding distance
#' 
#' Computes the encoding distance between two encodings.
#' @param a encoding (see \code{\link{validate_encoding}} for more information 
#' about the required structure of encoding).
#' @param b encoding to which \code{a} should be compared. Must have equal number 
#' of groups or less than \code{a}. Both \code{a} and {b} must have the the same 
#' number of elements.
#' @param prop \code{matrix} of physicochemical properties to normalize the 
#' encoding distance.  Each column should 
#' represent properties of the single amino acid/nucleotide. If \code{NULL},
#' encoding distance is not normalized.
#' @param measure \code{character} vector of length one specifying the measure. 
#' Currently avaible measures are \code{"pi"} (partition index) and 
#' \code{"si"} (similarity index).
#' If the parameter \code{prop} is supplied, the encoding distance is normalized by the 
#' factor equal to the sum of distances for each group in \code{a} and the closest group 
#' in \code{b}. The position of a group is defined as the mean value of properties of 
#' amino acids or nucleotides belonging the group.
#' 
#' See the package vignette for more details.
#' @seealso 
#' \code{\link{calc_si}}: compute the similarity index of two encodings.
#' \code{\link{encoding2df}}: converts an encoding to a data frame.
#' \code{\link{validate_encoding}}: validate a structure of an encoding.
#' @return an encoding distance.
#' @export
#' @examples
#' # calculate encoding distance between two encodings of amino acids
#' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#'            `2` = c("k", "h"), 
#'            `3` = c("d", "e"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' 
#' aa2 = list(`1` = c("g", "a", "p", "v", "m", "l", "q"), 
#'            `2` = c("k", "h", "d", "e", "i"), 
#'            `3` = c("f", "r", "w", "y", "s", "t", "c", "n"))
#' calc_ed(aa1, aa2, measure = "pi") 
#'     
#' # the encoding distance between two identical encodings is 0
#' calc_ed(aa1, aa1, measure = "pi") 

calc_ed <- function(a, b, prop = NULL, measure) {
  ed_function <- switch(measure, 
                        si = calc_si_hidden,
                        pi = calc_pi_hidden)
  
  if(!is.null(prop)) {
    if(!is.matrix(prop))
      stop("'prop' must have 'matrix' class.")
    
    if(!(ncol(prop) %in% c(4, 20)))
      stop("'prop' must have 4 or 20 columns.")
  }
  
  # compare temporary a to temporary b
  if(length(a) < length(b)) {
    warning("'a' must be longer than 'b'. Reverting a and b.") 
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
  
  if(!validate_encoding(ta, unlist(tb)))
    stop("Encodings ('a' and 'b') must contain the same elements.")
  
  ed <- ed_function(ta, tb)
  
  if(!is.null(prop)) {
    coords_a <- lapply(a, function(single_subgroup) rowMeans(prop[, single_subgroup, drop = FALSE]))
    coords_b <- lapply(b, function(single_subgroup) rowMeans(prop[, single_subgroup, drop = FALSE]))
    
    norm_factor <- sum(sapply(coords_a, function(single_coords_a) {
      distances <- sapply(coords_b, function(single_coords_b) 
        # vector of distances between groups
        sqrt(sum((single_coords_a - single_coords_b)^2))
      )
      # c(dist = min(distances), id = unname(which.min(distances)))
      min(distances)
    }))
    
    ed <- ed * norm_factor
  }
  
  unname(ed)
}


#' Calculate partition index
#' 
#' Computes the encoding distance between two encodings.
#' @param a encoding (see \code{\link{validate_encoding}} for more information 
#' about the required structure of encoding).
#' @param b encoding to which \code{a} should be compared. Must have equal number 
#' of groups or less than \code{a}. Both \code{a} and {b} must have the the same 
#' number of elements.
#' @details The encoding distance between \code{a} and \code{b} is defined as the 
#' minimum number of amino acids that have to be moved between subgroups of encoding 
#' to make \code{a} identical to \code{b} (order of subgroups in the encoding and amino 
#' acids in a group is unimportant).
#' 
#' If the parameter \code{prop} is supplied, the encoding distance is normalized by the 
#' factor equal to the sum of distances for each group in \code{a} and the closest group 
#' in \code{b}. The position of a group is defined as the mean value of properties of 
#' amino acids or nucleotides belonging the group.
#' 
#' See the package vignette for more details.
#' @seealso 
#' \code{\link{calc_si}}: compute the similarity index of two encodings.
#' \code{\link{encoding2df}}: converts an encoding to a data frame.
#' \code{\link{validate_encoding}}: validate a structure of an encoding.
#' @return an encoding distance.
#' @importFrom partitions listParts
#' @importFrom combinat permn
#' @export
#' @examples
#' # calculate encoding distance between two encodings of amino acids
#' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#'            `2` = c("k", "h"), 
#'            `3` = c("d", "e"), 
#'            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' 
#' aa2 = list(`1` = c("g", "a", "p", "v", "m", "l", "q"), 
#'            `2` = c("k", "h", "d", "e", "i"), 
#'            `3` = c("f", "r", "w", "y", "s", "t", "c", "n"))
#' calc_pi(aa1, aa2) 
#'     
#' # the encoding distance between two identical encodings is 0
#' calc_pi(aa1, aa1) 
#'  

calc_pi <- function(a, b) {
  calc_ed(a, b, measure = "pi")
}

calc_pi_hidden <- function(ta, tb) {
    
  # encoding distance - distance between encodings
  ed <- 0
  
  # exclude identical subgroups in a and b, because hey do not affect ed.
  ident_gr <- which(sapply(ta, function(single_subgroup_a)
    sapply(tb, function(single_subgroup_b)
      identical(sort(single_subgroup_a), sort(single_subgroup_b)))), arr.ind = TRUE)
  
  if(length(ident_gr) != 0) {
    ta <- ta[-ident_gr[, "col"]]
    tb <- tb[-ident_gr[, "row"]]
  }
  
  # if encodings are identical, ta and tb are 0
  if(length(ta) == 0 && length(tb) == 0)
    return(ed)
  
  len_a <- length(ta)
  len_b <- length(tb)
  
  ed <- if(len_a != len_b) {
    if(len_b == 1) {
      if(all(lengths(ta) == max(lengths(ta)))) {
        sum(lengths(ta)[-1])
      } else {
        sum(lengths(ta)[lengths(ta) != max(lengths(ta))])
      }
    } else {
      
      merges <- create_merges(len_a, len_b)
      
      # list of merged
      reduced_size <- lapply(merges, function(single_merge) {
        ed <- sum(sapply(single_merge, function(i) {
          all_lengths <- lengths(ta[i])
          sum(all_lengths[-which.max(all_lengths)])
        }))
        
        enc <- lapply(single_merge, function(i) {
          unlist(ta[i], use.names = FALSE)
        })
        
        list(enc = enc,
             ed = ed)
      })
      
      all_eds <- sapply(reduced_size, function(single_ta) 
        calc_ed_single(single_ta[["enc"]], tb, single_ta[["ed"]]))
      
      min(all_eds)
      
    }
  } else {
    calc_ed_single(ta, tb)
  }
  
  ed
}

calc_ed_single <- function(ta, tb, ed = 0) {
  
  len_b <- length(tb)
  # all permutations of assigning subgroups from encoding 'a' to subgroups of 'b'
  perms <- expand.grid(lapply(1L:len_b, function(dummy) 1L:len_b))
  perms <- perms[apply(perms, 1, function(single_permutation) 
    length(unique(single_permutation))) == len_b, ]
  bgr <- 1L:len_b
  
  # compute the ed for all permutations and get the minimum
  ed + min(apply(perms, 1, function(agr) {
    sum(lengths(lapply(1L:len_b, function(gr_id) {
      tb_group <- tb[[bgr[[gr_id]]]]
      # added amino acids from other groups
      unlist(lapply(agr[-gr_id], function(ida_from) {
        # id of the single amino acid
        ta[[ida_from]][which(ta[[ida_from]] %in% tb_group)]
      }))
    })))
  }))
}

create_merges <- function(x, n_groups) {
  all_merges <- listParts(x)
  # all partitions of x of the proper length
  chosen_merges <- all_merges[lengths(all_merges) == n_groups]
  
  # all permutations of a single partition
  perms <- permn(n_groups)
  
  unlist(lapply(chosen_merges, function(single_merge) {
    lapply(perms, function(single_perm) 
      unname(single_merge)[single_perm])
  }), recursive = FALSE)
}
