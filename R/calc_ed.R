#' Calculate encoding distance
#' 
#' Computes the encoding distance between two encodings.
#' @param a encoding (see \code{\link{validate_encoding}} for more information 
#' about the required structure of encoding).
#' @param b encoding to which \code{a} should be compared. Must have equal number 
#' of groups or less than \code{a}.
#' @param prop \code{matrix} of physicochemical properties to normalize the 
#' encoding distance.  Must have either four or 20 columns and each column should 
#' represent properties of the single amino acid/nucleotide. If \code{NULL},
#' encoding distance is not normalized.
#' @details The encoding distance between \code{a} and \code{b} is defined as the 
#' minimum number of amino acids that have to be moved between subgroups of encoding 
#' to make \code{a} identical to \code{b} (order of subgroups in the encoding and amino 
#' acids in a group is unimportant).
#' 
#' If the parameter \code{prop} is supplied, the encoding distance is normalized by the 
#' factor equal to the sum of distances for each group in \code{a} and the closest group 
#' in \code{b}. The position of a group is defined as the mean value of properties of 
#' amino acids or nucleotides belonging the group.
#' @seealso \code{\link{validate_encoding}}.
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
#' calc_ed(aa1, aa2) 
#'     
#' # the encoding distance between two identical encodings is 0
#' calc_ed(aa1, aa1) 
#'  

calc_ed <- function(a, b, prop = NULL) {
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
  
  if(!is.null(prop)) {
    if(!is.matrix(prop))
      stop("'prop' must have 'matrix' class.")
    
    if(!(ncol(prop) %in% c(4, 20)))
      stop("'prop' must have 4 or 20 columns.")
  }
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
      
#      # indices of groups to merge
#      merges <- unlist(lapply(1L:(len_a - len_b + 1), function(n_subgroups) {
#        combn(1L:len_a, n_subgroups, simplify = FALSE)
#      }), recursive = FALSE)
#      
#      # list of merged
#      reduced_size <- lapply(merges, function(single_merge) {
#        
#        all_lengths <- list(lengths(ta[single_merge], use.names = FALSE),
#                            lengths(ta[-single_merge], use.names = FALSE))
#        
#        list(enc = list(unlist(ta[single_merge], use.names = FALSE), 
#                        unlist(ta[-single_merge], use.names = FALSE)), # groups after merge
#             ed = sum(sapply(all_lengths, function(single_length) {
#               sum(single_length[-which.max(single_length)])
#             })))
#      })
      
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

calc_ed_single <- function(ta, tb, ed = 0) {
  
  # comparision table, for the so called 'future' when I have enough time to write something 
  # smarter
  #   comp_tab <- cbind(bgr = sort(rep(1L:length(ta), length(tb))),
  #                     do.call(rbind, lapply(tb, function(single_subgroup_b) {
  #                       counts <- do.call(rbind, lapply(ta, function(single_subgroup_a)
  #                         data.frame(counts = sum(single_subgroup_a %in% single_subgroup_b))))
  #                       
  #                       # order matrix
  #                       # first row - sequence of groups
  #                       # second row - groups' order
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


#' Validate encoding
#'
#' Checks if an encoding has the proper structure 
#' @param x encoding.
#' @param u \code{integer}, \code{numeric} or \code{character} vector of all
#' elements belonging to the encoding. See Details.
#' @keywords manip
#' @return \code{TRUE} if the \code{x} is a correctly reduced \code{u}, 
#' \code{FALSE} in any other cases.
#' @details The encoding is a reduced alphabet, the list of groups to which 
#' elements of sequence should be aggregated. All elements of the alphabet (all 
#' amino acids or all nucleotides) should appear in the encoding.
#' @export
#' @keywords manip
#' @seealso 
#' \code{\link{calc_ed}}: calculate the encoding distance between two encodings.
#' @examples
#' enc1 = list(`1` = c("a", "t"), 
#'             `2` = c("g", "c"))
#' # see if enc1 is the correctly reduced nucleotide (DNA) alphabet
#' validate_encoding(enc1, c("a", "c", "g", "t"))
#' 
#' # enc1 is not the RNA alphabet, so the results is FALSE
#' validate_encoding(enc1, c("a", "c", "g", "u"))
#' 
#' # validate_encoding works also on other notations
#' enc2 = list(a = c(1, 4),
#'             b = c(2, 3))
#' validate_encoding(enc2, 1L:4)

validate_encoding <- function(x, u) {
  if(!is.list(x))
    stop("'x' must have 'list' class.")
  all(sort(unlist(x)) == sort(u))
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

