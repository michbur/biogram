#' Generate single unigram 
#'
#' Assign randomly generated properties to a single unigram.
#' @param unigram_ranges list of ranges containing respective properties. If named, 
#' names are preserved.
#' @seealso \code{generate_single_unigram} is a helper function for 
#' \code{\link{generate_unigrams}}.
#' @export
#' @examples 
#' generate_single_unigram(list(P1 = c(0, 0.5), 
#'                              P2 = c(0.2, 0.4),
#'                              P3 = c(0.5, 1),
#'                              P4 = c(0, 0)))
#' 
generate_single_unigram <- function(unigram_ranges) {
  unigram_props <- lapply(unigram_ranges, function(i) {
    runif(1, min = i[1], max = i[2])
  })
  names(unigram_props) <- names(unigram_ranges)
  unlist(unigram_props)
}


#' Generate unigrams
#'
#' Generates an alphabet of unigrams based on given list of properties.
#' @param unigram_list a list of unigrams' parameters. See Details.
#' @param unigram_names names of unigrams. If not \code{NULL}, will
#' overwrite any existing unigram names.
#' @param prop_names names of properties. If not \code{NULL}, will 
#' overwrite any existing names.
#' @details Unigram parameters are represented as a list of intervals, where 
#' each interval corresponds to a different property. The function generate 
#' unigrams randomly choosing values of properties from given intervals 
#' using uniform distribution. All lists of ranges should have the same length, 
#' which equils to describing each unigram using the same properties. 
#' @export
#' @examples 
#' props1 <- list(P1 = c(0, 0.5), 
#'                P2 = c(0.2, 0.4),
#'                P3 = c(0.5, 1),
#'                P4 = c(0, 0))
#' 
#' props2 <- list(P1 = c(0.5, 1), 
#'                P2 = c(0.4, 1),
#'                P3 = c(0, 0.5),
#'                P4 = c(1, 1))
#' 
#' 
#' alph <- generate_unigrams(c(replicate(8, props1, simplify = FALSE),
#'                           replicate(12, props2, simplify = FALSE)),
#'                           unigram_names = letters[1L:20])
#'                           
generate_unigrams <- function(unigram_list,
                              unigram_names = NULL,
                              prop_names = NULL) {
  if(length(unique(lengths(unigram_list))) != 1)
    stop("All unigrams must be defined by the same number of properties (each element of \
         unigram_list must have the same length).")
  
  unigrams <- do.call(cbind, lapply(unigram_list, function(single_unigram) {
    generate_single_unigram(single_unigram)
  }))
  
  if(!is.null(prop_names)) 
    rownames(unigrams) <- prop_names
  
  if(!is.null(unigram_names)) 
    colnames(unigrams) <- unigram_names
  
  unigrams
}


#' Generate single region 
#'
#' Generate a region using an alphabet of unigrams and considering provided 
#' set of rules.
#' @param alphabet the unigram alphabet. Columns are equivalent to unigrams 
#' and rows to particular properties.
#' @param reg_len the number of unigrams inside the region.
#' @param prop_ranges required intervals of properties of unigrams in the region. 
#' See Details.
#' @param exactness a \code{numeric} value between 0 and 1 defining how stricly 
#' unigrams are kept within \code{prop_ranges}. If 1, only unigrams within 
#' \code{prop_ranges} are inside the region. if 0.9, there is 10% chance that 
#' unigrams that are not in the \code{prop_ranges} will be inside the region.
#' @export
#' @examples 
#' props1 <- list(P1 = c(0, 0.5), 
#'                P2 = c(0.2, 0.4),
#'                P3 = c(0.5, 1),
#'                P4 = c(0, 0))
#' 
#' props2 <- list(P1 = c(0.5, 1), 
#'                P2 = c(0.4, 1),
#'                P3 = c(0, 0.5),
#'                P4 = c(1, 1))
#' 
#' 
#' alph <- generate_unigrams(c(replicate(8, props1, simplify = FALSE),
#'                             replicate(12, props2, simplify = FALSE)),
#'                           unigram_names = letters[1L:20])
#' 
#' rules1 <- list(P1 = c(0.5, 1), 
#'                P2 = c(0.4, 1),
#'                P3 = c(0, 0.5),
#'                P4 = c(1, 1))
#' 
#' generate_single_region(alph, 10, rules1, 0.9)
generate_single_region <- function(alphabet, reg_len, prop_ranges, exactness) {
  
  min_range <- sapply(prop_ranges, function(i) i[1])
  max_range <- sapply(prop_ranges, function(i) i[2])
  
  regionality <- apply(alphabet >= min_range & alphabet <= max_range, 2, all)
  regional_unigrams <- names(which(regionality))
  nonregional_unigrams <- names(which(!regionality))
  
  rigorous_id <- runif(reg_len)
  region <- rep(NA, reg_len)
  region[rigorous_id <= exactness] <- sample(regional_unigrams, 
                                         size = sum(rigorous_id <= exactness), 
                                         replace = TRUE)
  
  region[rigorous_id > exactness] <- sample(nonregional_unigrams, 
                                        size = sum(rigorous_id > exactness), 
                                        replace = TRUE)
  
  region
}


#' Generate sequence
#'
#' Generate a sequences using an alphabet of unigrams and set of rules.
#' @param alphabet the unigram alphabet. Columns are equivalent to unigrams 
#' and rows to particular properties.
#' @param regions a list of rules describing regions.
generate_sequence <- function(alphabet, regions) {
  reg_list <- lapply(regions, function(single_region) {
    generate_single_region(alphabet, 
                           reg_len = single_region[["len"]],
                           prop_ranges = single_region[["prop"]],
                           exactness = single_region[["exactness"]])
  })
}

