# Check chosen criterion
# 
# Checks if the criterion is viable or matches it to the list of implemented 
# criterions.
# 
# @param input_criterion a \code{character} string, criterion from input.
# @param criterion_names list of implemented \code{\link{criterions}}, always in lowercase.
# @export
# @return a list of three: 
# \itemize{
# \item{criterion name,}
# \item{its function,}
# \item{nice name for outputs.}
# }
# @seealso
# All implemented \code{\link{criterions}}.

check_criterion <- function(input_criterion, criterion_names = c("ig")) {
  #think twice about grep
  valid_name <- criterion_names[grepl(tolower(input_criterion), criterion_names)]
  
  if (length(valid_name) == 0)
    stop("Name ", input_criterion, " cannot be associated with any available criterion.")
  
  if (length(valid_name) > 1)
    stop("Name ", input_criterion, " is too ambiguous. Rerun with more precise name.")
  
  criterion_data <- switch(valid_name,
                               ig = list(crit_function = calc_ig, nice_name = "Information Gain"))
  #TO DO - should also return the full name of criterion for purpose of summaries/plots
  c(crit_name = valid_name, criterion_data)
}
