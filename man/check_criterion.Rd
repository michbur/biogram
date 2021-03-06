% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criterions.R
\name{check_criterion}
\alias{check_criterion}
\title{Check chosen criterion}
\usage{
check_criterion(input_criterion, criterion_names = c("ig", "kl", "cs"))
}
\arguments{
\item{input_criterion}{\code{character} string, criterion from input.}

\item{criterion_names}{list of implemented criterions, always in lowercase.}
}
\value{
a list of three:
\itemize{
\item{criterion name,}
\item{its function,}
\item{nice name for outputs.}
}
}
\description{
Checks if the criterion is viable or matches it to the list of implemented
criterions.
}
\seealso{
Calculate the value of criterion: \code{\link{calc_criterion}}.
}
