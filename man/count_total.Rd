% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices_and_positions.R
\name{count_total}
\alias{count_total}
\title{Count total number of n-grams}
\usage{
count_total(seq, n, d)
}
\arguments{
\item{seq}{a vector or matrix describing sequence(s).}

\item{n}{\code{integer} size of n-gram.}

\item{d}{\code{integer} vector of distances between elements of n-gram (0 means 
consecutive elements). See Details.}
}
\value{
An \code{integer} rperesenting the total number of n-grams.
}
\description{
Computes total number of n-grams that can be extracted from sequences.
}
\details{
The maximum number of possible n-grams is limited by their length 
and the distance between elements of the n-gram.
}
\note{
A format of \code{d} vector is discussed in Details of 
\code{\link{count_ngrams}}. The maximum
}
\examples{
seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
# make several sequences shorter by replacing them partially with NA
seqs[8L:11, 46L:50] <- NA
seqs[1L, 31L:50] <- NA
count_total(seqs, 3, c(1, 0))
}
