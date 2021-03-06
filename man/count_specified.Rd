% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_specified.R
\name{count_specified}
\alias{count_specified}
\title{Count specified n-grams}
\usage{
count_specified(seq, ngrams)
}
\arguments{
\item{seq}{vector or matrix describing sequence(s).}

\item{ngrams}{vector of n-grams.}
}
\value{
A \code{\link[slam]{simple_triplet_matrix}} where columns represent
n-grams and rows sequences.
}
\description{
Counts specified n-grams in the input sequence(s).
}
\details{
\code{\link{count_specified}} counts only selected n-grams declared by
user in the \code{ngrams} parameter. Declared n-grams must be written using the
\code{biogram} notation.
}
\examples{
seqs <- matrix(c(1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 3, 4, 1, 2, 2, 4), nrow = 2)
count_specified(seqs, ngrams = c("1.1.1_0.0", "2.2.2_0.0", "1.1.2_0.0"))

seqs <- matrix(sample(1L:5, 200, replace = TRUE), nrow = 20)
count_specified(seqs, ngrams = c("2_4.2_0", "2_1.4_0", "3_1.3_0",
                                 "2_4.2_1", "2_1.4_1", "3_1.3_1",
                                 "2_4.2_2", "2_1.4_2", "3_1.3_2"))
}
\seealso{
Count all possible n-grams: \code{\link{count_ngrams}}.
}
