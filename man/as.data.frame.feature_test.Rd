% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_test_class.R
\name{as.data.frame.feature_test}
\alias{as.data.frame.feature_test}
\title{Coerce feature_test object to a data frame}
\usage{
\method{as.data.frame}{feature_test}(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)
}
\arguments{
\item{x}{object of class \code{\link{feature_test}}.}

\item{row.names}{ignored.}

\item{optional}{ignored.}

\item{stringsAsFactors}{logical: should the character vector be converted 
to a factor?.}

\item{...}{additional arguments to be passed to or from methods.}
}
\value{
a data frame with four columns: names of n-gram, p-values,
occurrences in positive and negative sequences.
}
\description{
Coerce results of \code{\link{test_features}} function to a 
\code{\link[base]{data.frame}}.
}
