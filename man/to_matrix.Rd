% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.R
\name{to_matrix}
\alias{to_matrix}
\alias{to_matrix.data.frame}
\title{to_matrix}
\usage{
to_matrix(x, rownames = NULL)

\method{to_matrix}{data.frame}(x, rownames = NULL)
}
\arguments{
\item{x}{data.frame to conver.}

\item{rownames}{column to be used as row names.}
}
\description{
Convert a data.frame into a matrix, optionally setting one of the columns
as the row names.
}
