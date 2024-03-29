% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dv.R
\name{dv}
\alias{dv}
\title{dv}
\usage{
dv(RT, Resp, color = "lightgray")
}
\arguments{
\item{RT}{A required \eqn{N × J} \code{matrix} containing the response
times of \eqn{N} individuals to \eqn{J} items. Zero and missing values need
to be coded as \code{NA}.}

\item{Resp}{A required \eqn{N × J} \code{matrix} or \code{data.frame}
containing the item responses of \eqn{N} individuals to \eqn{J} items. The
correct responses are marked by 1 and incorrect responses by 0. Missing values
need to be coded as \code{NA}.}

\item{color}{A non-required \code{character string} concerning the color to be used to fill the bars.}
}
\value{
a 3*4 matrix of plots containing the RT distribution of each item as well as the proportion of correctness.
}
\description{
dv gives histogram of RT distribution as well as the proportion of correctness
}
\examples{
x <- as.data.frame(matrix(c(1,2,3,4,5,6), nrow = 2))
dv(Resp,RT, color="grey")
}
