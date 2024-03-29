% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prim.R
\name{prim}
\alias{prim}
\title{prim}
\usage{
prim(ini)
}
\arguments{
\item{ini}{an object of class \code{ini} containing the initial values for the log-normal RT model}
}
\value{
an object of class \code{SecMin}. The method for \code{SecMin} objects
is \code{\link{prepftr}} for preparing the features to fit unsupervised or
supervised machine learning models.
\item{alpha}{estimates of alpha values}
\item{beta}{estimates of beta values}
\item{speed}{estimates of person speed values}
\item{LZ.cheating.cases}{predicted cheating cases based on LZ index}
\item{KLD.cheating.cases}{predicted cheating cases based on KLD index}
\item{KLD.value}{person fit estimates depend on RT-based KLD measure}
\item{LZ.value}{person fit estimates depend on RT-based LZ index}
\item{Ht.value}{person fit estimates depend on response-based Ht index}
\item{Lzstar.value}{person fit estimates depend on IRT-based Lzstar index}
\item{NCI.value}{person fit estimates depend on response-based NCI index}
\item{general}{results of the lognormal model}
}
\description{
prim fits the log-normal RT model and return parameter estimates as well as other
results of the model.
}
