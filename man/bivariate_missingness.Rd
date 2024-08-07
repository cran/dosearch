% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dosearch-package.R
\docType{data}
\name{bivariate_missingness}
\alias{bivariate_missingness}
\title{Systematic Analysis of Bivariate Missing Data Problems}
\format{
A data frame with 6144 rows and 8 variables:
\describe{
\item{graph}{the graph of the instance.}
\item{nedges}{number of edges in the graph (directed and bidirected).}
\item{arrowXtoY}{whether the graph contains an arrow
from \eqn{X} to \eqn{Y} or not.}
\item{jointXY}{identifiability of the joint distribution
of \eqn{X} and \eqn{Y}}
\item{marginX}{identifiability of the marginal distribution of \eqn{X}.}
\item{marginY}{identifiability of the marginal distribution of \eqn{Y}.}
\item{YcondX}{identifiability of the conditional distribution of \eqn{Y}
given \eqn{X}.}
\item{YdoX}{identifiability of the causal effect of \eqn{X} on \eqn{Y}.}
}
}
\source{
Tikka et al. \url{https://arxiv.org/abs/1902.01073}
}
\usage{
data(bivariate_missingness)
}
\description{
This data set contains the results of a systematic analysis of all missing
data problems of two variables. Each problem is associated with a graph
containing two vertices, \eqn{X} and \eqn{Y}, and their response indicators,
\eqn{R_X} and \eqn{R_Y}.
}
\keyword{datasets}
