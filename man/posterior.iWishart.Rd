% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wishart.R
\name{posterior.iWishart}
\alias{posterior.iWishart}
\title{Posterior based on an inverse Wishart distribution for covariance matrix}
\usage{
posterior.iWishart(n, S, m0, V0)
}
\arguments{
\item{n}{Sample size.}

\item{S}{Sample covariance matrix (MLE).}

\item{m0}{Degrees of freedom of the inverse Wishart prior.}

\item{V0}{Scale matrix of the inverse Wishart prior.}
}
\value{
The mean and variance of the inverse Wishart distribution.
}
\description{
Combine data and inverse Wishart prior of covariance matrix.
}
\seealso{
\code{\link{posterior.Wishart}} precision matrix.
}
