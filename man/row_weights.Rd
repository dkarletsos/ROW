\name{row_weights}
\alias{row_weights}
\title{ROW weights}
\description{Regularised inverse probability weights \eqn{w = 1/(e_A + 1/\lambda)}.}
\usage{row_weights(ps, A, lambda)}
\arguments{
  \item{ps}{matrix of generalised propensity scores (one column per treatment level).}
  \item{A}{integer vector of treatment levels received (column index of \code{ps}).}
  \item{lambda}{regularisation parameter on the absolute scale (\code{Inf} gives \eqn{1/e_A}).}
}
\value{A numeric vector of weights.}
\examples{
ps <- matrix(c(0.7, 0.2, 0.1, 0.1, 0.3, 0.6), 2, byrow = TRUE)
row_weights(ps, A = c(3, 1), lambda = 18)
}
