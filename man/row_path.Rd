\name{row_path}
\alias{row_path}
\title{ROW estimates along a grid of lambda values}
\description{Computes ROW estimates and sandwich standard errors for a grid of \eqn{\lambda/K^2}
values, from close to outcome regression (small \eqn{\lambda}) to AIPW (\eqn{\lambda = \infty}).
A stable path indicates agreement between the outcome model and the weights.}
\usage{
row_path(data, treatment, outcome, ps_formula, outcome_formula = ps_formula,
         grid = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf), reference = NULL)
}
\arguments{
  \item{data, treatment, outcome, ps_formula, outcome_formula, reference}{as in \code{\link{row_ate}}.}
  \item{grid}{values of \eqn{\lambda/K^2}.}
}
\value{A data frame with columns \code{lambda_K2}, \code{contrast}, \code{estimate}, \code{se}.}
\examples{
data(ncds)
row_path(ncds, "Dmult", "wage", ~ qmab + qvab + paed_u + maed_u, reference = "None",
         grid = c(0.5, 2, 8, Inf))
}
