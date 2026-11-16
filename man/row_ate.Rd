\name{row_ate}
\alias{row_ate}
\alias{print.row}
\alias{coef.row}
\alias{confint.row}
\title{Regularised outcome-weighted (ROW) estimate of the average treatment effect}
\description{
Estimates pairwise average treatment effects (each treatment versus a reference) for two or
more treatments. Per-arm outcome models are averaged over all patients and corrected with
residuals weighted by \eqn{w = 1/(e_A(X) + 1/\lambda)}, where \eqn{e_A(X)} is the generalised
propensity score of the treatment received. By default \eqn{\lambda} is chosen by an
outcome-blind balance rule: the smallest \eqn{\lambda} (on the \eqn{\lambda/K^2} scale, up to
\code{cap}) for which every weighted arm matches the full population on the covariates and the
squares of non-binary covariates (maximum standardised mean difference \eqn{\le}
\code{smd_threshold}). Standard errors come from an M-estimation sandwich.
}
\usage{
row_ate(data, treatment, outcome, ps_formula, outcome_formula = ps_formula,
        lambda = "balance", cap = 8, smd_threshold = 0.1,
        grid = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128),
        reference = NULL, level = 0.95)
\method{print}{row}(x, digits = 3, ...)
\method{coef}{row}(object, ...)
\method{confint}{row}(object, parm, level = object$level, ...)
}
\arguments{
  \item{data}{a data frame.}
  \item{treatment}{name of the treatment variable (two or more levels).}
  \item{outcome}{name of the continuous outcome variable.}
  \item{ps_formula}{one-sided formula for the multinomial generalised propensity score model, e.g. \code{~ x1 + x2}.}
  \item{outcome_formula}{one-sided formula for the per-arm linear outcome models (default: \code{ps_formula}).}
  \item{lambda}{\code{"balance"} (default) or a positive number on the \eqn{\lambda/K^2} scale; \code{Inf} gives AIPW.}
  \item{cap}{largest \eqn{\lambda/K^2} considered by the balance rule (default 8).}
  \item{smd_threshold}{balance target for the rule (default 0.1).}
  \item{grid}{candidate values of \eqn{\lambda/K^2} for the balance rule.}
  \item{reference}{reference treatment level (default: first level).}
  \item{level}{confidence level.}
  \item{x, object}{a \code{"row"} object.}
  \item{digits}{digits to print.}
  \item{parm}{ignored.}
  \item{...}{ignored.}
}
\value{An object of class \code{"row"}: a list with \code{estimates} (contrast, estimate, sandwich SE,
confidence limits, naive SE), arm \code{means}, the selected \code{lambda_K2} and \code{lambda}, the
\code{max_smd} of the weighted arms, the \code{balance_path}, the \code{weights} and the fitted \code{ps}.}
\references{
Karletsos D. When overlap weights miss the target: a regularised outcome-weighted estimator (ROW) for
the average treatment effect in multi-arm real-world comparisons. PHUSE EU Connect 2026, paper AS03.
}
\examples{
data(ncds)
f <- ~ white + maemp + factor(scht) + qmab + qmab2 + qvab + qvab2 +
       paed_u + maed_u + agepa + agema + sib_u
fit <- row_ate(ncds, treatment = "Dmult", outcome = "wage", ps_formula = f, reference = "None")
fit
confint(fit)
}
