% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{methods-unitquantreg}
\alias{methods-unitquantreg}
\alias{print.unitquantreg}
\alias{summary.unitquantreg}
\alias{coef.unitquantreg}
\alias{vcov.unitquantreg}
\alias{logLik.unitquantreg}
\alias{confint.unitquantreg}
\alias{fitted.unitquantreg}
\alias{terms.unitquantreg}
\alias{model.frame.unitquantreg}
\alias{model.matrix.unitquantreg}
\alias{update.unitquantreg}
\alias{print.unitquantregs}
\alias{summary.unitquantregs}
\title{Methods for \code{unitquantreg} and \code{unitquantregs} objects}
\usage{
\method{print}{unitquantreg}(x, digits = max(4, getOption("digits") - 3), ...)

\method{summary}{unitquantreg}(object, correlation = FALSE, ...)

\method{coef}{unitquantreg}(object, type = c("full", "quantile", "shape"), ...)

\method{vcov}{unitquantreg}(object, ...)

\method{logLik}{unitquantreg}(object, ...)

\method{confint}{unitquantreg}(object, parm, level = 0.95, ...)

\method{fitted}{unitquantreg}(object, type = c("quantile", "shape", "full"), ...)

\method{terms}{unitquantreg}(x, type = c("quantile", "shape"), ...)

\method{model.frame}{unitquantreg}(formula, ...)

\method{model.matrix}{unitquantreg}(object, type = c("quantile", "shape"), ...)

\method{update}{unitquantreg}(object, formula., ..., evaluate = TRUE)

\method{print}{unitquantregs}(x, digits = max(3, getOption("digits") - 3), ...)

\method{summary}{unitquantregs}(object, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{digits}{minimal number of \emph{significant} digits.}

\item{...}{additional argument(s) for methods. Currently not used.}

\item{object, x}{fitted model object of class \code{\link{unitquantreg}}.}

\item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed. Default is \code{FALSE}.}

\item{type}{character indicating type of fitted values to return.}

\item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{formula}{an R formula.}

\item{formula.}{Changes to the formula see \code{\link{update.formula}} for details.}

\item{evaluate}{If true evaluate the new call else return the call.}
}
\description{
Methods for extracting information from fitted unit quantile regression
objects of class \code{unitquantreg} and \code{unitquantregs}.
}
\author{
André F. B. Menezes
}
