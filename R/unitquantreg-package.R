#' @docType package
#' @name unitquantreg-package
#' @aliases unitquantreg-package
#'
#' @title Overview of the unitquantreg package
#'
#' @description The \pkg{unitquantreg} R package implements the probability density function, quantile function,
#' cumulative distribution function and random number generation function for 11 distributions with support on the unit-interval. Also,
#' \pkg{unitquantreg} also implements a general infrastructure for parametric quantile regression model based on these distributions.
#'
#' @details
#' \code{\link[unitquantreg]{johnsonsb}:} implements the \[d-q-p-r\]johnsonsb function for the Johnson SB distribution.
#'
#' \code{\link[unitquantreg]{kum}:} implements the \[d-q-p-r\]kum function for the Kumaraswamy distribution.
#'
#' \code{\link[unitquantreg]{leeg}:} implements the \[d-q-p-r\]leeg function for the Log-extended exponential-geometric distribution.
#'
#' \code{\link[unitquantreg]{ubs}:} implements the \[d-q-p-r\]ubs function for the unit-Birnbaum-Saunders distribution.
#'
#' \code{\link[unitquantreg]{uburrxii}:} implements the \[d-q-p-r\]uburrxii function for the unit-Burr-XII distribution.
#'
#' \code{\link[unitquantreg]{uchen}:} implements the \[d-q-p-r\]uchen function for the unit-Chen distribution.
#'
#' \code{\link[unitquantreg]{ughne}:} implements the \[d-q-p-r\]ughne function for the unit-Half-Normal-E distribution.
#'
#' \code{\link[unitquantreg]{ughnx}:} implements the \[d-q-p-r\]ughnx function for the unit-Half-Normal-X distribution.
#'
#' \code{\link[unitquantreg]{ugompertz}:} implements the \[d-q-p-r\]ugompertz function for the unit-Gompertz distribution.
#'
#' \code{\link[unitquantreg]{ulogistic}:} implements the \[d-q-p-r\]ulogistic function for the unit-Logistic distribution.
#'
#' \code{\link[unitquantreg]{uweibull}:} implements the \[d-q-p-r\]uweibull function for the unit-Weibull distribution.
#'
#' @author
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @useDynLib unitquantreg
#' @importFrom Rcpp sourceCpp
NULL
