#' @importFrom stats runif
#' @name kum
#' @aliases kum dkum pkum qkum rkum
#'
#' @title The Kumaraswamy distribution
#'
#' @description Density function, distribution function, quantile function and random number generation for the Kumaraswamy distribution  reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli
#'
#' André F. B. Menezes
#'
#' @references
#'
#' Kumaraswamy, P., (1980). A generalized probability density function for double-bounded random processes. \emph{Journal of Hydrology}, \bold{46}(1), 79--88.
#'
#' Jones, M. C., (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages. \emph{Statistical Methodology}, \bold{6}(1), 70-81.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#' @param theta nonnegative shape parameter.
#' @param tau the parameter to specify which quantile is to used.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.
#
#' @return \code{dkum} gives the density, \code{pkum} gives the distribution function,
#' \code{qkum} gives the quantile function and \code{rkum} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\alpha \theta y^{\theta -1}(1-y^{\theta })^{\alpha-1}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=1-\left( 1-y^{\theta }\right) ^{\alpha }}
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha ,\theta )=\left[ 1-\left( 1-\tau \right) ^{\frac{1}{\alpha }}\right] ^{\frac{1}{\theta }}}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )=\frac{\log (1-\tau )}{\log (1-\mu ^{\theta })}}
#'
#' @examples
#' set.seed(123)
#' x <- rkum(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'Kumaraswamy')
#' lines(S, dkum(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pkum(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qkum(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname kum
#' @export
#
dkum <- function (x, mu, theta, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dkum(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname kum
#' @export
#'
pkum <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pkum(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname kum
#' @export
#'
qkum <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qkum(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname kum
#' @export
#'
rkum <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qkum(runif(n), mu, theta, tau, TRUE, FALSE)
}
