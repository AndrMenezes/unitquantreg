#' @importFrom stats runif
#' @name ashw
#' @aliases ashw dashw pashw qashw rashw
#'
#' @title The arcsecant hyperbolic Weibull distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the arcsecant hyperbolic Weibull distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli
#'
#' Bruna Alves
#'
#' @references
#'
#' Korkmaz, M. C., Chesneau, C. and Korkmaz, Z. S., (2021). A new alternative quantile regression model for the bounded response with educational measurements applications of OECD countries. \emph{Journal of Applied Statistics}, 1--25.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#' @param theta shape parameter.
#' @param tau the parameter to specify which quantile use in the parametrization.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dashw} gives the density, \code{pashw} gives the distribution function,
#' \code{qashw} gives the quantile function and \code{rashw} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y;\alpha, \theta)=\frac{\alpha \theta}{y\sqrt{1-y^2}} \mathrm{arcsech}(y)^{\theta-1}\exp\left [ -\alpha \mathrm{arcsech}(y)^\theta \right ]}
#'
#' Cumulative distribution function
#' \deqn{F(y;\alpha, \theta)=\exp\left [ -\alpha \mathrm{arcsech}(y)^\theta \right ]}
#'
#' Quantile function
#' \deqn{Q(\tau;\alpha, \theta)= \mathrm{sech}\left \{ \left [ -\alpha^{-1} \log(\tau)\right ]^{\frac{1}{\theta}} \right \}}
#'
#' Reparameterization
#' \deqn{\alpha = g^{-1}(\mu) = -\frac{\log(\tau)}{\mathrm{arcsech}(\mu)^\theta}}
#'
#' where \eqn{\theta >0} is the shape parameter and \eqn{\mathrm{arcsech}(y)= \log\left[\left( 1+\sqrt{1-y^2} \right)/y \right]}.
#'
#' @examples
#' set.seed(6969)
#' x <- rashw(n = 1000, mu = 0.5, theta = 2.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1L], to = R[2L], by = 0.01)
#' hist(x, prob = TRUE, main = 'arcsecant hyperbolic Weibull')
#' lines(S, dashw(x = S, mu = 0.5, theta = 2.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pashw(q = S, mu = 0.5, theta = 2.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qashw(p = S, mu = 0.5, theta = 2.5, tau = 0.5), col = 2)

#' @rdname ashw
#' @export
#

dashw <- function (x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1);
  cpp_dashw(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ashw
#' @export
#'

pashw <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1);
  cpp_pashw(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ashw
#' @export
#'

qashw <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1);
  cpp_qashw(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ashw
#' @export
#'
rashw <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1);
  cpp_qashw(runif(n), mu, theta, tau, TRUE, FALSE)
}
