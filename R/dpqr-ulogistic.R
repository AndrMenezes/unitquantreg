#' @importFrom stats runif
#' @name ulogistic
#' @aliases ulogistic dulogistic pulogistic qulogistic rulogistic
#'
#' @title The unit-Logistic distribution
#'
#' @description Density function, distribution function, quantile function and random number generation for the unit-Logistic distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#'
#' Paz, R. F., Balakrishnan, N. and Bazán, J. L., 2019. L-Logistic regression models: Prior sensitivity analysis, robustness to outliers and applications. \emph{Brazilian Journal of Probability and Statistics}, \bold{33}(3), 455--479.
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
#' @return \code{dulogistic} gives the density, \code{pulogistic} gives the distribution function,
#' \code{qulogistic} gives the quantile function and \code{rulogistic} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\theta \exp \left( \alpha \right) \left(\frac{y}{1-y}\right) ^{\theta -1}}{\left[ 1+\exp \left( \alpha \right)\left( \frac{y}{1-y}\right) ^{\theta }\right] ^{2}}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=\frac{\exp \left( \alpha \right) \left( \frac{y}{1-y}\right) ^{\theta }}{1+\exp \left( \alpha \right) \left( \frac{y}{1-y}\right) ^{\theta }}}
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha ,\theta )=\frac{\exp \left( -\frac{\alpha }{\theta }\right) \left( \frac{\tau }{1-\tau }\right) ^{\frac{1}{\theta }}}{1+\exp\left( -\frac{\alpha }{\theta }\right) \left( \frac{\tau }{1-\tau }\right) ^{ \frac{1}{\theta }}} }
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )=\log \left( \frac{\tau }{1-\tau }\right) -\theta \log \left(   \frac{\mu }{1-\mu }\right) }
#'
#' @examples
#' set.seed(123)
#' x <- rulogistic(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Logistic')
#' lines(S, dulogistic(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pulogistic(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qulogistic(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname ulogistic
#' @export
#
dulogistic <- function(x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dulogistic(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ulogistic
#' @export
#'
pulogistic <- function(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pulogistic(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ulogistic
#' @export
#'
qulogistic <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qulogistic(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ulogistic
#' @export
#'
rulogistic <- function(n, mu, theta, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qulogistic(runif(n), mu, theta, tau, TRUE, FALSE)
}
