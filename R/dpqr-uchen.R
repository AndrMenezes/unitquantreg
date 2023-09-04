#' @importFrom stats runif
#' @name uchen
#' @aliases uchen duchen puchen quchen ruchen
#'
#' @title The unit-Chen distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Chen distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Korkmaz, M. C., Emrah, A., Chesneau, C. and Yousof, H. M., (2020). On the unit-Chen distribution with associated quantile regression and applications. \emph{Journal of Applied Statistics}, \bold{44}(1) 1--22.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#' @param theta nonnegative shape parameter.
#' @param tau the parameter to specify which quantile is to be used.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{duchen} gives the density, \code{puchen} gives the distribution function,
#' \code{quchen} gives the quantile function and \code{ruchen} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\alpha \theta }{y}\left[ -\log (y)\right]^{\theta -1}\exp \left\{ \left[ -\log \left( y\right) \right]^{\theta}\right\} \exp \left\{ \alpha \left\{ 1-\exp \left[ \left( -\log (y)\right)^{\theta }\right] \right\} \right\}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=\exp \left\{ \alpha \left\{ 1-\exp \left[ \left(-\log (y)\right)^{\theta }\right] \right\} \right\}}
#'
#' Quantile function
#' \deqn{Q\left( \tau \mid \alpha ,\theta \right) =\exp \left\{ -\left[ \log \left( 1-{\frac{\log \left( \tau \right) }{\alpha }}\right) \right]^{\frac{1}{\theta}}\right\}}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )={\frac{\log \left( \tau \right) }{1-\exp \left[ \left( -\log (\mu )\right)^{\theta }\right]}}}
#'
#' @examples
#' set.seed(123)
#' x <- ruchen(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Chen')
#' lines(S, duchen(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, puchen(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(quchen(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname uchen
#' @export
#
duchen <- function(x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_duchen(x, mu, theta, tau, log[1L])
}
##################################################
#' @rdname uchen
#' @export
#'
puchen <- function(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_puchen(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uchen
#' @export
#'
quchen <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quchen(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uchen
#' @export
#'
ruchen <- function(n, mu, theta, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quchen(runif(n), mu, theta, tau, TRUE, FALSE)
}
