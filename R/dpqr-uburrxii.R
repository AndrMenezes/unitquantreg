#' @importFrom stats runif
#' @name uburrxii
#' @aliases uburrxii duburrxii puburrxii quburrxii ruburrxii
#'
#' @title The unit-Burr-XII distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Burr-XII distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Korkmaz M. C. and Chesneau, C., (2021). On the unit Burr-XII distribution with the quantile regression modeling and applications. \emph{Computational and Applied Mathematics}, \bold{40}(29), 1--26.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#' @param theta nonnegative shape parameter.
#' @param tau the parameter to specify which quantile is to used.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{duburrxii} gives the density, \code{puburrxii} gives the distribution function,
#' \code{quburrxii} gives the quantile function and \code{ruburrxii} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha, \theta )=\frac{\alpha \theta }{y}\left[ -\log (y)\right]^{\theta -1}\left\{ 1+\left[ -\log (y)\right] ^{\theta }\right\} ^{-\alpha -1}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha, \theta )=\left\{ 1+\left[ -\log (y)\right] ^{\theta}\right\} ^{-\alpha }}
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha, \theta )=\exp \left[ -\left( \tau ^{-\frac{1}{\alpha }}-1\right)^{\frac{1}{\theta }} \right]}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu)=\frac{\log\left ( \tau^{-1} \right )}{\log\left [ 1+\log\left ( \frac{1}{\mu} \right )^\theta \right ]}}
#'
#' @examples
#' set.seed(123)
#' x <- ruburrxii(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Burr-XII')
#' lines(S, duburrxii(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, puburrxii(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(quburrxii(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#'
##################################################
#' @rdname uburrxii
#' @export
#
duburrxii <- function(x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_duburrxii(x, mu, theta, tau, log[1L])
}
##################################################
#' @rdname uburrxii
#' @export
#'
puburrxii <- function(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_puburrxii(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uburrxii
#' @export
#'
quburrxii <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quburrxii(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uburrxii
#' @export
#'
ruburrxii <- function(n, mu, theta, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quburrxii(runif(n), mu, theta, tau, TRUE, FALSE)
}

