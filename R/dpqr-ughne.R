#' @importFrom stats runif
#' @name ughne
#' @aliases ughne dughne pughne qughne rughne
#'
#' @title The unit-Half-Normal-E distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Half-Normal-E distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Korkmaz, M. C., (2020). The unit generalized half normal distribution: A new bounded distribution with inference and application. \emph{University Politehnica of Bucharest Scientific}, \bold{82}(2), 133--140.
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
#' @return \code{dughne} gives the density, \code{pughne} gives the distribution function,
#' \code{qughne} gives the quantile function and \code{rughne} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\sqrt{\frac{2}{\pi }}\frac{\theta }{y\left[ -\log\left( y\right) \right] }\left( -{\frac{\log \left( y\right) }{\alpha }} \right)^{\theta }\mathrm{\exp }\left\{ -\frac{1}{2}\left[ -{\frac{\log \left( y\right) }{\alpha }}\right]^{2\theta }\right\}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=2\Phi \left[ -\left( -{\frac{\log \left( y\right) }{\alpha }}\right)^{\theta }\right]}
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha ,\theta )=\exp \left\{ -\alpha \left[ -\Phi^{-1}\left(\frac{\tau }{2}\right) \right]^{\frac{1}{\theta }}\right\}}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )=-\log \left( \mu \right) \left[ -\Phi^{-1}\left( \frac{\tau }{2}\right) \right]^{-\frac{1}{\theta }}}
#'
#' @examples
#' set.seed(123)
#' x <- rughne(n = 1000, mu = 0.5, theta = 2, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Half-Normal-E')
#' lines(S, dughne(x = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pughne(q = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qughne(p = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#'
##################################################
#' @rdname ughne
#' @export
#
dughne <- function (x, mu, theta, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dughne(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ughne
#' @export
#'
pughne <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pughne(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ughne
#' @export
#'
qughne <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qughne(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ughne
#' @export
#'
rughne <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qughne(runif(n), mu, theta, tau, TRUE, FALSE)
}
