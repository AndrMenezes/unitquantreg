#' @name ugompertz
#' @aliases ugompertz dugompertz pugompertz qugompertz rugompertz
#' @importFrom stats runif
#'
#' @title The unit-Gompertz distribution
#'
#' @description Density function, distribution function, quantile function and random number deviates
#' for the unit-Gompertz distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Mazucheli, J., Menezes, A. F. and Dey, S., (2019). Unit-Gompertz Distribution with Applications. \emph{Statistica}, \bold{79}(1), 25-43.
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
#' @return \code{dugompertz} gives the density, \code{pugompertz} gives the distribution function,
#' \code{qugompertz} gives the quantile function and \code{rugompertz} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\alpha \theta }{x}\exp \left\{ \alpha -\theta \log \left( y\right) -\alpha \exp \left[ -\theta \log \left( y\right) \right] \right\} }
#'
#' Cumulative density function
#' \deqn{F(y\mid \alpha ,\theta )=\exp \left[ \alpha \left( 1-y^{\theta }\right) \right] }
#'
#' Quantile Function
#' \deqn{Q(\tau \mid \alpha ,\theta )=\left[ \frac{\alpha -\log \left( \tau \right) }{\alpha }\right] ^{-\frac{1}{\theta }} }
#'
#' Reparameterization
#' \deqn{\alpha =g^{-1}(\mu )=\frac{\log \left( \tau \right) }{1-\mu ^{\theta }}}
#'
#' @examples
#' set.seed(123)
#' x <- rugompertz(n = 1000, mu = 0.5, theta = 2, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Gompertz')
#' lines(S, dugompertz(x = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pugompertz(q = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qugompertz(p = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
##################################################
#' @rdname ugompertz
#' @export
#
dugompertz <- function (x, mu, theta, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dugompertz(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ugompertz
#' @export
#'
pugompertz <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pugompertz(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ugompertz
#' @export
#'
qugompertz <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qugompertz(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ugompertz
#' @export
#'
rugompertz <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qugompertz(runif(n), mu, theta, tau, TRUE, FALSE)
}
