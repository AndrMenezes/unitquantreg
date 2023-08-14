#' @importFrom stats runif
#' @name leeg
#' @aliases leeg dleeg pleeg qleeg rleeg
#'
#' @title The Log-extended exponential-geometric distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the Log-extended exponential-geometric distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references Jodrá, P. and Jiménez-Gamero, M. D., (2020). A quantile regression model for bounded responses based on the exponential-geometric distribution. \emph{Revstat - Statistical Journal}, \bold{18}(4), 415--436.
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
#' @return \code{dleeg} gives the density, \code{pleeg} gives the distribution function,
#' \code{qleeg} gives the quantile function and \code{rleeg} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\theta \left( 1+\alpha \right) y^{\theta -1}}{\left( 1+\alpha y^{\theta }\right) ^{2}}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=\frac{\left( 1+\alpha \right) y^{\theta }}{1+\alpha y^{\theta }}}
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha ,\theta )=\left[ \frac{\tau }{1+\alpha \left( 1-\tau\right) }\right] ^{\frac{1}{\theta }}}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )=-\frac{1-\tau \mu ^{\theta }}{\left( 1-\tau \right) }}
#'
#' @examples
#' set.seed(123)
#' x <- rleeg(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'Log-extended exponential-geometric')
#' lines(S, dleeg(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pleeg(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qleeg(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname leeg
#' @export
#
dleeg <- function(x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_dleeg(x, mu, theta, tau, log[1L])
}
##################################################
#' @rdname leeg
#' @export
#'
pleeg <- function(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pleeg(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname leeg
#' @export
#'
qleeg <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qleeg(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname leeg
#' @export
#'
rleeg <- function(n, mu, theta, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qleeg(runif(n), mu, theta, tau, TRUE, FALSE)
}
