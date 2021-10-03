#' @importFrom stats runif
#' @name ugumbel
#' @aliases ugumbel dugumbel pugumbel qugumbel rugumbel
#'
#' @title The unit-Gumbel distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Gumbel distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli
#'
#' Andre F. B. Menezes
#'
#' @references
#'
#' Mazucheli, J. and Alves, B., (2021). The unit-Gumbel Quantile Regression Model for Proportion Data. \emph{Under Review}.
#'
#' Gumbel, E. J., (1941). The return period of flood flows. \emph{The Annals of Mathematical Statistics}, \bold{12}(2), 163--190.
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#' @param theta nonnegative shape parameter.
#' @param tau the parameter to specify which quantile use in the parametrization.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dugumbel} gives the density, \code{pugumbel} gives the distribution function,
#' \code{qugumbel} gives the quantile function and \code{rugumbel} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\theta }{y(1-y)}\exp \left\{ -\alpha -\theta \log \left( \frac{y}{1-y}\right) -\exp \left[ -\alpha -\theta \log \left( \frac{y}{1-y}\right) \right] \right\}}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid\alpha,\theta)={\exp }\left[ -{{\exp }}\left( -\alpha \right)\left( \frac{1-y}{y}\right) ^{\theta } \right] }
#'
#' Quantile function
#' \deqn{Q(\tau \mid \alpha, \theta)= \frac{\left [-\frac{1}{\log(\tau)  }\right ]^{\frac{1}{\theta}}}{\exp\left ( \frac{\alpha}{\theta} \right )+\left [-\frac{1}{\log(\tau)  }\right ]^{\frac{1}{\theta}}}}
#'
#' Reparameterization
#' \deqn{\alpha =  g^{-1}(\mu ) =\theta \log \left( {\frac{1-\mu }{\mu }}\right) +\log \left( -\frac{1}{\log \left( \tau \right) }\right)}
#'
#' where \eqn{0<y<1} and \eqn{\theta >0} is the shape parameter.
#'
#' @examples
#' set.seed(6969)
#' x <- rugumbel(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-Gumbel')
#' lines(S, dugumbel(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pugumbel(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qugumbel(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname ugumbel
#' @export
#

dugumbel <- function (x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dugumbel(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ugumbel
#' @export
#'

pugumbel <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pugumbel(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ugumbel
#' @export
#'

qugumbel <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qugumbel(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ugumbel
#' @export
#'
rugumbel <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qugumbel(runif(n), mu, theta, tau, TRUE, FALSE)
}
