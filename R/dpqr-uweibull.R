#' @importFrom stats runif
#' @name uweibull
#' @aliases uweibull duweibull puweibull quweibull ruweibull
#'
#' @title The unit-Weibull distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Weibull distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli
#'
#' André F. B. Menezes
#'
#' @references
#'
#' Mazucheli, J., Menezes, A. F. B and Ghitany, M. E., (2018). The unit-Weibull distribution and associated inference. \emph{Journal of Applied Probability and Statistics}, \bold{13}(2), 1--22.
#'
#' Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P. and Ghitany, M. E., (2020). The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the modeling of quantiles conditional on covariates. \emph{Journal of Applied Statistics}, \bold{47}(6), 954--974.
#'
#' Mazucheli, J., Menezes, A. F. B., Alqallaf, F. and Ghitany, M. E., (2021). Bias-Corrected Maximum Likelihood Estimators of the Parameters of the Unit-Weibull Distribution. \emph{Austrian Journal of Statistics}, \bold{50}(3), 41--53.
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
#' @return \code{duweibull} gives the density, \code{puweibull} gives the distribution function,
#' \code{quweibull} gives the quantile function and \code{ruweibull} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{\alpha \theta }{y}\left[ -\log (y)\right]^{\theta -1}\exp \left\{ -\alpha \left[ -\log (y)\right]^{\theta }\right\} }
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=\exp \left\{ -\alpha \left[ -\log (y)\right]^{\theta }\right\}}
#'
#' Quantile function
#' \deqn{Q\left( \tau \mid \alpha ,\theta \right) =\exp \left\{ -\left[ -\frac{\log (\tau )}{\alpha }\right]^{\frac{1}{\theta }}\right\}}
#'
#' Reparameterization
#' \deqn{\alpha =g^{-1}(\mu )=-\frac{\log (\tau )}{[-\log (\mu )]^{\theta}}}
#'
#' @examples
#' set.seed(6969)
#' x <- ruweibull(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-Weibull')
#' lines(S, duweibull(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, puweibull(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(quweibull(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
##################################################
#' @rdname uweibull
#' @export
#

duweibull <- function (x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_duweibull(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname uweibull
#' @export
#'

puweibull <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_puweibull(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uweibull
#' @export
#'

quweibull <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quweibull(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname uweibull
#' @export
#'
ruweibull <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_quweibull(runif(n), mu, theta, tau, TRUE, FALSE)
}
