#' @importFrom stats runif
#' @name ubs
#' @aliases ubs dubs pubs qubs rubs
#'
#' @title The unit-Birnbaum-Saunders distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Birnbaum-Saunders distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#'
#' Birnbaum, Z. W. and Saunders, S. C., (1969). A new family of life distributions. \emph{Journal of Applied Probability}, \bold{6}(2), 637--652.
#
#' Mazucheli, J., Menezes, A. F. B. and Dey, S., (2018). The unit-Birnbaum-Saunders distribution with applications. \emph{Chilean Journal of Statistics}, \bold{9}(1), 47--57.
#'
#' Mazucheli, J., Alves, B. and Menezes, A. F. B., (2021). A new quantile regression for modeling bounded data under a unit Birnbaum-Saunders distribution with applications. \emph{Simmetry}, \bold{}(), 1--28.
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
#' @return \code{dubs} gives the density, \code{pubs} gives the distribution function,
#' \code{qubs} gives the quantile function and \code{rubs} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\frac{1}{2y\alpha \theta \sqrt{2\pi }}\left[\left( -\frac{\alpha }{\log (y)}\right) ^{\frac{1}{2}}+\left( -\frac{\alpha}{\log (y)}\right) ^{\frac{3}{2}}\right] \exp \left[ \frac{1}{2\theta ^{2}}\left( 2+\frac{\log (y)}{\alpha }+\frac{\alpha }{\log (y)}\right) \right]}
#'
#' Cumulative distribution function
#' \deqn{F(y\mid \alpha ,\theta )=1-\Phi \left\{ \frac{1}{\theta }\left[ \left( -\frac{\log (y)}{\alpha }\right) ^{\frac{1}{2}}-\left( -\frac{\alpha }{\log(y)}\right) ^{\frac{1}{2}}\right] \right\}}
#'
#' Quantile function
#' \deqn{Q\left( \tau \mid \alpha ,\theta \right) ={\exp }\left\{ -{\frac{2\alpha}{2+\left[ {\theta }\Phi ^{-1}\left( 1-\tau \right) \right] ^{2}-{\theta } \Phi ^{-1}\left( 1-\tau \right) \sqrt{4+\left[ {\theta }\Phi ^{-1}\left(1-\tau \right) \right] ^{2}}}}\right\}}
#'
#' Reparameterization
#' \deqn{\alpha=g^{-1}(\mu )=\log \left( \mu \right) g\left( \theta ,\tau \right)}
#' where \eqn{g\left( \theta ,\tau \right) =-\frac{1}{2}\left\{ 2+\left[ {\theta }\Phi^{-1}\left( 1-\tau \right) \right] ^{2}-{\theta }\Phi ^{-1}\left( 1-\tau\right) \sqrt{4+{\theta }\Phi ^{-1}\left( 1-\tau \right) }\right\} .}
#'
#' @examples
#' set.seed(123)
#' x <- rubs(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Birnbaum-Saunders')
#' lines(S, dubs(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pubs(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qubs(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
#'
##################################################
#' @rdname ubs
#' @export
#
dubs <- function(x, mu, theta, tau = 0.5, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_dubs(x, mu, theta, tau, log[1L])
}
##################################################
#' @rdname ubs
#' @export
#'
pubs <- function(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pubs(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ubs
#' @export
#'
qubs <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qubs(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ubs
#' @export
#'
rubs <- function(n, mu, theta, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qubs(runif(n), mu, theta, tau, TRUE, FALSE)
}
