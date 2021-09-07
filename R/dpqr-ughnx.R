#' @importFrom stats runif
#' @name ughnx
#' @aliases ughnx dughnx pughnx qughnx rughnx
#' @importFrom stats dnorm qnorm
#'
#' @title The unit-Half-Normal-X distribution
#'
#' @description Density function, distribution function, quantile function and random number generation function
#' for the unit-Half-Normal-X distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
#'
#' @author
#'
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#' Bakouch, H. S., Nik, A. S., Asgharzadeh, A. and Salinas, H. S., (2021). A flexible probability model for  proportion  data: Unit-Half-Normal  distribution. \emph{Communications  in  Statistics: CaseStudies, Data Analysis and Applications}, \bold{0}(0), 1--18.
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
#' @return \code{dughnx} gives the density, \code{pughnx} gives the distribution function,
#' \code{qughnx} gives the quantile function and \code{rughnx} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @details
#' Probability density function
#' \deqn{f(y\mid \alpha ,\theta )=\sqrt{\frac{2}{\pi }}\frac{\theta }{y\left(1-y\right) }\left( {\frac{y}{\alpha \left( 1-y\right) }}\right) ^{\theta }\mathrm{\exp }\left\{ -\frac{1}{2}\left[ {\frac{y}{\alpha \left( 1-y\right) }}\right] ^{2\theta }\right\}}
#'
#' Cumulative density function
#' \deqn{F(y\mid \alpha ,\theta )=2\Phi \left[ \left( \frac{y}{\alpha \left(1-y\right) }\right) ^{\theta }\right] -1}
#'
#' Quantile Function
#' \deqn{Q(\tau \mid \alpha )=\frac{\alpha \left[ \Phi ^{-1}\left( \frac{\tau +1}{2}\right) \right] ^{\frac{1}{\theta }}}{1+\alpha \left[ \Phi ^{-1}\left( \frac{ \tau +1}{2}\right) \right] ^{\frac{1}{\theta }}}}
#'
#' Reparametrization
#' \deqn{\alpha=g^{-1}(\mu )=\frac{\mu }{\left( 1-\mu \right) \left[ \Phi ^{-1}\left( \frac{\tau +1}{2}\right) \right] ^{\frac{1}{\theta }}}}
#'
#' @examples
#' set.seed(123)
#' x <- rughnx(n = 1000, mu = 0.5, theta = 2, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Half-Normal-X')
#' lines(S, dughnx(x = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pughnx(q = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qughnx(p = S, mu = 0.5, theta = 2, tau = 0.5), col = 2)
#'
##################################################
#' @rdname ughnx
#' @export
#
dughnx <- function (x, mu, theta, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0);
  cpp_dughnx(x, mu, theta, tau, log[1L]);
}
##################################################
#' @rdname ughnx
#' @export
#'
pughnx <- function (q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_pughnx(q, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ughnx
#' @export
#'
qughnx <- function(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qughnx(p, mu, theta, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname ughnx
#' @export
#'
rughnx <- function(n, mu, theta, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, tau > 0, tau < 1, theta > 0)
  cpp_qughnx(runif(n), mu, theta, tau, TRUE, FALSE)
}
