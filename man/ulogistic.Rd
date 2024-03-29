% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpqr-ulogistic.R
\name{ulogistic}
\alias{ulogistic}
\alias{dulogistic}
\alias{pulogistic}
\alias{qulogistic}
\alias{rulogistic}
\title{The unit-Logistic distribution}
\usage{
dulogistic(x, mu, theta, tau = 0.5, log = FALSE)

pulogistic(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

qulogistic(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

rulogistic(n, mu, theta, tau = 0.5)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{mu}{location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.}

\item{theta}{nonnegative shape parameter.}

\item{tau}{the parameter to specify which quantile is to used.}

\item{log, log.p}{logical; If TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\value{
\code{dulogistic} gives the density, \code{pulogistic} gives the distribution function,
\code{qulogistic} gives the quantile function and \code{rulogistic} generates random deviates.

Invalid arguments will return an error message.
}
\description{
Density function, distribution function, quantile function and random number generation for the unit-Logistic distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
}
\details{
Probability density function
\deqn{f(y\mid \alpha ,\theta )=\frac{\theta \exp \left( \alpha \right) \left(\frac{y}{1-y}\right) ^{\theta -1}}{\left[ 1+\exp \left( \alpha \right)\left( \frac{y}{1-y}\right) ^{\theta }\right] ^{2}}}

Cumulative distribution function
\deqn{F(y\mid \alpha ,\theta )=\frac{\exp \left( \alpha \right) \left( \frac{y}{1-y}\right) ^{\theta }}{1+\exp \left( \alpha \right) \left( \frac{y}{1-y}\right) ^{\theta }}}

Quantile function
\deqn{Q(\tau \mid \alpha ,\theta )=\frac{\exp \left( -\frac{\alpha }{\theta }\right) \left( \frac{\tau }{1-\tau }\right) ^{\frac{1}{\theta }}}{1+\exp\left( -\frac{\alpha }{\theta }\right) \left( \frac{\tau }{1-\tau }\right) ^{ \frac{1}{\theta }}} }

Reparameterization
\deqn{\alpha=g^{-1}(\mu )=\log \left( \frac{\tau }{1-\tau }\right) -\theta \log \left(   \frac{\mu }{1-\mu }\right) }
}
\examples{
set.seed(123)
x <- rulogistic(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
R <- range(x)
S <- seq(from = R[1], to = R[2], by =  0.01)
hist(x, prob = TRUE, main = 'unit-Logistic')
lines(S, dulogistic(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(ecdf(x))
lines(S, pulogistic(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(quantile(x, probs = S), type = "l")
lines(qulogistic(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
}
\references{
Paz, R. F., Balakrishnan, N. and Bazán, J. L., 2019. L-Logistic regression models: Prior sensitivity analysis, robustness to outliers and applications. \emph{Brazilian Journal of Probability and Statistics}, \bold{33}(3), 455--479.
}
\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}

André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
