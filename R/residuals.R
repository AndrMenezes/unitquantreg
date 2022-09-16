#' @title Residuals method for \code{unitquantreg} objects
#'
#' @description Extract various types of residuals from unit quantile regression models.
#'
#' @author André F. B. Menezes
#'
#' @param object fitted model object of class \code{\link[unitquantreg]{unitquantreg}}.
#' @param type character indicating type of residuals. The options are
#' \code{"quantile"}, \code{"cox-snell"}, \code{"working"} and \code{"partial"}.
#' @param ... currently not used.
#'
#' @details The \code{\link[stats]{residuals}} method can compute quantile
#' and Cox-Snell residuals. These residuals are defined, respectively, by
#'
#' \deqn{r_{Q} = \Phi^{-1}\left[ F(y_i \mid \widehat{\mu}_i, \widehat{\theta}_i)\right]}
#'
#' and
#'
#' \deqn{r_{CS} = -\log\left[1- F(y_i \mid \widehat{\mu}_i, \widehat{\theta}_i)\right]}
#' where \eqn{\widehat{\mu}_i} and \eqn{\widehat{\theta}_i} are the fitted values
#' of parameters \eqn{\mu} and \eqn{\theta}, \eqn{F(\cdot \mid \cdot, \cdot)} is
#' the cumulative distribution function (c.d.f.) and \eqn{\Phi(\cdot)} is the
#' c.d.f. of standard Normal distribution.
#'
#' Apart from the variability due the estimates of parameters,if the fitted
#' regression model is correctly specified then the quantile
#' residuals, \eqn{r_Q}, follow a standard Normal distribution and
#' the Cox-Snell residuals, \eqn{r_{CS}}, follow a standard exponential
#' distribution.
#'
#' @return Numeric vector of residuals extract from an object of class
#' \code{\link[unitquantreg]{unitquantreg}}.
#'
#'
#' @references
#' Cox, D. R. and Snell E. J., (1968). A general definition of residuals. \emph{Journal of the Royal Statistical Society - Series B}, \bold{30}(2), 248--265.
#'
#' Dunn, P. K. and Smyth, G. K., (1996). Randomized quantile residuals. \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3), 236--244.
#'
#' @rdname residuals.unitquantreg
#' @export
residuals.unitquantreg <- function(object, type = c("quantile", "cox-snell",
                                                    "working", "partial"),
                                   ...) {

  mu <- object$fitted.values$mu
  theta <- object$fitted.values$theta
  tau <- object$tau
  y <- object$y
  type <- match.arg(type)

  if (type %in% c("quantile", "cox-snell")) {
    abbrev <- .get_abbrev(object$family, fname = FALSE)
    pfun <- match.fun(paste0("p", abbrev))
    parms <- list(q = y, mu = mu, theta = theta, tau = tau)
    Fy <- do.call(pfun, parms)
    res <- switch (type,
      "cox-snell" = -log(1 - Fy),
      "quantile" =  qnorm(Fy))
  } else if (type %in% c("working", "partial")) {
    eta <- object$link$mu$linkfun
    d_eta_mu <- numDeriv::grad(func = eta, x = mu)
    r <- (y - mu) * d_eta_mu
    res <- switch (type,
      "working" = r,
      "partial" = r + predict(object, type = "terms"))
  }

  res
}
