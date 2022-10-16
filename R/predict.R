#' @title Prediction method for \code{unitquantreg} class
#'
#' @description Extract various types of predictions from unit quantile regression models.
#'
#' @param object fitted model object of class \code{\link{unitquantreg}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param type character indicating type of predictions. The options are
#' \code{link}, \code{quantile}, \code{shape} and \code{terms}.
#' @param interval type of interval desired. The options are \code{none} and
#' \code{confidence}. The "\code{terms}" option returns a matrix giving
#' the fitted values of each term in the model formula on the linear predictor
#' scale.
#' @param level coverage probability for the confidence intervals. Default is
#' \code{0.95}.
#' @param se.fit logical. If \code{TRUE} return the asymptotic standard errors.
#' @param ... currently not used.
#'
#' @return
#' If \code{se.fit = FALSE} then returns a \code{data.frame} with
#' predict values and confidence interval if \code{interval = TRUE}.
#'
#' If \code{se.fit = TRUE} returns a list with components:
#' \item{fit}{Predictions, as for \code{se.fit = FALSE}.}
#' \item{se.fit}{Estimated standard errors.}
#'
#' For \code{type = "terms"} the output is a \code{data.frame} with a columns
#' per term.
#'
#' @author
#' André F. B. Menezes
#'
#' @importFrom stats qnorm terms
#' @importFrom Formula as.Formula
#'
#' @rdname predict.unitquantreg
#' @export
#'
predict.unitquantreg <- function(object, newdata,
                                 type = c("link", "quantile", "shape", "terms"),
                                 interval = c("none", "confidence"),
                                 level = 0.95, se.fit = FALSE, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  vcov <- object$vcov
  p <- length(object$coefficients$mu)
  q <- length(object$coefficients$theta)
  ff <- object$formula

  if (type == "shape" & object$theta_const) {
    stop("It is not possible predict the shape parameter, because the model has constant shape")
  }

  # Get the design matrix
  if (missing(newdata)) {
    if (type == "shape") {
      Z <- if (is.null(object$x)) model.matrix(object, type = "shape") else object$x$shape
    }
    else {
      X <- if (is.null(object$x)) model.matrix(object, type = "quantile") else object$x$quantile
    }
  }
  else {
    if (type == "shape") {
      Z <- model.matrix(ff, newdata, rhs = 2L)
    }
    else {
      X <- model.matrix(ff, newdata, rhs = 1L)
    }
  }

  if (type %in% c("link", "quantile")) {
    beta <- object$coefficients$mu
    linkobj.mu <- object$link$mu
    vcov.beta <- vcov[1:p, 1:p]
    g.mu.hat <- drop(X %*% beta)
    if (type == "link") {
      out <- cbind(fit = g.mu.hat)
      J <- g.mu.hat * X ### check delta method
    }
    else if (type == "quantile") {
      mu.hat <- linkobj.mu$linkinv(g.mu.hat)
      out <- cbind(fit = mu.hat)
      J <- linkobj.mu$mu.eta(mu.hat) * X
    }
    if (se.fit | interval == "confidence") {
      variance <- tcrossprod(J %*% vcov.beta, J)
      stderror <- sqrt(diag(variance))
      if (interval == "confidence") {
        qa <- qnorm(1 - level / 2)
        out <- cbind(out, lower = out[, "fit"] - qa * stderror,
                     upper = out[, "fit"] + qa * stderror)
      }
    }
  }

  if (type == "shape") {
    gamma <- object$coefficients$theta
    linkobj.theta <- object$link$theta
    vcov.gamma <- vcov[1:p + q, 1:p + q]
    g.theta.hat <- drop(Z %*% gamma)
    theta.hat <- linkobj.theta$linkinv(g.theta.hat)
    out <- cbind(fit = theta.hat)
    if (se.fit | interval == "confidence") {
      J <- linkobj.theta$mu.eta(theta.hat) * Z
      variance <- tcrossprod(J %*% vcov.gamma, J)
      stderror <- sqrt(diag(variance))
      if (interval == "confidence") {
        qa <- qnorm(1 - level / 2)
        out <- cbind(out, lower = out[, "fit"] - qa * stderror,
                     upper = out[, "fit"] + qa * stderror)
      }
    }
  }

  if (type == "terms") {
    betas <- object$coefficients$mu
    X1 <- sweep(X, 2L, colMeans(X))
    terms_betas <- labels(terms(ff))
    terms_betas <- if (object$theta_const) terms_betas else terms_betas[1L:(p - 1)]
    out <- t(betas * t(X1))[, -1]
    out <- as.data.frame(out)
    names(out) <- terms_betas
    return(out)
  }

  if (!missing(newdata)) {
    out <- as.data.frame(cbind(newdata, out))
  }

  if (se.fit) {
    out <- list(fit = out, se.fit = stderror)
  }

  out
}
