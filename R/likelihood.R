#' @title Log-likelihood, score vector and hessian matrix.
#'
#' @description Internal functions using in \code{\link[unitquantreg]{unitquantreg.fit}}
#' to compute the negative log-likelihood function, the score vector and the hessian
#' matrix using analytic expressions implemented in \code{C++}.
#
#' @param par vector of regression model coefficients for \eqn{\mu} and/or
#' \eqn{\theta}.
#' @param family specify the distribution family name.
#' @param tau quantile level, value between 0 and 1.
#' @param linkobj,linkobj.theta a function, usually obtained from
#' \code{\link{make.link}} for link function of \eqn{\mu} and \eqn{\theta},
#' respectively.
#' @param X design matrix related to the \eqn{\mu} parameter.
#' @param Z design matrix related to the \eqn{\theta} parameter.
#' @param y vector of response variable.


#' @noMd
loglike_unitquantreg <- function(par, tau, family, linkobj, linkobj.theta, X, Z, y) {

  # Utils
  rc <- dim(X)
  n <- rc[1L]
  p <- rc[2L]

  # Location parameter (mu)
  beta <- par[seq.int(length.out = p)]
  mu <- linkobj$linkinv(drop(X %*% beta))

  # Shape parameter (theta)
  q <- ncol(Z)
  gamma <- par[-seq.int(length.out = p)]
  theta <- linkobj.theta$linkinv(drop(Z %*% gamma))

  # Get the minus log-likelihood function
  abbrev <- .get_abbrev(family)
  llfun <- paste0("cpp_loglike", abbrev)
  lny <- log(y)
  parms <- list(x = y, lnx = lny, n = n, mu = mu, theta = theta, tau = tau)

  # Compute the minus of log-likelihood
  ll <- do.call(llfun, parms)

  ll
}

#' @noMd
score_unitquantreg <- function(par, tau, family, linkobj, linkobj.theta, X, Z, y) {

  # Utils
  rc <- dim(X)
  n <- rc[1L]
  p <- rc[2L]

  # Location parameter (mu)
  beta <- par[seq.int(length.out = p)]
  eta_mu <- drop(X %*% beta)
  mu <- linkobj$linkinv(eta_mu)
  dmu_deta <- linkobj$mu.eta(eta_mu)

  # Shape parameter (theta)
  q <- ncol(Z)
  gamma <- par[-seq.int(length.out = p)]
  zeta_theta <- drop(Z %*% gamma)
  theta <- linkobj.theta$linkinv(zeta_theta)
  dtheta_dzeta <- linkobj.theta$mu.eta(zeta_theta)

  # Auxiliary
  U <- matrix(0, nrow = n, ncol = 2)

  # Get the gradient/score function
  abbrev <- .get_abbrev(family)
  gradfun <- paste0("cpp_gradient", abbrev)
  parms <- list(n = n, x = y, U = U, dmu_deta = dmu_deta,
                dtheta_dzeta = dtheta_dzeta, mu = mu, theta = theta, tau = tau)

  # Compute the gradient/score function
  score <- do.call(gradfun, parms)
  dbetas <- crossprod(X, score[, 1])
  dthetas <- crossprod(Z, score[, 2])

  -1L * c(dbetas, dthetas)
}


#' @noMd
hessian_unitquantreg <- function(par, tau, family, linkobj, linkobj.theta, X, Z, y) {

  # Utils
  rc <- dim(X)
  n <- rc[1L]
  p <- rc[2L]

  # Location parameter (mu)
  beta <- par[seq.int(length.out = p)]
  eta_mu <- drop(X %*% beta)
  mu <- linkobj$linkinv(eta_mu)

  # Shape parameter (theta)
  q <- ncol(Z)
  gamma <- par[-seq.int(length.out = p)]
  zeta_theta <- drop(Z %*% gamma)
  theta <- linkobj.theta$linkinv(zeta_theta)

  # Auxiliary to keep second derivatives
  W <- matrix(0, ncol = 3, nrow = n)

  # Get the hessian function
  abbrev <- .get_abbrev(family)
  hessfun <- paste0("cpp_hessian", abbrev)
  parms <- list(n = n, x = y, H = W, mu = mu, theta = theta, tau = tau)

  # Second derivatives of log-likelihood
  W <- do.call(hessfun, parms)

  # Diagonal matrix
  dmu_deta <- linkobj$mu.eta(eta_mu)
  dtheta_dzeta <- linkobj.theta$mu.eta(zeta_theta)
  w_bb <- dmu_deta^2 * W[, 1]
  w_bg <- dmu_deta * dtheta_dzeta * W[, 2]
  w_gg <- dtheta_dzeta^2 * W[, 3]

  # Hessian
  H_bb <- crossprod(X, w_bb * X)
  H_gg <- crossprod(Z, w_gg * Z)
  H_bg <- crossprod(X, w_bg * Z)
  H <- -1L * cbind(rbind(H_bb, t(H_bg)), rbind(H_bg, H_gg))

  H
}

