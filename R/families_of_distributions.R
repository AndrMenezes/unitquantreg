.linear_predictors_two_parms <- function(parms, X, Z) {
  p <- dim(X)[2L]
  betas <- parms[seq.int(length.out = p)]
  gammas <- parms[-seq.int(length.out = p)]
  list(eta_mu = drop(X %*% betas), zeta_theta = drop(Z %*% gammas))
}

uweibull <- function(link = "logit", link.theta = "identity") {

  family <- "uweibull"

  # Checking link function
  exp_links__mu <- c("logit", "probit", "cloglog", "cauchit")
  exp_links__theta <- c("identity", "log", "sqrt")

  if (link %in% exp_links__mu) {
    linkobj <- if (is.character(link)) make.link(link) else link
  } else {
    stop(gettextf(
      "link \"%s\" not available for %s family; available links are %s",
      link, family, paste(sQuote(exp_links__mu), collapse = ", ")),
      domain = NA)
  }
  if (link %in% exp_links__theta) {
    linkobj.theta <- if (is.character(link.theta)) make.link(link.theta) else link.theta
  } else {
    stop(gettextf(
      "link.theta \"%s\" not available for %s family; available links are %s",
      link, family, paste(sQuote(exp_links__theta), collapse = ", ")),
      domain = NA)
  }

  simfun <- function(object, nsim) {
    mu <- object$fitted.values$mu
    theta  <- object$fitted.values$theta
    tau <- object$tau
    ruweibull(nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  varfun <- function(mu, theta, tau, MC = FALSE) {
    if (MC) {
      x <- ruweibull(1e5, mu, theta, tau)
      stats::var(x)
    } else {
      alpha <- -log(tau) / (-log(mu))^theta
      n <- seq.int(0, 100)
      e <- sum((-1)^n / (factorial(n) * alpha^(n / theta)) * gamma(n / theta + 1))
      e2 <- sum((-2)^n / (factorial(n) * alpha^(n / theta)) * gamma(n / theta + 1))
      e2 - e^2
    }
  }

  # Minus log-likelihood
  loglike <- function(parms, tau, linkobj, linkobj.theta, X, Z, y) {

    # Computing the linear predictors (eta and zeta)
    l <- .linear_predictors_two_parms(parms, X, Z)
    mu <- linkobj$linkinv(l[["eta_mu"]])
    theta <- linkobj.theta$linkinv(l[["theta"]])

    # Compute the minus of log-likelihood
    cpp_loglikeuweibull(x = y, lnx = log(y), n = length(y),
                        mu = mu, theta = theta, tau = tau)
  }

  # Gradient function
  gradient <- score_unitquantreg
  gradient <- function(parms, tau, linkobj, linkobj.theta, X, Z, y) {
    n <- length(y)
    # Computing the linear predictors (eta and zeta)
    l <- .linear_predictors_two_parms(parms, X, Z)
    mu <- linkobj$linkinv(l[["eta_mu"]])
    dmu_deta <- linkobj$mu.eta(l[["eta_mu"]])
    theta <- linkobj.theta$linkinv(l[["theta"]])
    dtheta_dzeta <- linkobj.theta$mu.eta(l[["theta"]])

    # Compute the gradient/score function
    score <- cpp_gradientuweibull(
      n = n, x = y, U = matrix(0, nrow = n, ncol = 2),
      dmu_deta = dmu_deta, dtheta_dzeta = dtheta_dzeta,
      mu = mu, theta =  theta, tau = tau)
    dbetas <- crossprod(X, score[, 1L])
    dthetas <- crossprod(Z, score[, 2L])

    - 1L * c(dbetas, dthetas)
  }

  # Hessian function
  hessian <- function(parms, tau, linkobj, linkobj.theta, X, Z, y) {

    n <- length(y)
    # Computing the linear predictors (eta and zeta)
    l <- .linear_predictors_two_parms(parms, X, Z)
    mu <- linkobj$linkinv(l[["eta_mu"]])
    dmu_deta <- linkobj$mu.eta(l[["eta_mu"]])
    theta <- linkobj.theta$linkinv(l[["theta"]])
    dtheta_dzeta <- linkobj.theta$mu.eta(l[["theta"]])

    # Compute second derivatives of log-likelihood function
    W <- cpp_hessianuweibull(
      n = n, x = y, H = matrix(0, nrow = n, ncol = 3),
      mu = mu, theta =  theta, tau = tau)

    # Diagonal matrix
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

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         loglike = loglike,
         gradient = gradient,
         hessian = hessian,
         simulate = simfun,
         variance = varfun,
         class = "bounded_family")
  )
}
