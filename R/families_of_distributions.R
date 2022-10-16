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
  if (link.theta %in% exp_links__theta) {
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

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}

