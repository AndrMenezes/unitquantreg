

# arcsecant hyperbolic Weibull --------------------------------------------
ashw <- function(link = "logit", link.theta = "identity") {

  family <- "ashw"

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
    rashw(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# Johnson SB --------------------------------------------------------------
johnsonsb <- function(link = "logit", link.theta = "identity") {

  family <- "johnsonsb"

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
    rjohnsonsb(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# Kumaraswamy -------------------------------------------------------------
kum <- function(link = "logit", link.theta = "identity") {

  family <- "kum"

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
    rkum(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# Log-extended exponential-geometric --------------------------------------
leeg <- function(link = "logit", link.theta = "identity") {

  family <- "leeg"

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
    rleeg(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-BS -----------------------------------------------------------------

ubs <- function(link = "logit", link.theta = "identity") {

  family <- "ubs"

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
    rubs(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Burr-XII -----------------------------------------------------------

uburrxii <- function(link = "logit", link.theta = "identity") {

  family <- "uburrxii"

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
    ruburrxii(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Chen ---------------------------------------------------------------

uchen <- function(link = "logit", link.theta = "identity") {

  family <- "uchen"

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
    ruchen(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Half-Normal-E ------------------------------------------------------

ughne <- function(link = "logit", link.theta = "identity") {

  family <- "ughne"

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
    rughne(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Half-Normal-X ------------------------------------------------------

ughnx <- function(link = "logit", link.theta = "identity") {

  family <- "ughnx"

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
    rughnx(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Gompertz -----------------------------------------------------------

ugompertz <- function(link = "logit", link.theta = "identity") {

  family <- "ugompertz"

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
    rugompertz(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Gumbel -------------------------------------------------------------

ugumbel <- function(link = "logit", link.theta = "identity") {

  family <- "ugumbel"

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
    rugumbel(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}


# unit-Logistic -----------------------------------------------------------

ulogistic <- function(link = "logit", link.theta = "identity") {

  family <- "ulogistic"

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
    rulogistic(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}

# unit-Weibull ------------------------------------------------------------

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
    ruweibull(n = nsim * length(mu), mu = mu, theta = theta, tau = tau)
  }

  structure(
    list(family = family,
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         simulate = simfun,
         class = "bounded_family")
  )
}
