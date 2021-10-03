
uweibull <- function(link = "logit", link.theta = "identity") {

  linkobj <- if (is.character(link)) make.link(link) else link
  linkobj.theta <- if (is.character(link.theta)) make.link(link.theta) else link.theta

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

  structure(
    list(family = "uweibull",
         linkobj.mu = linkobj,
         linkobj.theta = linkobj.theta,
         gradient = cpp_gradientuweibull,
         hessian = cpp_hessianuweibull,
         simulate = simfun,
         variance = varfun,
         class = "bounded_family")
  )
}
