#' @name hnp
#'
#' @title (Half-)Normal probability plots with simulated envelopes for
#' \code{\link{unitquantreg}} objects
#'
#' @description Produces a (half-)normal probability plot from a fitted model
#' object of class \code{\link{unitquantreg}}.
#'
#' @author
#' André F. B. Menezes
#'
#' @param object fitted model object of class \code{\link{unitquantreg}}.
#' @param nsim number of simulations used to compute envelope. Default is 99.
#' @param resid.type type of residuals to be used. The default is \code{quantile}.
#' See \code{\link{residuals.unitquantreg}} for further details.
#' @param halfnormal logical. If \code{TRUE}, a half-normal plot is produced.
#' If \code{FALSE}, a normal plot is produced.
#' @param plot Should the (half-)normal plot be plotted? Default is \code{TRUE}.
#' @param output Should the output be returned? Default is \code{TRUE}.
#' @param level confidence level of the simulated envelope. Default is 0.95.
#' @param ... currently not used.
#'
#' @details Residuals plots with simulated envelope were proposed by Atkinson (1981)
#' and can be construct as follows:
#'
#' \enumerate{
#'    \item{generate sample set of \eqn{n} independent observations from the estimated
#'    parameters of the fitted model;}
#'    \item{fit the model using the generated sample, if \code{halfnormal} is
#'    \code{TRUE} compute the absolute values of the residuals and arrange them in order;}
#'    \item{repeat steps (1) and (2) \code{nsim} number of times;}
#'    \item{consider the \eqn{n} sets of the \code{nsim} ordered statistics
#'    of the residuals, then for each set compute the quantile \eqn{\code{level}/2},
#'    the median and the quantile \eqn{1 - \code{level}/2};}
#'    \item{plot these values and the ordered residuals of the original sample set
#'    versus the expected order statistics of a (half)-normal distribution,
#'    which is approximated as}
#'    \deqn{G^{-1} \left(\frac{i + n - 0.125}{2n + 0.5} \right)}
#'    for half-normal plots, i.e., \code{halfnormal=TRUE} or
#'    \deqn{G^{-1} \left(\frac{i - 0.375}{n + 0.25}\right)}
#'    for normal plots, i.e., \code{halfnormal=FALSE}, where \eqn{G(\cdot)} is the the
#'    cumulative distribution function of standard Normal distribution for
#'    \code{quantile} residuals or the standard exponential distribution for the
#'    \code{cox-snell} residuals.
#' }
#'
#' According to Atkinson (1981), if the model was correctly specified then no
#' more than \code{level}100% of the observations are expected to appear
#' outside the envelope bands. Additionally, if a large proportion of the
#' observations lies outside the envelope, thus one has evidence against
#' the adequacy of the fitted model.
#'
#' @return A list with the following components in ordered
#' (and absolute if \code{halfnormal} is \code{TRUE}) values:
#' \item{obs}{the observed residuals.}
#' \item{teo}{the theoretical residuals.}
#' \item{lower}{lower envelope band.}
#' \item{median}{median envelope band.}
#' \item{upper}{upper envelope band.}
#' \item{time_elapsed}{time elapsed to fit the \code{nsim} models.}
#'
#' @seealso \code{\link{residuals.unitquantreg}}
#'
#' @references
#' Atkinson, A. C., (1981). Two graphical displays for outlying and influential observations in regression. \emph{Biometrika} \bold{68}(1), 13--20.
#'
#'
#' @importFrom stats qnorm qexp residuals quantile median coef
#' @importFrom graphics plot lines
#'
#' @rdname hnp
#' @export
#'
hnp <- function(object, ...){
  UseMethod("hnp", object)
}

#' @rdname hnp
#' @export
#'
hnp.unitquantreg <- function(object, nsim = 99, halfnormal = TRUE, plot = TRUE,
                             output = TRUE, level = 0.95,
                             resid.type = c("quantile", "cox-snell"), ...) {

  resid.type <- match.arg(resid.type)
  mu <- object$fitted.values$mu
  theta  <- object$fitted.values$theta
  tau <- object$tau
  family <- object$family
  init <- coef(object)
  link <- object$link$mu$name
  link.theta <- object$link$theta$name

  control <- object$control
  control$hessian <- control$starttests <- control$dowarn <- FALSE
  control$trace <- 0L

  # Extract data
  X <- if (is.null(object$x)) model.matrix(object, type = "quantile") else object$x$quantile
  Z <- if (is.null(object$x)) model.matrix(object, type = "shape") else object$x$shape
  n <- nrow(X)

  # Get the random number generation
  fam_abbrev <- .get_abbrev(object$family, fname = FALSE)
  rfun <- match.fun(paste0("r", fam_abbrev))
  parms <- list(n = n * nsim, mu = mu, theta = theta, tau = tau)

  # Generating data
  ysim <- do.call(rfun, parms)
  ysim <- matrix(ysim, nrow = n, ncol = nsim)
  # apply(ysim, 2, median)
  # median(mu)

  # Computing residuals
  ini <- proc.time()
  res_sim <- sapply(seq_len(nsim), function(j) {
    out <- rep(NaN, n)
    ycur <- ysim[, j]
    try({
      fit <- suppressWarnings(
        unitquantreg.fit(y = ycur, X = X, Z = Z, tau = tau, family = fam_abbrev,
                         link = link, link.theta = link.theta, control = control,
                         start = init)
      )
      fit$y <- ycur
      class(fit) <- "unitquantreg"
      out <- residuals.unitquantreg(fit, type = resid.type)
    }, silent = TRUE)
    return(out)
  })
  fim <- proc.time()

  # Putting NA in infinite values
  is.na(res_sim) <- sapply(res_sim, is.infinite)

  if (halfnormal) {
    res_obs <- sort(abs(residuals(object, type = resid.type)))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n + n - 1/8) / (2 * n + 0.5))
    }
    else if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n + n - 1/8) / (2 * n + 0.5))
    }
  } else {
    res_obs <- sort(residuals(object, type = resid.type))
    res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n - 3 / 8) / (n + 1 / 4))
    }
    else if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n - 3 / 8) / (n + 1 / 4))
    }
  }

  alpha   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)

  if (plot) {
    yl <- ifelse(resid.type == "quantile", "Quantile residuals",
                 "Cox-Snell residuals")

    Ry <- c(min(res_lwr), max(res_upr))
    Rx <- range(res_teo)

    plot(x = res_teo, y = res_obs, xlab = 'Theoretical quantiles',
         ylab = yl, xlim = Rx, ylim = Ry, bty = 'o', pch = 3, ...)
    lines(x = res_teo, y = res_lwr)
    lines(x = res_teo, y = res_upr)
    lines(x = res_teo, y = res_mid, lty = 2)
  }

  if (output) {
    return(list(obs = res_obs, teo = res_teo, lower = res_lwr,
                median = res_mid, upper = res_upr,
                time_elapsed = unname(fim - ini)[3]))
  }
  invisible()
}
