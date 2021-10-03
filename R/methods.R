#' @name methods-unitquantreg
#'
#' @title Methods for \code{unitquantreg} and \code{unitquantregs} objects
#'
#' @description Methods for extracting information from fitted unit quantile regression
#' objects of class \code{unitquantreg} and \code{unitquantregs}.
#'
#' @author André F. B. Menezes
#'
#' @param object,x fitted model object of class \code{\link{unitquantreg}}.
#' @param digits  minimal number of _significant_ digits.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed. Default is \code{FALSE}.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type character indicating type of fitted values to return.
#' @param formula. Changes to the formula see \code{\link{update.formula}} for details.
#' @param formula an R formula.
#' @param evaluate If true evaluate the new call else return the call.
#' @param ... additional argument(s) for methods. Currently not used.
#'
#' @importFrom stats pnorm cov2cor coef vcov printCoefmat update formula
#'
#' @rdname methods-unitquantreg
#' @export

print.unitquantreg <- function(x, digits = max(4, getOption("digits") - 3), ...) {

  cat("\n", x$family, " quantile regression model \n", sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients (quantile model with ", x$link$mu$name, " link and tau = ", x$tau,  "): \n", sep = "")

  print.default(.FF(x$coefficients$mu, Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  if (x$theta_const) {
    cat("Model with constant shape parameter:", "\n", sep = "")
    print.default(.FF(x$coefficients$theta, Digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
  } else {
    names(x$coefficients$theta) <- gsub("(theta)_", "", names(x$coefficients$theta), fixed = TRUE)
    cat("Theta coefficients (shape model with ", x$link$theta$name, " link):", "\n", sep = "")
    print.default(.FF(x$coefficients$theta, Digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }

  invisible(x)
}

# Summary -----------------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
summary.unitquantreg <- function(object, correlation = FALSE, ...) {

  cf <- object$coefficients
  names(cf) <- NULL
  estimates <- unlist(cf)
  stderror <- sqrt(diag(object$vcov))
  zvalue <- estimates/stderror
  pvalue <- 2 * pnorm(-abs(zvalue))
  table <- cbind("Estimate" = estimates, "Std. Error" = stderror,
                 "Z value" = zvalue, "Pr(>|z|)" = pvalue)
  if (correlation) {
    correlation <- cov2cor(object$vcov)
  }

  out <- list(coeftable   = table,
              loglik      = object$loglik,
              df.residual = object$df.residual,
              correlation = correlation,
              call        = object$call,
              iterations  = object$iterations,
              tau         = object$tau,
              family      = object$family,
              link        = object$link,
              dims        = c(length(cf[[1L]]), length(cf[[2L]])),
              theta_const = object$theta_const)
  class(out) <- "summary.unitquantreg"
  out
}


# Print output summary ----------------------------------------------------

#' @export
print.summary.unitquantreg <- function(x, digits = max(4, getOption("digits") - 3), ...) {

  p <- x$dims[[1L]]
  q <- x$dims[[2L]]

  cat("\n Wald-tests for ", x$family, " quantile regression model", "\n",
      sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")

  cat("Mu coefficients: (quantile model with ", x$link$mu$name, " link and tau = ",
      x$tau,  "): \n", sep = "")
  printCoefmat(x$coeftable[1:p, , drop = FALSE], digits = digits,
               has.Pvalue = TRUE)
  cat("\n")

  if (x$theta_const) {
    cat("Model with constant shape:", "\n", sep = "")
    printCoefmat(x$coeftable[-(1:p), , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  } else {
    rownames(x$coeftable)[-(1:p)] <- gsub("(theta)_", "", rownames(x$coeftable)[-(1:p)], fixed = TRUE)
    cat("Theta coefficients (shape model with ", x$link$theta$name, " link):", "\n",
        sep = "")
    printCoefmat(x$coeftable[-(1:p), , drop = FALSE], digits = digits,
                 has.Pvalue = TRUE)
    cat("\n")
  }

  if (is.matrix(x$correlation)) {
    cat("Correlation of coefficients:", "\n", sep = "")
    corr <- x$correlation
    corr <- format(round(corr, 2L), nsmall = 2L, digits = digits)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -ncol(corr), drop = FALSE], quote = FALSE)
    cat("\n")
  }

  cat("Residual degrees of freedom: ", x$df.residual, "\n", sep = "")
  cat("Log-likelihood: ", x$loglik, "\n", sep = "")
  cat("Number of iterations: ", x$iterations, "\n", sep = "")

  invisible(x)
}

# coef function -----------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
coef.unitquantreg <- function(object, type = c("full", "quantile", "shape"), ...) {
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  cf <- object$coefficients
  type <- match.arg(type)
  out <- switch (type,
    "full" = {
      names(cf) <- NULL
      unlist(cf)
    },
    "quantile" = cf$mu,
    "shape" = cf$theta
  )
  out
}

# vcov function -----------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
vcov.unitquantreg <- function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  object$vcov
}


# logLik function ---------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
logLik.unitquantreg <- function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  ll <- object$loglik
  attr(ll, "df")   <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}

# confint function --------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
confint.unitquantreg <- function(object, parm, level = 0.95, ...)
{
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  pnames <- names(ses)
  if (missing(parm)) {
    parm <- pnames
  }
  else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level) * 0.5
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- .format_perc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}


# Fitted values -----------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
fitted.unitquantreg <- function(object, type = c("quantile", "shape", "full"),  ...) {

  if (!missing(...)) {
    warning("Extra arguments discarded")
  }

  type <- match.arg(type)

  switch (type,
          "all" = object$fitted.values,
          "quantile"  = object$fitted.values$mu,
          "shape" = object$fitted.values$theta
  )
}


# Terms -------------------------------------------------------------------


#' @rdname methods-unitquantreg
#' @export
terms.unitquantreg <- function(x, type = c("quantile", "shape"), ...) {
  x$terms[[match.arg(type)]]
}

#' @rdname methods-unitquantreg
#' @export
model.frame.unitquantreg <- function(formula, ...) {
#  if (!is.null(formula$model)) return(formula$model)
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]

  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    ## need stats:: for non-standard evaluation
    fcall[[1L]] <- quote(unitquantreg::unitquantreg)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model

}

#' @rdname methods-unitquantreg
#' @export
model.matrix.unitquantreg <- function(object, type = c("quantile", "shape"), ...) {
  type <- match.arg(type)
  if (is.null(object$x[[type]])) {
    model.matrix(object$terms[[type]], model.frame(object))
  } else {
    object$x[[type]]
  }
}

# Update ------------------------------------------------------------------

#' @rdname methods-unitquantreg
#' @export
update.unitquantreg <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if (is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    call$formula <- formula(update(Formula::Formula(formula(object)), formula.))
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) eval(call, parent.frame())
  else call
}


# Simulate --------------------------------------------------------------------------

# simulate.uniquantreg <- function(object, nsim = 1, seed = NULL) {
#   val <- object$family$simulate(object, nsim)
#   dim(val) <- c(length(object$y), nsim)
#   val <- as.data.frame(val)
#   names(val) <- paste("sim", seq_len(nsim), sep="_")
#   val
# }

# Methods unitquantregs (tau vectorized) ----------------------------------

#' @rdname methods-unitquantreg
#' @export

print.unitquantregs <- function(x, digits = max(3, getOption("digits") - 3), ...) {


  cat("\n", x[[1L]]$family, " quantile regression model \n", sep = "")

  cat("\nCall:  ", paste(deparse(x[[1L]]$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  quant_coefs <- sapply(x, coef, type = "quantile")
  shape_coefs <- sapply(x, coef, type = "shape")
  taus <- .FF(sapply(x, "[[", "tau"), Digits = digits)

  colnames(quant_coefs) <- paste0("tau = ", taus)

  cat("Mu coefficients (quantile model with ", x[[1L]]$link$mu$name, " link): \n", sep = "")

  print.default(.FF(quant_coefs, Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  cat("\n")

  if (x[[1L]]$theta_const) {
    shape_coefs <- matrix(shape_coefs, nrow = 1)
    colnames(shape_coefs) <- paste0("tau = ", taus)
    rownames(shape_coefs) <- names(coef(x[[1L]], type = "shape"))

    cat("Model with constant shape parameter:", "\n", sep = "")
    print.default(.FF(shape_coefs, Digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")

  } else {
    colnames(shape_coefs) <- paste0("tau = ", taus)

    cat("Theta coefficients (shape model with ", x[[1L]]$link$theta$name, " link):", "\n", sep = "")
    print.default(.FF(shape_coefs, Digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }

  invisible(x)
}

#' @rdname methods-unitquantreg
#' @export

summary.unitquantregs <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  out <- lapply(object, summary)
  names(out) <- .FF(sapply(object, "[[", "tau"), Digits = digits)
  out
}

