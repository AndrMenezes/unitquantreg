#' @name likelihood_stats
#'
#' @title Likelihood-based statistics of fit for \code{unitquantreg} objects.
#'
#' @description Computes the likelihood-based statistics (Neg2LogLike, AIC, BIC and HQIC)
#' from \code{unitquantreg} objects.
#'
#' @param ... \code{\link{unitquantreg}} objects separated by commas.
#' Not use in \code{print} method.
#' @param lt a list with one or more \code{\link{unitquantreg}} objects.
#' @param x object of class \code{likelihood_stats} obtained from \code{likelihood_stats}
#' function.
#'
#' @details
#' Neg2LogLike: The log-likelihood is reported as \deqn{Neg2LogLike= -2\log(L)}
#'
#' AIC: The Akaike's information criterion (AIC) is defined as \deqn{AIC = -2\log(L)+2p}
#'
#' BIC: The Schwarz Bayesian information criterion (BIC) is defined as \deqn{BIC =  -2\log(L) + p\log(n)}
#'
#' HQIC: The Hannan and Quinn information criterion (HQIC)  is defined as \deqn{HQIC =  -2\log(L) + 2p\log[\log(n)]}
#' where \eqn{L} is the likelihood function.
#'
#' @return A list with class \code{"likelihood_stats"} containing the following
#' components:
#' \item{call}{the matched call.}
#' \item{stats}{ordered matrix according AIC value containg the likelihood
#' based statistics.}
#'
#'
#' @author
#' André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification. \emph{IEEE Transaction on Automatic Control}, \bold{19}(6), 716--723.
#'
#' Hannan, E. J. and Quinn, B. G. (1979). The determination of the order of an autoregression. \emph{Journal of the Royal Statistical Society, Series B}, \bold{41}(2), 190--195.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, \bold{6}(2), 461--464.
#'
#'
#' @examples
#' data(water, package = "unitquantreg")
#'
#' models  <- c("johnsonsb", "kum", "leeg", "ubs", "uburrxii", "uchen", "ughne", "ughnx",
#' "ugompertz", "ulogistic", "uweibull")
#' fits <- lapply(models, function(M) unitquantreg(formula = phpws ~ mhdi + incpc +
#' region + log(pop),tau = 0.5, data = water, family = M))
#'
#' ans <- likelihood_stats(lt = fits)
#' ans
#'
#' @importFrom stats AIC BIC logLik
#'
#' @rdname likelihood_stats
#' @export


likelihood_stats <- function(..., lt = NULL) {

  if (is.null(lt)) lt <-  list(...)

  n <- length(lt[[1L]]$y)
  lt_mat <- lapply(lt, function(x) {
    matrix(c(-2 * logLik(x), AIC(x), BIC(x), AIC(x, k = 2 * log(log(n)))), nrow = 1,
           dimnames = list(x$family, c("Neg2LogLike", "AIC", "BIC", "HQIC")))
  })
  mat <- do.call(rbind, lt_mat)
  mat <- mat[order(mat[, 2]), ]

  out <- list(call = match.call(), stats = mat)
  class(out) <- "likelihood_stats"
  out
}

#' @rdname likelihood_stats
#' @export
#'
print.likelihood_stats <- function(x, ...) {

  cr <- .FF(x[["stats"]], 3)

  cat("\n Likelihood-based statistics of fit for unit quantile regression models \n", sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  # cat("Likelihood-based statistics: \n")
  print(cr, quote = FALSE)
}
