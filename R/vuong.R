#' @name vuong.test
#'
#' @title Vuong test
#'
#' @description Performs Vuong test between two fitted objects of class
#' \code{\link[unitquantreg]{unitquantreg}}
#'
#'
#' @param object1,object2 objects of class \code{\link[unitquantreg]{unitquantreg}}
#' containing the fitted models.
#' @param alternative indicates the alternative hypothesis and must be one
#' of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}. You can
#' specify just the initial letter of the value, but the argument name must be
#' given in full. See ‘Details’ for the meanings of the possible values.
#'
#'
#' @details The statistic of Vuong likelihood ratio test for compare two
#' non-nested regression models is defined by
#' \deqn{T = \frac{1}{\widehat{\omega}^2\,\sqrt{n}}\,\sum_{i=1}^{n}\,
#' \log\frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{
#' g(y_i \mid \boldsymbol{x}_i,\widehat{\boldsymbol{\gamma}})}}
#' where
#' \deqn{\widehat{\omega}^2 = \frac{1}{n}\,\sum_{i=1}^{n}\,\left(\log \frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}\right)^2 - \left[\frac{1}{n}\,\sum_{i=1}^{n}\,\left(\log \frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{ g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}\right)\right]^2}
#' is an estimator for the variance of
#' \eqn{\frac{1}{\sqrt{n}}\,\displaystyle\sum_{i=1}^{n}\,\log\frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}},
#' \eqn{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})} and
#' \eqn{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}
#' are the corresponding rival densities evaluated at the maximum likelihood estimates.
#'
#' When \eqn{n \rightarrow \infty} we have that \eqn{T \rightarrow N(0, 1)} in distribution.
#' Therefore, at \eqn{\alpha\%} level of significance  the null hypothesis of
#' the equivalence of the competing models is rejected if \eqn{|T| > z_{\alpha/2}},
#' where \eqn{z_{\alpha/2}} is the \eqn{\alpha/2} quantile of standard normal distribution.
#'
#' In practical terms, \eqn{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}
#' is better (worse) than \eqn{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}
#' if \eqn{T>z_{\alpha/2}} (or \eqn{T< -z_{\alpha/2}}).
#'
#'
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string with the method used.}
#' \item{data.name}{a character string ginven the name of families models under comparison.}
#'
#' @author
#' André F. B. Menezes
#'
#' Josmar Mazucheli
#'
#' @references
#' Vuong, Q. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, \bold{57}(2), 307--333.
#'
#' @examples
#' data(sim_bounded, package = "unitquantreg")
#' sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]
#'
#' fit_uweibull <- unitquantreg(formula = y1 ~ x, tau = 0.5,
#'                              data = sim_bounded_curr,
#'                              family = "uweibull")
#' fit_kum <- unitquantreg(formula = y1 ~ x, tau = 0.5,
#'                              data = sim_bounded_curr,
#'                              family = "kum")
#'
#' ans <- vuong.test(object1 = fit_uweibull, object2 = fit_kum)
#' ans
#' str(ans)
#'
#'
#' @importFrom stats var pnorm
#'
#' @rdname vuong.test
#' @export

vuong.test <- function(object1, object2, alternative = c("two.sided", "less", "greater")) {

  alternative <- match.arg(alternative)

  y <- object1$y
  n <- length(y)

  # Fitted parms
  mu_1 <- object1$fitted.values$mu
  mu_2 <- object2$fitted.values$mu
  theta_1 <- object1$fitted.values$theta
  theta_2 <- object2$fitted.values$theta

  # Quantile
  tau_1 <- object1$tau
  tau_2 <- object2$tau
  if (tau_1 != tau_2)
    warning("Comparison is done between fitted models for different quantiles!")

  # Get pdf of models
  dfun_1 <- match.fun(paste0("d", .get_abbrev(object1$family, fname = FALSE)))
  parms_1 <- list(x = y, mu = mu_1, theta = theta_1, tau = tau_1, log = TRUE)
  dfun_2 <- match.fun(paste0("d", .get_abbrev(object2$family, fname = FALSE)))
  parms_2 <- list(x = y, mu = mu_2, theta = theta_2, tau = tau_2, log = TRUE)

  # Compute log-pdfs
  ll_1 <- do.call(dfun_1, parms_1)
  ll_2 <- do.call(dfun_2, parms_2)

  # Compute vuong statistics
  om2 <- (n - 1) / n * var(ll_1 - ll_2)
  lr <- sum(ll_1 - ll_2)
  Tstat <- (1/sqrt(n)) * lr/sqrt(om2)
  names(Tstat) <- "T_LR"

  # Compute p-value according to hypothesis
  pvalue <- switch(alternative,
                   "two.sided" = 2 * pnorm(abs(Tstat), lower.tail = FALSE),
                   "greater" = pnorm(abs(Tstat), lower.tail = FALSE),
                   "less" = pnorm(abs(Tstat), lower.tail = TRUE)
  )

  # Output
  method <- paste0("Vuong likelihood ratio test for non-nested models (",
                   object1$family, " versus ", object2$family, ")")
  out <- list(statistic = Tstat, p.value = pvalue,
              method = method,
              data.name = paste0(object1$family, " versus ", object2$family))
  class(out) <- "htest"
  out
}

#' @name pairwise.vuong.test
#'
#' @title Pairwise vuong test
#'
#' @description Calculate pairwise comparisons between fitted models performing
#' vuong test for objects of class \code{\link[unitquantreg]{unitquantreg}}.
#'
#' @param ... \code{\link[unitquantreg]{unitquantreg}} objects separated by commas.
#' @param lt a list with one or more \code{\link[unitquantreg]{unitquantreg}} objects.
#' @param p.adjust.method a character string specifying the method for multiple
#' testing adjustment; almost always one of
#' \code{p.adjust.methods}. Can be abbreviated.
#' @param alternative indicates the alternative hypothesis and must be one
#' of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' Can be abbreviated.
#'
#'
#' @return Object of class \code{"pairwise.htest"}
#'
#' @seealso \code{\link{vuong.test}}, \code{\link{p.adjust}}
#'
#' @examples
#' data(sim_bounded, package = "unitquantreg")
#' sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]
#'
#' models <- c("uweibull", "kum", "ulogistic")
#' lt_fits <- lapply(models, function(fam) {
#'   unitquantreg(formula = y1 ~ x, tau = 0.5, data = sim_bounded_curr,
#'                family = fam)
#' })
#'
#' ans <- pairwise.vuong.test(lt = lt_fits)
#' ans
#'
#' @importFrom stats pairwise.table p.adjust.methods
#'
#' @rdname pairwise.vuong.test
#' @export
#'


pairwise.vuong.test <- function(..., lt, p.adjust.method = p.adjust.methods,
                                alternative = c("two.sided", "less", "greater")) {

  if (is.null(lt)) lt <-  list(...)

  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)

  families <- sapply(lt, "[[", "family")
  nfam <- length(families)

  compare.levels <- function(i, j) {
    xi <- lt[[as.integer(i)]]
    xj <- lt[[as.integer(j)]]
    vuong.test(object1 = xi, object2 = xj, alternative = alternative)$p.value
  }
  pval <- pairwise.table(compare.levels = compare.levels,
                         level.names = seq_along(families),
                         p.adjust.method = p.adjust.method)
  colnames(pval) <- families[1L:(nfam - 1)]
  rownames(pval) <- families[2L:nfam]

  # Output
  method <- "Vuong likelihood ratio test for non-nested models"
  dname <- paste(deparse(lt[[1L]]$call, width.cutoff = 60L), collapse = "\n")
  out <- list(method = method, data.name = dname, p.value = pval,
              p.adjust.method = p.adjust.method)
  class(out) <- "pairwise.htest"
  out
}
