#' @title Plot method for \code{unitquantreg} objects
#'
#' @description Provide diagnostic plots to check model assumptions for fitted model
#' of class \code{unitquantreg}.
#'
#' @param x fitted model object of class \code{unitquantreg}.
#' @param which integer. if a subset of the plots is required, specify a subset
#' of the numbers 1 to 4, see below for further details.
#' @param caption character. Captions to appear above the plots.
#' @param sub.caption character. Common title-above figures if there are multiple.
#' @param main character. Title to each plot in addition to the above caption.
#' @param ask logical. If \code{TRUE}, the user is asked before each plot.
#' @param ... other parameters to be passed through to plotting functions.
#' @param add.smooth logical. Indicates if a smoother should be added to most plots
#' @param type character. Indicates type of residual to be used, see
#' \code{\link{residuals.unitquantreg}}.
#' @param nsim integer. Number of simulations in half-normal plots, see
#' \code{\link{hnp.unitquantreg}}.
#' @param level numeric. Confidence level of the simulated envelope, see
#' \code{\link{hnp.unitquantreg}}.
#'
#' @details
#' The \code{plot} method for \code{unitquantreg} objects produces four types
#' of diagnostic plot.
#'
#' The \code{which} argument can be used to select a subset of currently four
#' supported plot, which are: Residuals versus indices of observations
#' (\code{which = 1}); Residuals versus linear predictor (\code{which = 2});
#' Working response versus linear predictor (\code{which = 3}) to
#' check possible misspecification of link function; Half-normal plot of
#' residuals (\code{which = 4}) to check distribution assumption.
#'
#' @return
#' No return value, called for side effects.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (2018) Generalized Linear Models With Examples in R, Springer, New York.
#'
#' @seealso \code{\link{residuals.unitquantreg}},
#' \code{\link{hnp.unitquantreg}},
#' \code{\link{unitquantreg}}.
#'
#' @importFrom graphics par plot mtext abline lines title
#' @importFrom grDevices dev.interactive
#' @importFrom stats predict lowess
#'
#' @author
#' André F. B. Menezes

#' @rdname plot.unitquantreg
#' @export
#'
plot.unitquantreg <- function(x, which = 1L:4L,
                              caption = c(
                                "Residuals vs. indices of obs.",
                                "Residuals vs. linear predictor",
                                "Working response vs. linear predictor",
                                "Half-normal plot of residuals"),
                              sub.caption = paste(
                                deparse(x$call), collapse = "\n"),
                              main = "",
                              ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                              ...,
                              add.smooth = getOption("add.smooth"),
                              type = "quantile",
                              nsim = 99L,
                              level = 0.95) {

  if (!is.numeric(which) || any(which < 1L) || any(which > 4L))
    stop("'which' must be in 1:4")

  types <- c("quantile", "cox-snell")
  Types <- c("Quantile residuals", "Cox-Snell residuals")
  type <- match.arg(type, types)
  Type <- Types[type == types]

  res <- residuals(x, type = type)
  eta <- predict(x, type = "link")
  n <- length(res)
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  Main <- rep("", 4)
  Main[which] <- rep(main, length.out = sum(show))
  one_fig <- prod(par("mfcol")) == 1

  # Keeping user's graphs options
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1L]) {
    plot(seq_len(n), res, xlab = "Obs. number", ylab = Type,
         main = Main[1L], ...)
    if (one_fig) title(sub = sub.caption, ...)
    mtext(caption[1L], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
    if (add.smooth) lines(lowess(1:n, res), col = "blue", lwd = 2)
  }
  if (show[2L]) {
    plot(eta, res,
         xlab = "Linear predictor", ylab = Type, main = Main[2L], ...)
    if (one_fig) title(sub = sub.caption, ...)
    mtext(caption[2L], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
    if (add.smooth) lines(lowess(eta, res), col = "blue", lwd = 2)
  }
  if (show[3L]) {
    z <- eta + residuals(x, type = "working")
    plot(z ~ eta, xlab = "Linear predictor", ylab = "Working responses",
         main = Main[3L], ...)
    if (one_fig) title(sub = sub.caption, ...)
    mtext(caption[3L], 3, 0.25)
    if (add.smooth) lines(lowess(z ~ eta), lwd = 2, col = "blue")
  }
  if (show[4L]) {
    hn <- hnp(x, resid.type = type, nsim = nsim, halfnormal = TRUE,
              plot = FALSE, level = level)
    r_y <- c(min(hn[["lower"]]), max(hn[["upper"]]))
    r_x <- range(hn[["teo"]])

    plot(hn[["teo"]], hn[["obs"]], ylim = r_y, xlim = r_x, main = Main[5],
         xlab = "Theoretical quantiles",
         ylab = paste(Type, "(absolute values)"), ...)
    lines(hn[["teo"]], y = hn[["lower"]])
    lines(hn[["teo"]], y = hn[["upper"]])
    lines(hn[["teo"]], y = hn[["median"]], lty = 2)
    if (one_fig) title(sub = sub.caption, ...)
    mtext(caption[4L], 3, 0.25)
  }
  if (!one_fig && par("oma")[3] >= 1)
    mtext(sub.caption, outer = TRUE, cex = 1.25)

  invisible()
}


#' @title Plot method for \code{unitquantregs} objects
#'
#' @description Provide two type of plots for \code{unitquantregs} objects.
#'
#' @param x fitted model object of class \code{unitquantregs}.
#' @param which character. Indicate the type of plot. Currently supported are \code{"coef"} which
#' provide the estimated coefficients for several quantiles and \code{"conddist"} which provide
#' the conditional distribution (cdf or pdf) at specific values of covariates.
#' @param output_df logical. Should \code{data.frame} used to plot be returned?
#' @param parm a specification of which parameters are to be plotted, either a vector
#' of numbers or a vector of names. By default, all parameters are considered.
#' @param level level of significance for the confidence interval of parameters.
#' @param mean_effect logical. Should a line for the mean effect coefficients be added?
#' @param mfrow,mar,ylim,main,col,border,cex,pch,type,xlab,ylab graphical parameters.
#' @param dist_type character. Which conditional distribution should be plotted?
#' The options are \code{"density"} or \code{"cdf"}.
#' @param at_avg logical. Should consider the conditional distribution at average values of covariates?
#' @param at_obs list. List with name and values for each covariate.
#' @param legend_position character. The legend position argument used in \code{legend} function.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @details
#' The plot method for \code{unitquantregs} objects is inspired in PROC QUANTREG of SAS/STAT.
#' This plot method provide two type of visualizations.
#'
#' If \code{which = "coef"} plot the estimated coefficients for several quantiles.
#'
#' If \code{which = "conddist"} plot the conditional distribution at specific values of
#' covariates. The conditional distribution could be the cumulative distribution function
#' if \code{dist_type = "cdf"} or the probability density function if \code{dist_type = "pdf"}.
#'
#' @return
#' If \code{output_df = TRUE} then returns a data.frame used to plot.
#' Otherwise, no return value, called for side effects.
#'
#' @seealso \code{\link{plot.unitquantreg}}.
#'
#' @importFrom graphics par plot mtext abline lines title grid polygon points legend
#' @importFrom grDevices n2mfrow gray
#' @importFrom stats predict coef confint
#'
#' @author
#' André F. B. Menezes
#'
#' @rdname plot.unitquantregs
#' @export
#'

plot.unitquantregs <- function(x, which = c("coef", "conddist"), output_df = FALSE,
                               parm = NULL, level = 0.95, mean_effect = FALSE,
                               mfrow = NULL, mar = NULL, ylim = NULL, main = NULL,
                               col = gray(c(0, 0.75)), border = NULL, cex = 1, pch = 20,
                               type = "b", xlab = bquote("Quantile level ("*tau*")"),
                               ylab = "Estimate effect", dist_type = c("density", "cdf"),
                               at_avg = TRUE, at_obs = NULL, legend_position = "topleft", ...) {
  if (which == "coef") {
    .plot_coef(x = x, output_df = output_df, parm = parm, level = level,
               mean_effect = mean_effect, mfrow = mfrow, mar = mar, ylim = ylim,
               main = main, col = col, border = border, cex = cex, pch = pch,
               type = type, xlab = xlab, ylab = ylab, ...)
  }
  if (which == "conddist") {
    .plot_conddist(x = x, output_df = output_df, dist_type = dist_type,
                   at_avg = at_avg, at_obs = at_obs,
                   legend_position = legend_position, ...)
  }
}
