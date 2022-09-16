# Possible family names ---------------------------------------------------

.families <- c("uweibull", "unitweibull",
               "kum", "kumaraswamy",
               "ulogistic", "unitlogistic",
               "uchen", "unitchen",
               "ubs", "unitbs",
               "leeg", "logexpgeom",
               "ughne", "unitghne",
               "ughnx", "unitghnx",
               "ugompertz", "unitgompertz",
               "uburrxii", "unitburrxii",
               "jsb", "johnsonsb",
               "ashw", "arcsechwei",
               "ugumbel", "unitgumbel")

# Get family names --------------------------------------------------------

.get_family_name <- function(family) {

  if (family %in% c("uweibull", "unitweibull")) {
    "unit-Weibull"
  }
  else if (family %in% c("kum", "kumaraswamy")) {
    "Kumaraswamy"
  }
  else if (family %in% c("ulogistic", "unitlogistic")) {
    "unit-Logistic"
  }
  else if (family %in% c("uchen", "unitchen")) {
    "unit-Chen"
  }
  else if (family %in% c("ubs", "unitbs")) {
    "unit-Birnbaum-Saunders"
  }
  else if (family %in% c("leeg", "logexpgeom")) {
    "log-extended Exponential-Geometric"
  }
  else if (family %in% c("ughne", "unitghne")) {
    "unit-Generalized Half-Normal-E"
  }
  else if (family %in% c("ughnx", "unitghnx")) {
    "unit-Generalized Half-Normal-X"
  }
  else if (family %in% c("ugompertz", "unitgompertz")) {
    "unit-Gompertz"
  }
  else if (family %in% c("uburrxii", "unitburrxii")) {
    "unit-Burr-XII"
  }
  else if (family %in% c("johnsonsb", "jsb")) {
    "Johnson-SB"
  }
  else if (family %in% c("ashw", "arcsechwei")) {
    "arc-secant hyperbolic Weibull"
  }
  else if (family %in% c("ugumbel", "unitgumbel")) {
    "unit-Gumbel"
  }
  else {
    NULL
  }
}


# Get abbreviation for d,p,q,r --------------------------------------------

.get_abbrev <- function(family, fname = TRUE) {
  if (fname) family <- .get_family_name(family)
  out <- switch (family,
                 "unit-Weibull" = "uweibull",
                 "Kumaraswamy" = "kum",
                 "unit-Logistic" = "ulogistic",
                 "unit-Chen" = "uchen",
                 "unit-Birnbaum-Saunders" = "ubs",
                 "log-extended Exponential-Geometric" = "leeg",
                 "unit-Generalized Half-Normal-E" = "ughne",
                 "unit-Generalized Half-Normal-X" = "ughnx",
                 "unit-Gompertz" = "ugompertz",
                 "unit-Burr-XII" = "uburrxii",
                 "Johnson-SB" = "johnsonsb",
                 "arc-secant hyperbolic Weibull" = "ashw",
                 "unit-Gumbel" = "ugumbel"
                 )
  out
}

# Format output ------------------------------------------------------------
.format_perc <- function(probs, digits) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")
}

.FF <- function(x,Digits = 4, Width = 4){
  formatC(x, digits = Digits, width = Width, format = "f")
}



# Function to plot estimated coefficients versus quantile levels
#' @importFrom graphics par plot mtext abline lines title grid polygon points legend
#' @importFrom grDevices n2mfrow gray
#' @importFrom stats predict coef confint

.plot_coef <- function(x, output_df = FALSE, parm = NULL, level = 0.95,
                       mean_effect = FALSE, mfrow = NULL, mar = NULL, ylim = NULL,
                       main = NULL, col = gray(c(0, 0.75)), border = NULL, cex = 1,
                       pch = 20, type = "b", xlab = bquote("Quantile level ("*tau*")"),
                       ylab = "Estimate effect", ...) {

  # Get point and interval estimates
  lt_cfs <- lapply(x, function(u) {
    cf <- coef(u, type = "quantile")
    bd <- confint(u, level = level)
    bd <- bd[-grep("(theta)_*", rownames(bd)), ]
    cf <- cbind(cf, bd, u[["tau"]])
    colnames(cf) <- c("coef", "lower_bound", "upper_bound", "tau")
    cf
  })
  mat_cfs <- do.call("rbind", lt_cfs)

  if (mean_effect) {
    cf_mean <- tapply(mat_cfs[, "coef"], rownames(mat_cfs), mean)
    mat_mean_effect <- matrix(c(
      cf_mean,
      tapply(mat_cfs[, "lower_bound"], rownames(mat_cfs), mean),
      tapply(mat_cfs[, "upper_bound"], rownames(mat_cfs), mean)),
      ncol = 3
    )
    rownames(mat_mean_effect) <- names(cf_mean)
    colnames(mat_mean_effect) <- c("coef", "lower_bound", "upper_bound")
  }

  # Get parameters names
  if (is.null(parm))
    parm <- rownames(lt_cfs[[1L]])
  if (is.numeric(parm))
    parm <- rownames(lt_cfs[[1L]])[parm]
  parm <- parm[parm != "theta"]

  # Keeping user's graphs options
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Graphics options
  if (is.null(mfrow))
    mfrow <- n2mfrow(length(parm))
  if (is.null(mar))
    mar <- c(4, 4, 1, 1)
  if (is.null(border))
    border <- col[2]
  if (is.null(main))
    main <- parm
  main <- rep(main, length.out = length(parm))

  # Plot variables effects according to quantile
  par(mfrow = mfrow, mar = mar)
  for (i in seq_along(parm)) {
    df_cur <- mat_cfs[rownames(mat_cfs) == parm[i], ]
    yl <- range(c(df_cur[, "lower_bound"], df_cur[, "upper_bound"]))
    if (mean_effect)
      yl <- range(c(yl,
                    mat_mean_effect[rownames(mat_mean_effect) == parm[i], 2],
                    mat_mean_effect[rownames(mat_mean_effect) == parm[i], 3]))

    plot(x = df_cur[, "tau"], y = df_cur[, "coef"], type = "n", ylim = yl,
         xlab = xlab, ylab = ylab, main = main[i])
    polygon(x = c(df_cur[, "tau"], rev(df_cur[, "tau"])),
            y = c(df_cur[, "lower_bound"], rev(df_cur[, "upper_bound"])),
            border = TRUE, col = col[2L])
    points(x = df_cur[, "tau"], y = df_cur[, "coef"], cex = cex, pch = pch,
           type = type, col = col[1L], ...)
    abline(h = 0, col = gray(0.3), lty = 2)
    if (mean_effect) {
      abline(h = mat_mean_effect[rownames(mat_mean_effect) == parm[i], 1], col = "red",
             lty = 1)
      abline(h = mat_mean_effect[rownames(mat_mean_effect) == parm[i], 2], col = "red",
             lty = 2)
      abline(h = mat_mean_effect[rownames(mat_mean_effect) == parm[i], 3], col = "red",
             lty = 2)
    }
    grid()
  }

  if (output_df) {
    df2output <- as.data.frame(mat_cfs)
    df2output$parm <- rownames(mat_cfs)
    rownames(df2output) <- NULL
    return(df2output)
  }

  invisible()
}

# Function to plot conditional distribution at specific values of covariates
.plot_conddist <- function(x, dist_type = c("density", "cdf"), at_avg = TRUE,
                           at_obs = NULL, output_df = FALSE,
                           legend_position = "topleft", ...) {

  type <- match.arg(dist_type)
  distr_name <- .get_abbrev(x[[1L]]$family, fname = FALSE)

  # Get pdf or cdf
  cond_foo <- if (type == "cdf") match.fun(paste0("p", distr_name)) else match.fun(paste0("d", distr_name))

  # Create the data.frame with observed value of covariates
  if (!is.null(at_obs)) {
    newdata <- expand.grid(at_obs)
    newdata$avg <- FALSE
  } else newdata <- data.frame()
  if (at_avg) {
    df_avg <- as.data.frame(
        t(colMeans(model.matrix(x[[1L]]))))[, -1, drop = FALSE]
    df_avg$avg <- TRUE
    newdata <- rbind(newdata, df_avg)
  }
  covariates_names <- colnames(newdata)
  covariates_names <- covariates_names[covariates_names != "avg"]

  # Create identifier (id) for each value of covariate
  for (nm in covariates_names) {
    newdata[[paste0("id_", nm)]] <- paste0(nm, "=", round(newdata[[nm]], 2))
  }
  id_nm <- colnames(newdata)
  id_nm <- id_nm[grep(pattern = "id_", id_nm)]
  if (length(id_nm) > 1) {
    newdata[["id"]] <- apply(newdata[, id_nm], 1, paste, collapse = "; ")
  } else {
    newdata[["id"]] <- newdata[, id_nm]
  }
  newdata <- newdata[, -which(colnames(newdata) %in% id_nm)]

  # Get predict value for parameters (mu and theta)
  df_mu <- do.call(rbind, lapply(x, function(z) {
    tmp <- predict(z, type = "quantile", newdata = newdata)
    tmp$tau <- z$tau
    tmp
  }))
  if (x[[1L]]$theta_const) {
    df_theta <- do.call(rbind, lapply(x, function(z) {
      tmp <- newdata
      tmp$theta <- z$fitted.values$theta[1L]
      tmp$tau <- z$tau
      tmp
    }))
  } else {
    df_theta <- do.call(rbind, lapply(x, function(z) {
      tmp <- predict(z, type = "shape", newdata = newdata)
      tmp$tau <- z$tau
      tmp
    }))
  }
  rownames(df_mu) <- rownames(df_theta) <- NULL
  colnames(df_mu)[which(colnames(df_mu) == "fit")] <- "mu"
  colnames(df_theta)[which(colnames(df_theta) == "fit")] <- "theta"
  df_parms <- cbind(df_mu, theta = df_theta[, "theta"])

  # Compute cdf or pdf
  lt_out <- lapply(split(df_parms, df_parms$id), function(z) {
    z$predict <- sapply(seq_along(z$tau), function(i) {
      cond_foo(z$tau[i], mu = z$mu[i], theta = z$theta[i], tau = z$tau[i])
    })
    z
  })
  df2output <- do.call(rbind, lt_out)
  df2output$type <- type
  rownames(df2output) <- NULL

  # Plotting
  yl <- range(df2output$predict)
  xl <- range(df2output$tau)
  ylab <- ifelse(type == "cdf", "Cumulative density function",
                 "Probability density function")
  type_line <- ifelse(type == "cdf", "S", "l")
  plot(x = 0, y = 0, ylim = yl, xlim = xl, type = "n", ylab = ylab,
       xlab = "y", ...)
  counter <- 1
  for (id in unique(df2output$id)) {
    df_cur <- df2output[df2output$id == id, ]
    lines(x = df_cur$tau, y = df_cur$predict, col = counter, type = type_line)
    counter <- counter + 1L
  }
  legend(legend_position, col = seq_len(counter), lty = 1,
         legend = unique(df2output$id), bty = "n", ...)

  if (output_df) {
    return(df2output)
  }

}
