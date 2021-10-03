test_that("residuals method works for all family of distributions", {

  data(water)
  lt_fits <- lapply(lt_families, function(fam){
    cat(fam, "\n")
    unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                 tau = 0.5, data = water, family = fam)
  })
  sort(sapply(lt_fits, AIC))

  lt_qres <- lapply(lt_fits, residuals, type = "quantile")
  # x11()
  # par(mfrow = c(2, 4))
  # invisible(lapply(seq_along(lt_qres), function(i){
  #   qqnorm(lt_qres[[i]][is.finite(lt_qres[[i]])], main = lt_fits[[i]]$family)
  # }))

  expect_equal(length(lt_qres), length(lt_families),
               label = "quantile residuals works for all family of distributions")


  lt_qcs <- lapply(lt_fits, residuals, type = "cox-snell")
  # x11()
  # par(mfrow = c(2, 4))
  # invisible(lapply(seq_along(lt_qres), function(i){
  #   teo <- lt_qcs[[i]][is.finite(lt_qcs[[i]])]
  #   n <- length(teo)
  #   qqplot(x = qexp((1:n - 3 / 8) / (n + 1 / 4)),
  #          y = teo, main = lt_fits[[i]]$family)
  # }))

  expect_equal(length(lt_qcs), length(lt_families),
               label = "cox-snell residuals works for all family of distributions")

})


test_that("working residuals to check link function", {

  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1 = runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 2.0
  data_simulated <- data.frame(x1 = X[, 2])
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = 0.5))
    data_simulated$y[data_simulated$y == 0] <- 0.00001
    unitquantreg(formula = y ~ x1, tau = 0.5, data = data_simulated,
                 family = lt_families[[i]])
  })
  names(lt_fits) <- names(lt_families)
  t(sapply(lt_fits, coef))

  lt_working <- lapply(lt_fits, residuals, type = "working")
  expect_equal(length(lt_working), length(lt_families),
               label = "working residuals works for all family of distributions")
  # x11()
  # par(mfrow = c(2, 6))
  # invisible(lapply(lt_fits, function(x) {
  #   eta <- x$linear.predictors$mu
  #   z <- eta + residuals(x, type = "working")
  #   plot(z ~ eta, main = x$family)
  #   lines(lowess(z ~ eta), lwd = 2, col = "blue")
  # }))

  lt_partial <- lapply(lt_fits, residuals, type = "partial")
  expect_equal(length(lt_partial), length(lt_families),
               label = "partial residuals works for all family of distributions")

})
