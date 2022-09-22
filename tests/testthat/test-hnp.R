test_that("hnp method works for all distribution families in simulated data", {

  set.seed(1212)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1=runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 2.0
  data_simulated <- data.frame(x1 = X[,2])

  # Theta fixed
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = 0.5))
    data_simulated$y[data_simulated$y == 0] <- 0.0001
    unitquantreg(formula = y ~ x1, tau = 0.5, data = data_simulated,
                 family = lt_families[[i]])
  })
  names(lt_fits) <- names(lt_families)
  # t(sapply(lt_fits, coef))
  # lapply(lt_fits, summary)

  # x11()
  # par(mfrow = c(2, 6))
  set.seed(6969)
  system.time(
    out_hnp <- invisible(lapply(seq_along(lt_fits), function(i) {
      cat(lt_families[[i]], '\n')
      hnp(lt_fits[[i]], main = names(lt_fits)[i], plot = FALSE, nsim = 10)
    }))
  )
  expect_equal(length(out_hnp), length(lt_families))

  # x11()
  # invisible(lapply(seq_along(lt_fits), function(i) {
  #   hnp(lt_fits[[i]], main = names(lt_fits)[i], halfnormal = FALSE)
  # }))

  # Theta varying
  set.seed(6969)
  Z <- cbind(1, z1 = rexp(n))
  gammas <- c(1, 0.5)
  theta <- exp(c(Z %*% gammas))
  data_simulated$z1 <- Z[, 2]

  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = 0.5))
    data_simulated$y[data_simulated$y == 0] <- 0.00001
    data_simulated$y[data_simulated$y == 1] <- 0.99999
    unitquantreg(formula = y ~ x1 | z1, tau = 0.5, data = data_simulated,
                 link.theta = "log", family = lt_families[[i]])
  })
  names(lt_fits) <- names(lt_families)
  # t(sapply(lt_fits, coef))
  # lapply(lt_fits, confint)
  # lapply(lt_fits, summary)

  # x11()
  # par(mfrow = c(2, 6))
  set.seed(6969)
  system.time(
    out_hnp <- invisible(lapply(seq_along(lt_fits), function(i) {
      cat(lt_families[[i]], '\n')
      hnp(lt_fits[[i]], main = names(lt_fits)[i], plot = FALSE, nsim = 10)
    }))
  )
  expect_equal(length(out_hnp), length(lt_families))

})
