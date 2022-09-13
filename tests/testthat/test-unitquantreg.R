test_that("unitquanreg fits with simulated data works", {

  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1 = runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 2.0
  data_simulated <- data.frame(x1 = X[, 2])

  # Theta fixed
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = 0.5))
    data_simulated$y[data_simulated$y == 0] <- 0.00001
    unitquantreg(formula = y ~ x1, tau = 0.5, data = data_simulated,
                 family = lt_families[[i]])
  })
  names(lt_fits) <- names(lt_families)

  expect_equal(length(lt_fits), length(lt_families))

  m_parms <- t(sapply(lt_fits, coef))
  expect_equal(unname(m_parms[, 1]), rep(betas[1], length(lt_families)), tol = 2e-1,
               label = "estimated beta_0 is equal to true beta_0 with 0.2 of tolerance")

  expect_equal(unname(m_parms[, 2]), rep(betas[2], length(lt_families)), tol = 2e-1,
               label = "estimated beta_1 is equal to true beta_1 with 0.2 of tolerance")

  expect_equal(unname(m_parms[, 3]), rep(theta, length(lt_families)), tol = 2e-1,
               label = "estimated theta is equal to true theta with 0.2 of tolerance")


  # Theta varying
  set.seed(1212)
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

  m_parms <- t(sapply(lt_fits, coef))
  expect_equal(unname(m_parms[, 1]), rep(betas[1], length(lt_families)), tol = 1e-1,
               label = "estimated beta_0 is equal to true beta_0 with 0.1 of tolerance")

  expect_equal(unname(m_parms[, 2]), rep(betas[2], length(lt_families)), tol = 1e-1,
               label = "estimated beta_1 is equal to true beta_1 with 0.1 of tolerance")

  expect_equal(unname(m_parms[, 3]), rep(gammas[1], length(lt_families)), tol = 3e-1,
               label = "estimated gamma_0 is equal to true gamma_0 with 0.3 of tolerance")

  expect_equal(unname(m_parms[, 4]), rep(gammas[2], length(lt_families)), tol = 6e-1,
               label = "estimated gamma_1 is equal to true gamma_1 with 0.6 of tolerance")


})


test_that("comparison unitquantreg fits and SAS NLMIXED", {
  data(water)
  fit <- unitquantreg(formula = phpws ~ mhdi + region + incpc + log(pop),
                      tau = 0.5, data = water, family = "uweibull")

  SAS_coef <- c(-6.5145, 11.8262, -0.2699, 0.00026, 0.1055, 1.1883)
  expect_equal(unname(coef(fit)), SAS_coef, tol = 1e-5)

  SAS_se <- c(0.3272, 0.5562, 0.04944, 0.000149, 0.0155, 0.0142)
  R_se <- unname(sqrt(diag(vcov(fit))))
  expect_equal(R_se, SAS_se, tol = 1e-1)

  SAS_grad <- c(0.000892, 0.000701, -0.00105, 0.905682, 0.008567, -0.00087)
  cbind(round(fit$gradient, 6), SAS_grad)
  expect_equal(unname(fit$gradient), SAS_grad, tol = 9e-1)

})
