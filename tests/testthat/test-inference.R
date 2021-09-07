test_that("lmtest wrapper functions for inference works", {

  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1 = runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 2.0
  tau <- 0.5
  data_simulated <- data.frame(x1 = X[, 2])
  data_simulated$y <- ruweibull(n, mu = mu, theta = theta, tau = tau)

  fit1 <- unitquantreg(formula = y ~ x1, tau = tau, data = data_simulated,
                       family = "uweibull")
  fit2 <- unitquantreg(formula = y ~ 1, tau = tau, data = data_simulated,
                       family = "uweibull")

  tmp <- lmtest::coeftest(fit1)
  expect_equal(length(tmp), 12)

  tmp <- lmtest::coefci(fit1)
  expect_equal(length(tmp), 6)

  tmp <- lmtest::waldtest(fit1, fit2)
  expect_equal(length(tmp), 4)

  tmp <- lmtest::lrtest(fit1, fit2)
  expect_equal(length(tmp), 5)

})
