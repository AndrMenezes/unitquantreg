test_that("methods works for unitquanreg objects", {

  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  gamas <- c(-1, 1)
  X <- cbind(1, x1 = runif(n))
  Z <- cbind(1, z1 = rexp(n))
  eta <- drop(X %*% betas)
  zeta <- drop(Z %*% gamas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- exp(zeta)
  data_simulated <- data.frame(x1 = X[, 2], z1 = Z[, 2])
  data_simulated$y <- ruweibull(n, mu, theta)

  fit_1 <- unitquantreg(formula = y ~ x1 + z1 + I(x1^2) | z1 + x1,
                        data = data_simulated, family = "uweibull",
                        tau = 0.5, link.theta = "log")
  fit_1
  summary(fit_1)
  vcov(fit_1)
  coef(fit_1)
  confint(fit_1)
  terms(fit_1)
  model.frame(fit_1)
  model.matrix(fit_1)
  update(fit_1, . ~ . -x1)
  update(fit_1, . ~ . -z1)
  update(fit_1, . ~ . -I(x1^2))
  update(fit_1, . ~ . | . -z1)
  update(fit_1, . ~ . | . -x1)
})
