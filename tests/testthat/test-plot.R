test_that("plot method for unitquantreg(s) class works", {

  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1 = runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 1.0
  data_simulated <- data.frame(x1 = X[, 2])
  data_simulated$y <- ruweibull(n, mu = mu, theta = theta, tau = 0.5)

  # For one quantile (tau)
  fit <- unitquantreg(formula = y ~ x1 | x1, tau = 0.5, data = data_simulated,
                      family = "uweibull")
  plot(fit, which = 1)
  plot(fit, which = 2)
  plot(fit, which = 3)
  plot(fit, which = 4)
  par(mfrow = c(2, 2))
  plot(fit)

  # For several quantiles (taus)
  fits <- unitquantreg(formula = y ~ x1 + I(x1^2) | x1, tau = 1:49 / 50, data = data_simulated,
                       family = "uweibull")
  plot(fits, which = "coef")
  plot(fits, which = "coef", mean_effect = TRUE)
  plot(fits, which = "coef", parm = "x1")
  plot(fits, which = "coef", parm = "x1", mean_effect = FALSE)

  plot(fits, which = "conddist", at_obs = list(x1 = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x1^2)" = c(0.4, 0.8)),
       dist_type = "cdf")
  plot(fits, which = "conddist", at_obs = list(x1 = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x1^2)" = c(0.4, 0.8)),
       dist_type = "density")

  data(water)
  fits <- unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                       data = water, tau = 1:49/50,
                       family = "uweibull", link = "logit",
                       link.theta = "log")
  plot(fits, which = "conddist",
       at_obs = list(mhdi = c(0.5, 0.7), incpc = mean(water$incpc),
                     region = c(1, 0), pop = mean(water$pop)),
       at_avg = FALSE, dist_type = "density")
  plot(fits, which = "conddist",
       at_obs = list(mhdi = c(0.5, 0.7), incpc = mean(water$incpc),
                     region = c(1, 0), pop = mean(water$pop)),
       at_avg = FALSE, dist_type = "cdf")

})
