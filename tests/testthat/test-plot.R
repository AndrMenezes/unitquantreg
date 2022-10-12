test_that("plot method for unitquantreg(s) class works", {

  skip_on_cran()

  data(sim_bounded, package = "unitquantreg")
  sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]

  # For one quantile (tau)
  fit <- unitquantreg(formula = y1 ~ x,
                      data = sim_bounded_curr,
                      family = "uweibull",
                      tau = 0.5, link.theta = "log")

  plot(fit, which = 1)
  plot(fit, which = 2)
  plot(fit, which = 3)
  plot(fit, which = 4, nsim = 5)

  # For several quantiles (taus)
  fits <- unitquantreg(formula = y1 ~ x + I(x^2) | x, tau = 1:49 / 50,
                       data = sim_bounded_curr, family = "uweibull",
                       link.theta = "log")
  plot(fits, which = "coef")
  plot(fits, which = "coef", mean_effect = TRUE)
  plot(fits, which = "coef", parm = "x")
  plot(fits, which = "coef", parm = "x", mean_effect = FALSE)

  plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x^2)" = c(0.4, 0.8)),
       dist_type = "cdf")
  plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x^2)" = c(0.4, 0.8)),
       dist_type = "density")
})
