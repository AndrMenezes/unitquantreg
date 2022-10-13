test_that("plot method for unitquantreg(s) class works", {

  skip_on_cran()

  data(sim_bounded, package = "unitquantreg")
  sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]

  # For one quantile (tau)
  fit <- unitquantreg(formula = y1 ~ x,
                      data = sim_bounded_curr,
                      family = "uweibull",
                      tau = 0.5, link.theta = "log")

  expect_no_error(plot(fit, which = 1))
  expect_no_error(plot(fit, which = 2))
  expect_no_error(plot(fit, which = 3))
  expect_no_error(plot(fit, which = 4, nsim = 5))


  # For several quantiles (taus)
  fits <- unitquantreg(formula = y1 ~ x + I(x^2), tau = 1:49 / 50,
                       data = sim_bounded_curr, family = "uweibull",
                       link.theta = "log")
  expect_no_error(plot(fits, which = "coef"))
  expect_no_error(plot(fits, which = "coef", mean_effect = TRUE))
  expect_no_error(plot(fits, which = "coef", parm = "x"))
  expect_no_error(plot(fits, which = "coef", parm = "x", mean_effect = FALSE))

  expect_error(
    plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7)),
         dist_type = "cdf"))
  expect_no_error(
    plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x^2)" = c(0.4, 0.8)),
         dist_type = "cdf"))
  expect_no_error(
    plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7),
                                               "I(x^2)" = c(0.4, 0.8)),
         dist_type = "density"))

  expect_type(
    plot(fits, which = "conddist", at_obs = list(x = c(0.2, 0.4, 0.6, 0.7),
                                                 "I(x^2)" = c(0.4, 0.8)),
         dist_type = "density", output_df = TRUE),
    "list")

})
