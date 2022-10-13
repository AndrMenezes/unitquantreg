test_that("lmtest wrapper functions for inference works", {

  data(sim_bounded, package = "unitquantreg")
  sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]

  fit_1 <- unitquantreg(formula = y1 ~ x,
                        data = sim_bounded_curr,
                        family = "uweibull",
                        tau = 0.5, link.theta = "log")
  fit_2 <- unitquantreg(formula = y1 ~ 1,
                        data = sim_bounded_curr,
                        family = "uweibull",
                        tau = 0.5, link.theta = "log")

  tmp <- lmtest::coeftest(fit_1)
  expect_equal(length(tmp), 12)

  tmp <- lmtest::coefci(fit_1)
  expect_equal(length(tmp), 6)

  tmp <- lmtest::waldtest(fit_1, fit_2)
  expect_equal(length(tmp), 4)

  tmp <- lmtest::lrtest(fit_1, fit_2)
  expect_equal(length(tmp), 5)

})
