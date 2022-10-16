test_that("methods works for unitquanreg objects", {

  data(sim_bounded, package = "unitquantreg")
  sim_bounded_curr <- sim_bounded[sim_bounded$family == "uweibull", ]
  fit_1 <- unitquantreg(formula = y1 ~ x + z + I(x^2) | z + x,
                        data = sim_bounded_curr,
                        family = "uweibull",
                        tau = 0.5, link.theta = "log")
  fit_1
  expect_s3_class(summary(fit_1), "summary.unitquantreg")
  expect_type(vcov(fit_1), "double")
  expect_type(coef(fit_1), "double")
  expect_type(confint(fit_1), "double")
  expect_s3_class(terms(fit_1), "terms")
  expect_s3_class(model.frame(fit_1), "data.frame")
  expect_type(model.matrix(fit_1), "double")
  expect_s3_class(update(fit_1, . ~ . -x), "unitquantreg")
  expect_s3_class(update(fit_1, . ~ . -z), "unitquantreg")
  expect_s3_class(update(fit_1, . ~ . -I(x^2)), "unitquantreg")
  expect_s3_class(update(fit_1, . ~ . | . -z), "unitquantreg")
  expect_s3_class(update(fit_1, . ~ . | . -x), "unitquantreg")
})
