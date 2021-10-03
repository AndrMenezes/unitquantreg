test_that("analitical and numerical hessian works in water data", {

  data(water)

  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    fit_numeric <- unitquantreg(
      formula = phpws ~ mhdi + incpc + region + log(pop),
      tau = 0.5, data = water, family = lt_families[[i]], link.theta = "log",
      control = unitquantreg.control(hessian = TRUE))
    fit_analitic <- unitquantreg(
      formula = phpws ~ mhdi + incpc + region + log(pop),
      tau = 0.5, data = water, family = lt_families[[i]], link.theta = "log",
      control = unitquantreg.control(hessian = FALSE))

    se_numeric <- round(sqrt(diag(vcov(fit_numeric))), 5)
    se_analitic <- round(sqrt(diag(vcov(fit_analitic))), 5)
    df_se <- data.frame(se_numeric, se_analitic)
    df_se$family <- lt_families[[i]]
    list(
      elapsed_time = c(numeric = fit_numeric$elapsed_time,
                       analitic = fit_analitic$elapsed_time),
      df_se = df_se
    )
  })

  expect_equal(length(lt_fits), length(lt_families))
  expect_equal(length(lt_fits[[1L]]), 2)
  expect_equal(class(lt_fits[[1L]]$df_se), "data.frame")

})
