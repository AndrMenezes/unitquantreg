test_that("testing control options in the fitting", {

  skip_on_cran()

  data(sim_bounded, package = "unitquantreg")

  # Gradient argument
  lt_fits <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]

    fit_nogradient <- unitquantreg(
      formula = y1 ~ x, data = sim_bounded_curr,
      tau = 0.5, family = fam, link.theta = "log",
      control = unitquantreg.control(gradient = FALSE))
    fit_gradient <- unitquantreg(
      formula = y1 ~ x, data = sim_bounded_curr,
      tau = 0.5, family = fam, link.theta = "log",
      control = unitquantreg.control(gradient = TRUE))

    se_nogradient <- sqrt(diag(vcov(fit_nogradient)))
    se_gradient <- sqrt(diag(vcov(fit_gradient)))
    coef_nogradient <- coef(fit_nogradient)
    coef_gradient <- coef(fit_gradient)

    df_se <- data.frame(se_nogradient, se_gradient, coef_nogradient,
                        coef_gradient)
    df_se <- data.frame(apply(df_se, 2, function(x) round(x, 6)))
    df_se$family <- fam
    df_se$parms <- rownames(df_se)
    rownames(df_se) <- NULL
    list(
      elapsed_time = c(nogradient = fit_nogradient$elapsed_time,
                       gradient = fit_gradient$elapsed_time),
      df_se = df_se
    )
  })
  df_res_parms <- do.call("rbind", lapply(lt_fits, "[[", "df_se"))
  npar <- length(unique(df_res_parms$parms))
  expect_equal(nrow(df_res_parms), length(lt_families) * npar,
               label = "gradient argument in control works")


  # Hessian argument
  lt_fits_2 <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]

    fit_numeric <- unitquantreg(
      formula = y1 ~ x, data = sim_bounded_curr,
      tau = 0.5, family = fam, link.theta = "log",
      control = unitquantreg.control(hessian = TRUE))

    fit_analitic <- unitquantreg(
      formula = y1 ~ x, data = sim_bounded_curr,
      tau = 0.5, family = fam, link.theta = "log",
      control = unitquantreg.control(hessian = FALSE))

    se_numeric <- round(sqrt(diag(vcov(fit_numeric))), 5)
    se_analitic <- round(sqrt(diag(vcov(fit_analitic))), 5)
    df_se <- data.frame(se_numeric, se_analitic)
    df_se$family <- fam
    list(
      elapsed_time = c(numeric = fit_numeric$elapsed_time,
                       analitic = fit_analitic$elapsed_time),
      df_se = df_se
    )
  })
  expect_equal(length(lt_fits_2), length(lt_families))
  expect_equal(length(lt_fits_2[[1L]]), 2)
  expect_equal(class(lt_fits_2[[1L]]$df_se), "data.frame")



})

