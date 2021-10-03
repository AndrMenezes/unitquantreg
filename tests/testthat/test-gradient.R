test_that("analitical and numerical gradient works in water data", {

  data(water)

  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    fit_nogradient <- unitquantreg(
      phpws ~ mhdi + incpc + region + log(pop), data = water,
      tau = 0.5, family = lt_families[[i]], link.theta = "log",
      control = unitquantreg.control(gradient = FALSE))
    fit_gradient <- unitquantreg(
      phpws ~ mhdi + incpc + region + log(pop), data = water,
      tau = 0.5, family = lt_families[[i]], link.theta = "log",
      control = unitquantreg.control(gradient = TRUE))
    se_nogradient <-sqrt(diag(vcov(fit_nogradient)))
    se_gradient <- sqrt(diag(vcov(fit_gradient)))
    coef_nogradient <- coef(fit_nogradient)
    coef_gradient <- coef(fit_gradient)

    df_se <- data.frame(se_nogradient, se_gradient, coef_nogradient, coef_gradient)
    df_se <- data.frame(apply(df_se, 2, function(x) round(x, 6)))
    df_se$family <- lt_families[[i]]
    df_se$parms <- rownames(df_se)
    rownames(df_se) <- NULL
    list(
      elapsed_time = c(nogradient = fit_nogradient$elapsed_time,
                       gradient = fit_gradient$elapsed_time),
      df_se = df_se
    )
  })
  df_res_parms <- do.call("rbind", lapply(lt_fits, "[[", "df_se"))
  npar <- 6
  expect_equal(ncol(df_res_parms), 6)
  expect_equal(nrow(df_res_parms), 11 * npar)

})

test_that("analitical and numerical gradient works in bodyfat data", {

  data(bodyfat)
  head(bodyfat)
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], '\n')
    fit_nogradient <- unitquantreg(arms ~ age + sex + ipaq, data = bodyfat, tau = 0.5,
                                   family = lt_families[[i]], link.theta = "log",
                                   control = unitquantreg.control(gradient = FALSE))
    fit_gradient <- unitquantreg(arms ~ age + sex + ipaq, data = bodyfat, tau = 0.5,
                                 family = lt_families[[i]], link.theta = "log",
                                 control = unitquantreg.control(gradient = TRUE))
    se_nogradient <-sqrt(diag(vcov(fit_nogradient)))
    se_gradient <- sqrt(diag(vcov(fit_gradient)))
    coef_nogradient <- coef(fit_nogradient)
    coef_gradient <- coef(fit_gradient)

    df_se <- data.frame(se_nogradient, se_gradient, coef_nogradient, coef_gradient)
    df_se <- data.frame(apply(df_se, 2, function(x) round(x, 6)))
    df_se$family <- lt_families[[i]]
    df_se$parms <- rownames(df_se)
    rownames(df_se) <- NULL
    list(
      elapsed_time = c(nogradient = fit_nogradient$elapsed_time,
                       gradient = fit_gradient$elapsed_time),
      df_se = df_se
    )
  })
  df_res_parms <- do.call("rbind", lapply(lt_fits, "[[", "df_se"))
  npar <- 6
  expect_equal(ncol(df_res_parms), 6)
  expect_equal(nrow(df_res_parms), 11 * npar)

})

