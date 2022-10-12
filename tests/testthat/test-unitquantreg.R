test_that("testing fit and their methods with simulated data for theta fixed", {

  # Load data
  data(sim_bounded, package = "unitquantreg")

  # True values
  betas <- c(1, 2)
  theta <- 2

  # Fitting
  lt_fits_1 <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]
    unitquantreg(formula = y1 ~ x, data = sim_bounded_curr, tau = 0.5,
                 family = fam)
  })
  names(lt_fits_1) <- names(lt_families)

  expect_equal(length(lt_fits_1), length(lt_families),
               label = "fit procedure works for all families of distributions")

  m_parms <- t(sapply(lt_fits_1, coef))
  expect_equal(unname(m_parms[, 1]), rep(betas[1], length(lt_families)), tol = 2e-1,
               label = "estimated beta_0 is equal to true beta_0 with 0.2 of tolerance")

  expect_equal(unname(m_parms[, 2]), rep(betas[2], length(lt_families)), tol = 2e-1,
               label = "estimated beta_1 is equal to true beta_1 with 0.2 of tolerance")

  expect_equal(unname(m_parms[, 3]), rep(theta, length(lt_families)), tol = 2e-1,
               label = "estimated theta is equal to true theta with 0.2 of tolerance")

  # Residuals
  lt_qres <- lapply(lt_fits_1, residuals, type = "quantile")
  expect_equal(length(lt_qres), length(lt_families),
               label = "quantile residuals works for all family of distributions")

  lt_qcs <- lapply(lt_fits_1, residuals, type = "cox-snell")
  expect_equal(length(lt_qcs), length(lt_families),
               label = "cox-snell residuals works for all family of distributions")

  lt_working <- lapply(lt_fits_1, residuals, type = "working")
  expect_equal(length(lt_working), length(lt_families),
               label = "working residuals works for all family of distributions")

  lt_partial <- lapply(lt_fits_1, residuals, type = "partial")
  expect_equal(length(lt_partial), length(lt_families),
               label = "partial residuals works for all family of distributions")

  # Likelihood stats
  like_stats <- likelihood_stats(lt = lt_fits_1)
  expect_equal(class(like_stats), "likelihood_stats")

  # Voung
  tmp_vuong <- vuong.test(lt_fits_1[[1L]], lt_fits_1[[2L]])
  expect_s3_class(tmp_vuong, "htest")
  tmp_pairwise_vuong <- pairwise.vuong.test(lt = lt_fits_1)
  expect_s3_class(tmp_pairwise_vuong, "pairwise.htest")

  # Predict quantile

  ## Point predict
  tmp <- do.call("cbind", lapply(lt_fits_1, predict, type = "quantile"))
  colnames(tmp) <- names(lt_families)
  expect_equal(ncol(tmp), length(lt_families), label = "quantile predict")

  ## Point and standard error
  tmp <- lapply(lt_fits_1, predict, type = "quantile", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "quantile predict with se.fit")

  ## Point and confidence interval
  tmp <- do.call("cbind", lapply(lt_fits_1, predict, type = "quantile",
                                 interval = "confidence"))
  expect_equal(ncol(tmp), length(lt_families) * 3,
               label = "quantile predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits_1, predict, type = "quantile", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(ncol(check_1), length(lt_families) * 3,
               label = "quantile predict with se.fit and confidence interval")
  expect_equal(ncol(check_2), length(lt_families),
               label = "quantile predict with se.fit and confidence interval")

  # Predict link

  ## Just point predict
  tmp <- do.call("cbind", lapply(lt_fits_1, predict, type = "link"))
  colnames(tmp) <- names(lt_families)
  expect_equal(ncol(tmp), length(lt_families), label = "link predict")

  ## Point and standard error
  tmp <- lapply(lt_fits_1, predict, type = "link", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "link predict with se.fit")

  ## Point and confidence
  tmp <- do.call("cbind", lapply(lt_fits_1, predict, type = "link",
                                 interval = "confidence"))
  expect_equal(ncol(tmp), length(lt_families) * 3,
               label = "link predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits_1, predict, type = "link", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(ncol(check_1), length(lt_families) * 3,
               label = "link predict with se.fit and confidence interval")
  expect_equal(ncol(check_2), length(lt_families),
               label = "link predict with se.fit and confidence interval")

  # Predict theta
  expect_error(lapply(lt_fits_1, predict, type = "theta"),
               label = "try predict theta with constant shape model")

  # Predict terms
  tmp <- do.call("cbind", lapply(lt_fits_1, predict, type = "terms"))
  expect_equal(ncol(tmp), length(lt_families), label = "terms predict")


})

test_that("testing fit and related methods with simulated data for theta varying", {

  # Load data
  data(sim_bounded, package = "unitquantreg")

  # True values
  betas <- c(1, 2)
  gammas <- c(-1, 1)

  # Fitting
  lt_fits_2 <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]
    unitquantreg(formula = y2 ~ x | z, data = sim_bounded_curr,
                 tau = 0.5, family = fam, link.theta = "log")
  })
  names(lt_fits_2) <- names(lt_families)

  m_parms <- t(sapply(lt_fits_2, coef))
  cbind(m_parms[, 1], betas[1], m_parms[, 1] - betas[1])
  expect_equal(unname(m_parms[, 1]), rep(betas[1], length(lt_families)), tol = 6e-1,
               label = "estimated beta_0 is equal to true beta_0 with 0.6 of tolerance")

  expect_equal(unname(m_parms[, 2]), rep(betas[2], length(lt_families)), tol = 6e-1,
               label = "estimated beta_1 is equal to true beta_1 with 0.6 of tolerance")

  expect_equal(unname(m_parms[, 3]), rep(gammas[1], length(lt_families)), tol = 6e-1,
               label = "estimated gamma_0 is equal to true gamma_0 with 0.6 of tolerance")

  expect_equal(unname(m_parms[, 4]), rep(gammas[2], length(lt_families)), tol = 6e-1,
               label = "estimated gamma_1 is equal to true gamma_1 with 0.6 of tolerance")

  # Residuals
  lt_qres <- lapply(lt_fits_2, residuals, type = "quantile")
  expect_equal(length(lt_qres), length(lt_families),
               label = "quantile residuals works for all family of distributions")

  lt_qcs <- lapply(lt_fits_2, residuals, type = "cox-snell")
  expect_equal(length(lt_qcs), length(lt_families),
               label = "cox-snell residuals works for all family of distributions")

  lt_working <- lapply(lt_fits_2, residuals, type = "working")
  expect_equal(length(lt_working), length(lt_families),
               label = "working residuals works for all family of distributions")

  lt_partial <- lapply(lt_fits_2, residuals, type = "partial")
  expect_equal(length(lt_partial), length(lt_families),
               label = "partial residuals works for all family of distributions")

  # Likelihood stats
  like_stats <- likelihood_stats(lt = lt_fits_2)
  expect_equal(class(like_stats), "likelihood_stats")

  # Voung
  tmp_vuong <- vuong.test(lt_fits_2[[1L]], lt_fits_2[[2L]])
  expect_s3_class(tmp_vuong, "htest")
  tmp_pairwise_vuong <- pairwise.vuong.test(lt = lt_fits_2)
  expect_s3_class(tmp_pairwise_vuong, "pairwise.htest")

  # Predict
  tmp <- do.call("cbind", lapply(lt_fits_2, predict, type = "shape"))
  colnames(tmp) <- names(lt_families)
  expect_equal(ncol(tmp), length(lt_families), label = "shape predict")

  ## Point and standard error
  tmp <- lapply(lt_fits_2, predict, type = "shape", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "shape predict with se.fit")

  ## Point and confidence
  tmp <- do.call("cbind", lapply(lt_fits_2, predict, type = "shape",
                                 interval = "confidence"))
  expect_equal(ncol(tmp), length(lt_families) * 3,
               label = "shape predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits_2, predict, type = "shape", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(ncol(check_1), length(lt_families) * 3,
               label = "shape predict with se.fit and confidence interval")
  expect_equal(ncol(check_2), length(lt_families),
               label = "shape predict with se.fit and confidence interval")


})
