test_that("predict unitquantreg works", {

  # Theta fixed -----------------------------------------------------------------------
  set.seed(6669)
  n <- 200
  betas <- c(1, 2)
  X <- cbind(1, x1 = runif(n))
  eta <- drop(X %*% betas)
  mu <- exp(eta) / (1 + exp(eta))
  theta <- 2.0
  data_simulated <- data.frame(x1 = X[, 2])
  tau <- 0.5
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = tau))
    data_simulated$y[data_simulated$y == 0] <- 0.0001
    data_simulated$y[data_simulated$y == 1] <- 0.9999
    unitquantreg(formula = y ~ x1, tau = tau, data = data_simulated,
                 family = lt_families[[i]])
  })

  # Predict quantile

  ## Point predict
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "quantile"))
  colnames(tmp) <- names(lt_families)
  cbind(mu, tmp)[1:10, ]
  expect_equal(dim(tmp), c(n, length(lt_families)), label = "quantile predict")

  ## Point and standard error
  tmp <- lapply(lt_fits, predict, type = "quantile", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "quantile predict with se.fit")

  ## Point and confidence interval
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "quantile",
                                 interval = "confidence"))
  expect_equal(dim(tmp), c(n, length(lt_families) * 3),
               label = "quantile predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits, predict, type = "quantile", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), c(n, length(lt_families) * 3),
               label = "quantile predict with se.fit and confidence interval")
  expect_equal(dim(check_2), c(n, length(lt_families)),
               label = "quantile predict with se.fit and confidence interval")

  # Predict link

  ## Just point predict
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "link"))
  colnames(tmp) <- names(lt_families)
  cbind(eta, tmp)[1:10, ]
  expect_equal(dim(tmp), c(n, length(lt_families)), label = "link predict")

  ## Point and standard error
  tmp <- lapply(lt_fits, predict, type = "link", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "link predict with se.fit")

  ## Point and confidence
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "link",
                                 interval = "confidence"))
  expect_equal(dim(tmp), c(n, length(lt_families) * 3),
               label = "link predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits, predict, type = "link", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), c(n, length(lt_families) * 3),
               label = "link predict with se.fit and confidence interval")
  expect_equal(dim(check_2), c(n, length(lt_families)),
               label = "link predict with se.fit and confidence interval")

  # Predict theta
  expect_error(lapply(lt_fits, predict, type = "theta"),
               label = "try predict theta with constant shape model")

  # Predict terms
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "terms"))
  expect_equal(dim(tmp), c(n, length(lt_families)), label = "terms predict")



  # Theta varying ---------------------------------------------------------------------

  Z <- cbind(1, z1 = runif(n))
  gammas <- c(1, 0.2)
  theta <- exp(drop(Z %*% gammas))
  data_simulated$z1 <- Z[, 2]

  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], "\n")
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = tau))
    data_simulated$y[data_simulated$y == 0] <- 0.0001
    data_simulated$y[data_simulated$y == 1] <- 0.9999
    unitquantreg(formula = y ~ x1 | z1, tau = tau, data = data_simulated,
                 link.theta = "log", family = lt_families[[i]])
  })

  # Theta predict
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "shape"))
  colnames(tmp) <- names(lt_families)
  cbind(theta, tmp)[1:10, ]
  expect_equal(dim(tmp), c(n, length(lt_families)), label = "shape predict")

  ## Point and standard error
  tmp <- lapply(lt_fits, predict, type = "shape", se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), dim(check_2), label = "shape predict with se.fit")

  ## Point and confidence
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "shape",
                                 interval = "confidence"))
  expect_equal(dim(tmp), c(n, length(lt_families) * 3),
               label = "shape predict with confidence interval")

  ## Point, standard error and confidence
  tmp <- lapply(lt_fits, predict, type = "shape", interval = "confidence",
                se.fit = TRUE)
  check_1 <- do.call("cbind", lapply(tmp, "[[", "fit"))
  check_2 <- do.call("cbind", lapply(tmp, "[[", "se.fit"))
  expect_equal(dim(check_1), c(n, length(lt_families) * 3),
               label = "shape predict with se.fit and confidence interval")
  expect_equal(dim(check_2), c(n, length(lt_families)),
               label = "shape predict with se.fit and confidence interval")

  #############################################################################
  # Generating new covariate
  set.seed(1212)
  data_simulated$x2 <- rbinom(n, size = 1, prob = 0.5)
  lt_fits <- lapply(seq_along(lt_families), function(i) {
    cat(lt_families[[i]], "\n")
    rfun <- match.fun(paste0("r", lt_families[[i]]))
    data_simulated$y <- do.call(rfun, list(n, mu = mu, theta = theta, tau = tau))
    data_simulated$y[data_simulated$y == 0] <- 0.0001
    data_simulated$y[data_simulated$y == 1] <- 0.9999
    unitquantreg(formula = y ~ x1 + x2, tau = tau, data = data_simulated,
                 family = lt_families[[i]])
  })
  tmp <- do.call("cbind", lapply(lt_fits, predict, type = "terms"))
  expect_equal(dim(tmp), c(n, 2 * length(lt_families)),
               label = "terms predict with two covariate")

  # Check predict with newdata
  data2pred <- data.frame(x1 = runif(4), x2 = rbinom(n = 4, size = 1, prob = 0.5))
  m <- lt_fits[[1]]
  pred_1 <- predict(m, newdata = data2pred, type = "quantile")
  pred_2 <- predict(m, newdata = data2pred, type = "quantile",
                             interval = "confidence")
  pred_3 <- predict(m, newdata = data2pred, type = "quantile",
                             interval = "confidence", se.fit = TRUE)

  expect_equal(ncol(pred_1), ncol(data2pred) + 1,
               label = "predict with newdata")
  expect_equal(ncol(pred_2), ncol(data2pred) + 3,
               label = "predict with newdata and interval")
  expect_equal(length(pred_3), 2,
               label = "predict with se.fit and confidence interval")

})
