test_that("hnp method works for all distribution families in simulated data", {

  skip_on_cran()

  lt_fits_1 <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]
    unitquantreg(formula = y1 ~ x, data = sim_bounded_curr, tau = 0.5,
                 family = fam)
  })
  names(lt_fits_1) <- names(lt_families)

  # x11()
  # par(mfrow = c(2, 6))
  set.seed(6969)
  system.time(
    out_hnp <- invisible(lapply(seq_along(lt_fits_1), function(i) {
      # cat(lt_families[[i]], '\n')
      hnp(lt_fits_1[[i]], main = names(lt_fits_1)[i], plot = FALSE, nsim = 10L)
    }))
  )
  expect_equal(length(out_hnp), length(lt_families))

  # x11()
  # invisible(lapply(seq_along(lt_fits), function(i) {
  #   hnp(lt_fits[[i]], main = names(lt_fits)[i], halfnormal = FALSE)
  # }))

  # Theta varying
  lt_fits_2 <- lapply(lt_families, function(fam) {
    # cat(fam, '\n')
    sim_bounded_curr <- sim_bounded[sim_bounded$family == fam, ]
    unitquantreg(formula = y2 ~ x | z, data = sim_bounded_curr,
                 tau = 0.5, family = fam, link.theta = "log")
  })
  names(lt_fits_2) <- names(lt_families)

  # x11()
  # par(mfrow = c(2, 6))
  set.seed(6969)
  system.time(
    out_hnp <- invisible(lapply(seq_along(lt_fits_2), function(i) {
      # cat(lt_families[[i]], '\n')
      hnp(lt_fits_2[[i]], main = names(lt_fits_2)[i], plot = FALSE, nsim = 10L)
    }))
  )
  expect_equal(length(out_hnp), length(lt_families))

})
