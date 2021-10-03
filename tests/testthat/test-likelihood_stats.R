test_that("likelihood_stats works for water data", {
  data(water)
  system.time(
    lt_fits <- lapply(lt_families, function(fam) {
      cat(fam, "\n")
      unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                   tau = 0.5, data = water, family = fam)
  }))

  tmp <- likelihood_stats(lt = lt_fits)
  expect_equal(class(tmp), "likelihood_stats")
  expect_equal(length(tmp), 2)

})

test_that("likelihood_stats works for bodyfat data", {
  data(bodyfat)
  system.time(
    lt_fits <- lapply(lt_families, function(fam) {
      cat(fam, "\n")
      unitquantreg(formula = arms ~ age + sex + ipaq, tau = 0.5, data = bodyfat,
                   family = fam)
  }))

  tmp <- likelihood_stats(lt = lt_fits)
  expect_equal(class(tmp), "likelihood_stats")
  expect_equal(length(tmp), 2)

})
