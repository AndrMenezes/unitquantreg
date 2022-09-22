test_that("likelihood_stats works for bodyfat data", {
  data(bodyfat)
  lt_fits <- lapply(lt_families, function(fam) {
    cat(fam, "\n")
    unitquantreg(formula = arms ~ bmi, tau = 0.5, data = bodyfat,
                 family = fam)
  })
  tmp <- likelihood_stats(lt = lt_fits)
  expect_equal(class(tmp), "likelihood_stats")
  expect_equal(length(tmp), 2)
})
