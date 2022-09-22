test_that("vuong and pairwise vuong tests works", {
  data(bodyfat)
  lt_fits <- lapply(lt_families, function(fam) {
    cat(fam, "\n")
    unitquantreg(formula = arms ~ bmi, tau = 0.5, data = bodyfat, family = fam)
  })
  ans <- pairwise.vuong.test(lt = lt_fits)
  ans
})
