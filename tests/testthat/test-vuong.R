test_that("vuong and pairwise vuong tests works", {
  data(water)
  lt_fits <- lapply(lt_families, function(fam) {
    cat(fam, "\n")
    unitquantreg(formula = phpws ~ mhdi,
                 tau = 0.5, data = water, family = fam)
  })
  ans <- pairwise.vuong.test(lt = lt_fits)
  ans
})
