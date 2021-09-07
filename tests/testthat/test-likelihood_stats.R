lt_families <- list("unit-Weibull" = "uweibull",
                    "Kumaraswamy" = "kum",
                    "unit-Logistic" = "ulogistic",
                    "unit-Birnbaum-Saunders" = "ubs",
                    "log-extended Exponential-Geometric" = "leeg",
                    "unit-Chen" = "uchen",
                    "unit-Generalized Half-Normal-E" = "ughne",
                    "unit-Generalized Half-Normal-X" = "ughnx",
                    "unit-Gompertz" = "ugompertz",
                    "Johnson-SB" = "johnsonsb",
                    "unit-Burr-XII" = "uburrxii")

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
