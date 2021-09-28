test_that("multiplication works", {
  data(water)
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
                      "unit-Burr-XII" = "uburrxii"
  )
  lt_fits <- lapply(lt_families, function(fam) {
    cat(fam, "\n")
    unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                 tau = 0.5, data = water, family = fam)
  })

  ans <- pairwise.vuong.test(lt = lt_fits)
  ans
})
