library(testthat)
library(unitquantreg)

lt_families <- list(
  "unit-Weibull" = "uweibull",
  "Kumaraswamy" = "kum",
  "unit-Logistic" = "ulogistic",
  "unit-Birnbaum-Saunders" = "ubs",
  "log-extended Exponential-Geometric" = "leeg",
  "unit-Chen" = "uchen",
  "unit-Generalized Half-Normal-E" = "ughne",
  "unit-Generalized Half-Normal-X" = "ughnx",
  "unit-Gompertz" = "ugompertz",
  "Johnson-SB" = "johnsonsb",
  "unit-Burr-XII" = "uburrxii",
  "arc-secant hyperbolic Weibull" = "ashw"
  )

test_check("unitquantreg")
