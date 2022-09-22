
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to `unitquantreg` R package <img src="man/figures/unitquantreg_hex.png" align="right" alt="" width="200">

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/unitquantreg)](https://CRAN.R-project.org/package=unitquantreg)
[![R-CMD-check](https://github.com/AndrMenezes/unitquantreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AndrMenezes/unitquantreg/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/AndrMenezes/unitquantreg/branch/main/graph/badge.svg)](https://app.codecov.io/gh/AndrMenezes/unitquantreg?branch=main)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/unitquantreg?)](https://github.com/metacran/cranlogs.app)
<!-- badges: end -->

The goal of `unitquantreg` is to provide tools for estimation and
inference on parametric quantile regression models for bounded data.

We developed routines with similar interface as `stats::glm` function,
which contains estimation, inference, residual analysis, prediction, and
model comparison.

For more computation efficient the \[`dpqr`\]’s, likelihood, score and
hessian functions are vectorized and written in `C++`.

You can install the stable version from
[CRAN](https://cran.r-project.org/web/packages/unitquantreg/index.html)
with:

``` r
install.packages("unitquantreg")
```

Or you can install the development version from
[GitHub](https://github.com/) with:

``` r
if(!require(remotes)) install.packages('remotes')
remotes::install_github("AndrMenezes/unitquantreg", build_vignettes = TRUE)
```

You can then load the package

``` r
library(unitquantreg)
```

and look at user manuals typing:

``` r
vignette("unitquantreg")
vignette("structure_functionality")
```

## Citation

``` r
citation("unitquantreg")
#> 
#> To cite unitquantreg in publications use:
#> 
#>   Menezes A, Mazucheli J (2021). _unitquantreg: Parametric quantile
#>   regression models for bounded data_. R package version 0.0.3,
#>   <https://andrmenezes.github.io/unitquantreg/>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {unitquantreg: {P}arametric quantile regression models for bounded data},
#>     author = {Andr{'}e F. B. Menezes and Josmar Mazucheli},
#>     note = {R package version 0.0.3},
#>     url = {https://andrmenezes.github.io/unitquantreg/},
#>     year = {2021},
#>   }
```

## License

The `unitquantreg` package is released under the Apache License, Version
2.0. Please, see file
[`LICENSE.md`](https://github.com/AndrMenezes/unitquantreg/blob/master/LICENSE.md).
