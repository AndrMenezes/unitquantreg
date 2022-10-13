# `unitquantreg` 0.0.5

- Reduced unit testing timings.
- Reduced example timings.
- Used `tryCatch()` to handle when `chol()` fails in the hessian matrix.
- Fixed `.plot_conddist()` when shape parameter is not constant.
- Included a `stop` in `.plot_conddist()` to check if all covariates have values in `at_obs` argument.
- Included the `sim_bounded` data set, a new toy simulated data set.


# `unitquantreg` 0.0.4

- Fixed unit tests for ATLAS linear algebra operation library.
- Suppressed output of `unitquantreg` when computing the numerical hessian.

# `unitquantreg` 0.0.3

- Included unit-Gumbel distribution.
- Fixed `print.vuong` method.
- Corrected Vuong test name .
