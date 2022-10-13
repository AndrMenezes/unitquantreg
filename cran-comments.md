## `unitquantreg` 0.0.5

- Used `tryCatch()` to handle when `chol()` fails in the hessian matrix.
- Unit testing and examples have been re-written in order reduce the timings.


## `unitquantreg` 0.0.4

> Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_unitquantreg.html>.
Please correct before 2022-10-03 to safely retain your package on CRAN.

Minor corrections were made and tested in ATLAS linear algebra operation library.

> Please change http --> https, add trailing slashes, or follow moved
content as appropriate.

Fixed.

## `unitquantreg` 0.0.3

> Please add more details about the package functionality and implemented
methods in your Description text. If there are references describing the
methods in your package, please add these.

The DESCRIPTION file was re-written and a reference was included following the 
specified format.

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation.

The missing Rd-tags were updated.

> Please make sure that you do not change the user's options, par or working
directory.

Checked, please see the following lines:
  - 88-89 in `R/plot.R`;
  - 142-143 in `R/utils.R`.

> Please always make sure to reset to user's options().

Corrected, please see the following lines:
  - 232-235 in `inst/doc/structure_functionality.R`;
  - 20-23 in `tests/testthat/test-plot.R`.

## Test environments

* ubuntu 20.04, devel and release
* windows, release
* macOS, release

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is 10.1Mb
    sub-directories of 1Mb or more:
      libs   9.0Mb

The NOTE refers to large size of `C++` codes.
