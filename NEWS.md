# AHR 1.4

## Changes

* Fixed variance estimation for restricted mean method `rmeanDiff`.

* Removed `start` parameter from function `wkm` (and all functions
  passing arguments to `wkm` such as `ahrWKM` and
  `rmeanDiff`). Replaced with `rr.subset` parameter (logical vector
  indicating which rows from the data.frame should be used for
  response rate estimation. This is more flexible than before and does
  not require `wkm` to know about recruitment times.

* Test cases are now run automatically when `R CMD check` is called.

