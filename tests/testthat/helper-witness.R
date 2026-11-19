# Exact-output witnesses (inst/fixtures/statistics-*.rds) were recorded on
# aarch64 macOS, where R has no extended-precision long double: sum() and the
# native kernels that reproduce it accumulate in 64-bit doubles. On platforms
# with 80-bit long double accumulation (x86_64 Linux/Windows) the same ordered
# sums round differently by a few ulp (observed up to 1.2e-15 relative). The
# signed arc-length tent estimator sums cancelling kernel contributions and
# locates zero crossings from them, which amplifies those ulp differences
# (observed up to ~5e-14 relative within a bin; ~6.5e-16 of the row total).
# Bitwise identity is therefore required only where accumulator precision
# matches the recording; elsewhere values must agree to 1e-12 (mean relative
# difference), with identical structure and attributes. Any genuine estimator
# or binning change moves values by many orders of magnitude more.
witness_accumulation_matches <- function() !isTRUE(capabilities("long.double"))

expect_witness <- function(actual, expected, ...) {
  if (witness_accumulation_matches()) {
    expect_identical(actual, expected, ...)
  } else {
    expect_equal(actual, expected, tolerance = 1e-12, ...)
  }
}
