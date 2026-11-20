# Surgical refresh of exact-output witnesses after the MAT default change
# (aof = "centered", morphology = "hmt8"; track B, docs/evidence/mlvref/mat).
# Witness cases that carry a serialized feature_spec hold the 0.2.x definitions
# (morphology "octave_image_pratt_2.20.0", no aof); those specifications are now
# rejected as noncanonical by design. This script rebuilds each stored
# specification from its own settings (gaining only the new mat$aof/morphology
# defaults and their definition strings), recomputes the case, and stops unless
# the result equals the stored expectation after substituting the rebuilt
# specification into the embedded schema/metadata. No numerical value may change.
# Usage: Rscript tools/mlvref/mat_refresh_witnesses.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(normalizePath(lib), .libPaths()))
library(mlvtoolbox)
swap_spec <- function(x, spec) {
  fs <- attr(x, "feature_schema"); if (!is.null(fs$measurement_spec)) { fs$measurement_spec <- spec; attr(x, "feature_schema") <- fs }
  mm <- attr(x, "measurement_metadata")
  if (!is.null(mm$schema$measurement_spec)) { mm$schema$measurement_spec <- spec; attr(x, "measurement_metadata") <- mm }
  x
}
path <- "inst/fixtures/statistics-0.2.0.rds"; w <- readRDS(path); n <- 0L
for (i in seq_along(w$cases)) {
  case <- w$cases[[i]]; k <- which(vapply(case$args, inherits, TRUE, "mlv_feature_spec"))
  if (!length(k)) next
  stopifnot(length(k) == 1L)
  old <- case$args[[k]]
  spec <- do.call(feature_spec, old$settings)
  keep <- setdiff(names(old$definitions), c("aof", "morphology", "branch_rule"))
  stopifnot(identical(spec$definitions[keep], old$definitions[keep]))
  case$args[[k]] <- spec
  actual <- do.call(getExportedValue("mlvtoolbox", case$fun), case$args)
  # Embedded schemas carry the calibration-free semantics of the specification.
  expected <- swap_spec(case$expected, mlvtoolbox:::.spec_semantics(spec))
  if (!identical(actual, expected)) stop(path, " case ", i, ": differences beyond the specification definitions")
  w$cases[[i]] <- case; w$cases[[i]]$expected <- expected; n <- n + 1L
}
note <- paste("serialized feature specifications rebuilt for the MAT default changes (mat aof='centered',",
  "morphology='minimal', branch_rule='invariant'; definitions centered-tie-mean-1 /",
  "hmt8-then-sequential-simple-pixel-removal-1 / direction-symmetric-tiemax-graph-2); every value unchanged;",
  "see docs/evidence/mlvref/mat (track B)")
w$errata <- c(w$errata[!grepl("serialized feature specifications rebuilt for the MAT default", w$errata)], note)
saveRDS(w, path); message(path, ": ", n, " cases with embedded specifications refreshed; values unchanged")
