# Surgical refresh of the exact-output witnesses after MLV-R005/R006 (junction
# arm definition). Only rows whose property is "junction_angle" or
# "junction_type" may change; the script stops unless every other value, row,
# column and attribute of every case is identical to the current output.
# Usage: Rscript tools/mlvref/refresh_junction_witnesses.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(normalizePath(lib), .libPaths()))
library(mlvtoolbox)
refresh <- function(expected, actual, label) {
  if (identical(expected, actual)) return(list(value = expected, changed = 0L))
  stopifnot(is.data.frame(expected), identical(dim(expected), dim(actual)),
            identical(names(expected), names(actual)), identical(attributes(expected)$row.names, attributes(actual)$row.names))
  for (v in setdiff(names(expected), "value")) stopifnot(identical(expected[[v]], actual[[v]]))
  rows <- which(!(expected$value %in% actual$value & expected$value == actual$value) | is.na(expected$value) != is.na(actual$value))
  if (!all(expected$property[rows] %in% c("junction_angle", "junction_type"))) stop(label, ": non-junction rows differ")
  out <- expected; out$value[rows] <- actual$value[rows]
  if (!identical(out, actual)) stop(label, ": differences outside junction values")
  list(value = out, changed = length(rows))
}
path <- "inst/fixtures/statistics-0.2.0.rds"; w <- readRDS(path); n <- 0L
for (i in seq_along(w$cases)) {
  case <- w$cases[[i]]
  a <- do.call(getExportedValue("mlvtoolbox", case$fun), case$args)
  r <- refresh(case$expected, a, paste(path, i)); w$cases[[i]]$expected <- r$value; n <- n + r$changed
}
w$errata <- c(w$errata, "junction_angle/junction_type rows refreshed for MLV-R005/R006 (arms leave the junction neighbourhood; direction of the leaving segment; see docs/evidence/mlvref/junctions/REPORT.md)")
saveRDS(w, path); message(path, ": ", n, " junction values refreshed")
path <- "inst/fixtures/statistics-before-optimization.rds"; w <- readRDS(path); n <- 0L
for (i in seq_along(w$cases)) {
  case <- w$cases[[i]]
  a <- do.call(contour_statistics, c(list(w$drawing), case$args)); attr(a, "measurement_metadata") <- NULL
  r <- refresh(case$expected, a, paste(path, i)); w$cases[[i]]$expected <- r$value; n <- n + r$changed
}
w$errata <- c(w$errata, "junction_angle/junction_type rows refreshed for MLV-R005/R006 (arms leave the junction neighbourhood; direction of the leaving segment; see docs/evidence/mlvref/junctions/REPORT.md)")
saveRDS(w, path); message(path, ": ", n, " junction values refreshed")
