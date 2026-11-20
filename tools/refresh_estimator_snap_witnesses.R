# Surgical refresh of the exact-output witnesses for two 0.3.0 changes:
#  1. the specification-free curvature default became estimator = "turning";
#  2. histogram binning snaps values within 1e-9 * max(1, |edge|) of a bin edge
#     onto the edge (R/statistics.R, src/statistics.cpp).
# These witnesses are regression checks of the published (mlv2023) estimator
# and of the fused native kernel that reproduces it, so their specification-free
# cases are pinned to estimator = "mlv2023" rather than re-recorded under the new
# default. Serialized specifications that no longer validate (definition bumps
# made elsewhere) are rebuilt only to verify their values; those cases are left
# untouched for tools/mlvref/mat_refresh_witnesses.R.
# The script stops unless, for every case:
#  * the pinned call differs from the stored output only in histogram values
#    whose change is reproduced independently here: each moved mass belongs to
#    a segment whose axial orientation lies within the tolerance of an
#    orientation bin edge, below it, and moves from the lower to the upper bin
#    of that edge (only orientation changes are attributable; any other changed
#    property stops the script);
#  * all other values, rows, columns and attributes are identical (up to the
#    rebuilt specification's own definitions).
# Usage: Rscript tools/refresh_estimator_snap_witnesses.R <lib> [--write]
args <- commandArgs(TRUE)
lib <- args[1]; write <- "--write" %in% args
.libPaths(c(normalizePath(lib), .libPaths()))
library(mlvtoolbox)
tol <- function(edge) 1e-9 * pmax(1, abs(edge))

rebuild_spec <- function(spec) {
  ok <- tryCatch({validate_feature_spec(spec); TRUE}, error = function(e) FALSE)
  if (ok) return(list(spec = spec, changed = character()))
  rebuilt <- do.call(feature_spec, spec$settings)
  if (!identical(spec$settings, rebuilt$settings)) stop("spec rebuild changed settings")
  changed <- names(spec$definitions)[spec$definitions != rebuilt$definitions]
  list(spec = rebuilt, changed = paste0(changed, ": ", spec$definitions[changed], " -> ", rebuilt$definitions[changed]))
}
pin <- function(fun_args) {
  if (is.null(fun_args$spec)) fun_args$estimator <- "mlv2023"
  fun_args
}
# Replace every embedded copy of the old spec semantics by the rebuilt ones so
# that attributes can be compared exactly.
swap_spec <- function(x, old, new) {
  if (identical(x, old)) return(new)
  if (is.list(x) && !is.data.frame(x)) {
    at <- attributes(x); x <- lapply(x, swap_spec, old = old, new = new); attributes(x) <- at
  }
  if (!is.null(attributes(x))) for (a in setdiff(names(attributes(x)), c("names", "row.names", "class", "dim", "dimnames")))
    attr(x, a) <- swap_spec(attr(x, a), old, new)
  x
}
strip <- function(s) {s$settings$pixels_per_degree <- NULL; s}

# Independent attribution of orientation-bin moves caused by snapping.
orientation_moves <- function(drawing, bins, by, normalize) {
  st <- mlvtoolbox:::.segment_table(drawing)
  axial <- st$orientation %% 180
  scaled <- axial * bins / 180 + .5
  old <- (floor(scaled) %% bins) + 1
  nearest <- floor(scaled + .5); edge <- (nearest - .5) * 180 / bins
  near <- abs(axial - edge) <= tol(edge)
  new <- old; new[near] <- (nearest[near] %% bins) + 1
  g <- if (by == "drawing") rep(1L, nrow(st)) else st$contour_index
  ng <- if (by == "drawing") 1L else length(drawing$ids)
  w <- st$length
  if (normalize) {
    tot <- tapply(st$length, factor(g, levels = seq_len(ng)), sum); w <- w / tot[g] * 10000
  }
  delta <- matrix(0, ng, bins)
  for (k in which(new != old)) {
    delta[g[k], old[k]] <- delta[g[k], old[k]] - w[k]; delta[g[k], new[k]] <- delta[g[k], new[k]] + w[k]
  }
  list(delta = delta, moved = which(new != old), orientation = st$orientation[new != old])
}
check_values <- function(expected, actual, drawing, a, label) {
  # Returns a description of the change, after verifying its attribution.
  if (is.data.frame(expected) && "property" %in% names(expected)) {
    for (v in setdiff(names(expected), "value")) if (!identical(expected[[v]], actual[[v]])) stop(label, ": column ", v, " differs")
    rows <- which(expected$value != actual$value | is.na(expected$value) != is.na(actual$value))
    if (!length(rows)) return(NULL)
    if (!all(expected$property[rows] == "orientation")) stop(label, ": non-orientation rows differ: ", paste(unique(expected$property[rows]), collapse = ","))
    o <- expected$property == "orientation"
    mold <- matrix(expected$value[o], ncol = a$bins, byrow = TRUE); mnew <- matrix(actual$value[o], ncol = a$bins, byrow = TRUE)
    m <- orientation_moves(drawing, a$bins, a$by, a$normalize)
  } else {
    cols <- grep("^orientation_", names(expected), value = TRUE)
    for (v in setdiff(names(expected), cols)) if (!identical(expected[[v]], actual[[v]])) stop(label, ": column ", v, " differs")
    if (identical(expected[cols], actual[cols])) return(NULL)
    rows <- NULL
    mold <- as.matrix(expected[cols]); mnew <- as.matrix(actual[cols])
    norm <- if (is.null(a$normalize)) TRUE else a$normalize
    m <- orientation_moves(drawing, length(cols), "contour", norm)
  }
  if (!length(m$moved)) stop(label, ": orientation changed without a segment near a bin edge")
  if (!isTRUE(all.equal(unname(mnew), unname(mold + m$delta), tolerance = 1e-12)))
    stop(label, ": orientation change is not explained by snapping")
  changed <- which(mold != mnew)
  sprintf("orientation %s: %s -> %s (segments at %s deg)",
          if (is.null(rows)) "columns" else "rows",
          paste(signif(mold[changed], 10), collapse = ","), paste(signif(mnew[changed], 10), collapse = ","),
          paste(unique(sprintf("%.17g", m$orientation)), collapse = ", "))
}
refresh_case <- function(fun, fun_args, expected, label, drawing = fun_args$x) {
  rebuilt <- NULL; fun_args_original <- fun_args
  if (!is.null(fun_args$spec)) {
    r <- rebuild_spec(fun_args$spec)
    if (length(r$changed)) rebuilt <- list(old = fun_args$spec, new = r$spec, changed = r$changed)
    fun_args$spec <- r$spec
  }
  pinned <- pin(fun_args)
  actual <- do.call(fun, pinned)
  if (identical(expected, actual)) return(list(args = pinned, value = expected, note = NULL, spec = NULL))
  if (!is.null(rebuilt) && identical(unname(unclass(expected)), unname(unclass(actual))))
    return(list(args = fun_args_original, value = expected, note = NULL, spec = rebuilt$changed))
  cmp <- expected
  if (!is.null(rebuilt)) cmp <- swap_spec(swap_spec(cmp, rebuilt$old, rebuilt$new), strip(rebuilt$old), strip(rebuilt$new))
  if (!identical(attributes(cmp)[setdiff(names(attributes(cmp)), "row.names")],
                 attributes(actual)[setdiff(names(attributes(actual)), "row.names")]))
    stop(label, ": attributes differ beyond the rebuilt specification")
  note <- check_values(cmp, actual, drawing, pinned, label)
  if (!is.null(rebuilt)) {
    # Definition drift from other changes (e.g. a MAT definition bump) is not
    # folded into this refresh: such a case is verified under the rebuilt spec
    # and left untouched, and only an unchanged value set is accepted.
    if (!is.null(note)) stop(label, ": values changed under a rebuilt specification; refresh after the definition bump is resolved")
    return(list(args = fun_args_original, value = expected, note = NULL, spec = rebuilt$changed))
  }
  list(args = pinned, value = actual, note = note, spec = NULL)
}

report <- character(); spec_notes <- character()
path <- "inst/fixtures/statistics-0.2.0.rds"; w <- readRDS(path)
for (i in seq_along(w$cases)) {
  case <- w$cases[[i]]
  r <- refresh_case(getExportedValue("mlvtoolbox", case$fun), case$args, case$expected, paste(path, "case", i))
  w$cases[[i]]$args <- r$args; w$cases[[i]]$expected <- r$value
  if (!is.null(r$note)) report <- c(report, paste0("case ", i, " (", case$fun, ", bins = ", case$args$bins, "): ", r$note))
  spec_notes <- union(spec_notes, r$spec)
}
cat(path, "\n", paste(" ", report, collapse = "\n"), "\n", sep = "")
if (length(spec_notes)) cat("  left untouched (embedded specs differ only by other definition changes, values verified unchanged): ",
  paste(spec_notes, collapse = "; "), "\n", sep = "")
errata <- "specification-free cases pinned to estimator = 'mlv2023' (0.3.0 default became 'turning'; these cases witness the published estimator and its fused kernel); every pinned value unchanged"
if (length(report)) errata <- c(errata, paste0("orientation values refreshed for bin-edge snapping (a value within 1e-9 * max(1, |edge|) of a bin edge is binned as if on it): ",
  "segments at 89.999999999998991 and 269.99999999999898 deg, vertical up to rounding residue (x offsets ~9e-14 px), now enter the upper bin of the 90-degree edge with the exact 90s, as their counterparts at 90.000000000001037 deg already did (odd bin counts). ",
  paste(report, collapse = "; ")))
w$errata <- c(w$errata, setdiff(errata, w$errata))
if (write) saveRDS(w, path)

path <- "inst/fixtures/statistics-before-optimization.rds"; b <- readRDS(path); n <- 0L
for (i in seq_along(b$cases)) {
  case <- b$cases[[i]]
  stats <- function(...) {z <- contour_statistics(b$drawing, ...); attr(z, "measurement_metadata") <- NULL; z}
  r <- refresh_case(stats, case$args, case$expected, paste(path, "case", i), drawing = b$drawing)
  b$cases[[i]]$args <- r$args[names(r$args) != "x"]; b$cases[[i]]$expected <- r$value
  if (!is.null(r$note)) {n <- n + 1L; cat(path, "case", i, r$note, "\n")}
}
features <- feature_table(b$drawing, c("length", "orientation", "curvature", "parallelism"), 13, estimator = "mlv2023")
attr(features, "feature_schema") <- NULL; attr(features, "measurement_metadata") <- NULL
if (!identical(features, b$features)) stop(path, ": pinned features differ")
b$errata <- (function(e) c(b$errata, setdiff(e, b$errata)))(paste0("specification-free cases pinned to estimator = 'mlv2023' (0.3.0 default became 'turning'); ",
  if (n) paste(n, "cases refreshed for bin-edge snapping") else "every value unchanged, including under bin-edge snapping"))
if (write) saveRDS(b, path)
message(if (write) "witnesses written" else "dry run: no files written (use --write)")
