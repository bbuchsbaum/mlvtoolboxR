# Track A (P2): compare mlvtoolbox segment geometry and junctions with stored
# MATLAB MLV outputs exported by tools/mlvref/export.m, and with the pinned
# MLV source executed under Octave by tools/mlvref/drift_pinned.m.
#
# Usage:
#   Rscript tools/mlvref/compare_junctions.R [lib] [reference] [cores]
#     lib        library holding the mlvtoolbox build under test
#     reference  "stored" (default; .agent-work/mlvref) or "pinned"
#                (.agent-work/track-a/pinned, drift sample only) or
#                "pinned-vs-stored" (P1 drift: no R involved)
#
# Natural data stay under .agent-work/. Only the aggregate receipts written to
# docs/evidence/mlvref/junctions/ are intended for Git; per-junction detail is
# written to .agent-work/track-a/detail/.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else ".agent-work/track-a/lib"
reference <- if (length(args) >= 2) args[2] else "stored"
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L
suppressPackageStartupMessages(library(mlvtoolbox, lib.loc = lib))
stopifnot(reference %in% c("stored", "pinned", "pinned-vs-stored"))

export_root <- ".agent-work/mlvref"
pinned_root <- ".agent-work/track-a/pinned"
detail_root <- file.path(".agent-work/track-a/detail", reference)
receipt_root <- "docs/evidence/mlvref/junctions"
dir.create(detail_root, recursive = TRUE, showWarnings = FALSE)
dir.create(receipt_root, recursive = TRUE, showWarnings = FALSE)

# Tolerances. Inputs are integer pixel coordinates (stored natively as double),
# so every intersection is a ratio of integers and both implementations compute
# it in IEEE double from identical operands; 1e-6 px is about 1e8 ulp at
# coordinate 1000 and far below any geometric effect. Angles: 1e-6 degree.
tol_pos <- 1e-6
tol_ang <- 1e-6
tol_geom <- 1e-9
diag_radius <- 2 # loose radius used only to diagnose near-matches

# R's own decimal parser (R_strtod, also used by read.csv, data.table and
# vroom here) is not correctly rounded for 17-digit input on this platform:
# "870.99999999999977" parses one ulp away from the nearest double. The
# export writes %.17g, so numeric columns are re-parsed with the C library's
# correctly rounded strtod to recover the MATLAB doubles exactly.
exact_num <- Rcpp::cppFunction('NumericVector exact_num(CharacterVector s) {
  NumericVector out(s.size());
  for (R_xlen_t i = 0; i < s.size(); ++i) {
    if (CharacterVector::is_na(s[i]) || !std::strlen(s[i])) { out[i] = NA_REAL; continue; }
    out[i] = std::strtod(s[i], nullptr);
  }
  return out;
}', includes = "#include <cstdlib>\n#include <cstring>")
read_gz <- function(path) {
  gz <- paste0(path, ".gz")
  z <- if (file.exists(gz)) utils::read.csv(gzfile(gz), colClasses = "character")
       else utils::read.csv(path, colClasses = "character")
  for (v in intersect(names(z), c("drawing", "contour", "segment", "junction", "width", "height", "num_contours")))
    z[[v]] <- as.integer(z[[v]])
  for (v in intersect(names(z), c("x1", "y1", "x2", "y2", "x", "y", "orientation", "length", "curvature",
                                  "min_angle", "max_angle", "parallelism", "separation", "mirror")))
    z[[v]] <- exact_num(z[[v]])
  z
}
complete_sets <- function(root) basename(dirname(Sys.glob(file.path(root, "*", "COMPLETE"))))
split_ids <- function(s) lapply(strsplit(as.character(s), ";", fixed = TRUE), exact_num)

# Zero-length under the constructor's 64-ulp relative tolerance (MLV-R009).
zero_seg <- function(p) {
  tol <- function(a, b) 64 * .Machine$double.eps * pmax(abs(a), abs(b))
  abs(p$x2 - p$x1) <= tol(p$x1, p$x2) & abs(p$y2 - p$y1) <= tol(p$y1, p$y2)
}
# Build an R drawing from stored contours. Zero-length segments are dropped
# explicitly (MLV-M015); `map` keeps original (contour, segment) indices so all
# comparisons refer to the stored MATLAB numbering.
build_drawing <- function(s, width, height) {
  s <- s[order(s$contour, s$segment), ]
  pieces <- split(s, s$contour)
  keep <- lapply(pieces, function(p) p[!zero_seg(p), ])
  keep <- keep[vapply(keep, nrow, 1L) > 0L]
  mats <- lapply(keep, function(p) as.matrix(p[, c("x1", "y1", "x2", "y2")]))
  d <- mlv_drawing(unname(mats), width, height, ids = names(keep))
  map <- lapply(keep, function(p) p$segment)
  list(drawing = d, map = map, contour_map = as.integer(names(keep)))
}

# R junctions translated into stored numbering, with the MLV-style
# participant list: intersecting segments plus their adjacent segment when the
# junction lies within endpoint radius of a shared vertex (cf.
# computeJunctionAnglesTypes.m).
r_junctions <- function(b, ...) {
  j <- detect_junctions(b$drawing, ...)
  if (!nrow(j)) return(data.frame(x = double(), y = double(), type = character(),
    contours = character(), segments = character(), angles = character()))
  out <- data.frame(x = j$x, y = j$y, type = j$type,
    contours = vapply(j$members, function(m) paste(b$contour_map[m$contour_index], collapse = ";"), ""),
    segments = vapply(j$members, function(m) paste(mapply(function(ci, si) b$map[[ci]][si],
      m$contour_index, m$segment_index), collapse = ";"), ""),
    angles = vapply(j$angles, function(a) paste(sprintf("%.17g", a), collapse = ";"), ""),
    stringsAsFactors = FALSE)
  attr(out, "raw") <- j
  out
}

# Greedy one-to-one matching by distance (ascending), radius diag_radius.
match_points <- function(P, Q) {
  if (!nrow(P) || !nrow(Q)) return(list(i = integer(), j = integer(), d = double()))
  D <- sqrt(outer(P[, 1], Q[, 1], "-")^2 + outer(P[, 2], Q[, 2], "-")^2)
  cand <- which(D <= diag_radius, arr.ind = TRUE)
  if (!nrow(cand)) return(list(i = integer(), j = integer(), d = double()))
  cand <- cand[order(D[cand], cand[, 1], cand[, 2]), , drop = FALSE]
  ui <- logical(nrow(P)); uj <- logical(nrow(Q)); keep <- logical(nrow(cand))
  for (k in seq_len(nrow(cand))) {
    a <- cand[k, 1]; b <- cand[k, 2]
    if (!ui[a] && !uj[b]) { ui[a] <- uj[b] <- keep[k] <- TRUE }
  }
  cand <- cand[keep, , drop = FALSE]
  list(i = cand[, 1], j = cand[, 2], d = D[cand])
}

pair_key <- function(cs, ss) {
  vapply(seq_along(cs), function(k) paste(sort(unique(paste(cs[[k]], ss[[k]], sep = ":"))), collapse = " "), "")
}
set_key <- function(v) vapply(v, function(x) paste(sort(unique(x)), collapse = " "), "")
angle_equal <- function(a, b) mapply(function(x, y) length(x) == length(y) &&
  all(abs(sort(x) - sort(y)) <= tol_ang), a, b)

# Compare junction tables A (reference) and B (candidate) for one drawing.
compare_one <- function(A, B, drawing) {
  m <- match_points(as.matrix(A[, c("x", "y")]), as.matrix(B[, c("x", "y")]))
  rows <- list()
  if (length(m$i)) {
    a <- A[m$i, ]; b <- B[m$j, ]
    ac <- split_ids(a$contours); bc <- split_ids(b$contours)
    as <- split_ids(a$segments); bs <- split_ids(b$segments)
    aa <- split_ids(a$angles); ba <- split_ids(b$angles)
    rows[[1]] <- data.frame(drawing = drawing, ref = m$i, cand = m$j, dist = m$d,
      pos_ok = m$d <= tol_pos, type_a = a$type, type_b = b$type,
      contours_ok = set_key(ac) == set_key(bc),
      pairs_ok = pair_key(ac, as) == pair_key(bc, bs),
      pairs_sub = mapply(function(x, y) all(strsplit(y, " ")[[1]] %in% strsplit(x, " ")[[1]]),
                         pair_key(ac, as), pair_key(bc, bs)),
      n_arms_a = lengths(aa), n_arms_b = lengths(ba),
      angles_ok = angle_equal(aa, ba),
      minmax_ok = abs(a$min_angle - vapply(ba, min, 0)) <= tol_ang &
                  abs(a$max_angle - vapply(ba, max, 0)) <= tol_ang,
      stringsAsFactors = FALSE)
  }
  ua <- setdiff(seq_len(nrow(A)), m$i); ub <- setdiff(seq_len(nrow(B)), m$j)
  list(matched = if (length(rows)) rows[[1]] else NULL,
       ref_only = if (length(ua)) data.frame(drawing = drawing, ref = ua, type = A$type[ua]) else NULL,
       cand_only = if (length(ub)) data.frame(drawing = drawing, cand = ub, type = B$type[ub]) else NULL)
}

with_minmax <- function(J) {
  if (!nrow(J)) { J$min_angle <- double(); J$max_angle <- double(); return(J) }
  a <- split_ids(J$angles)
  if (is.null(J$min_angle)) J$min_angle <- vapply(a, min, 0)
  if (is.null(J$max_angle)) J$max_angle <- vapply(a, max, 0)
  J
}

geometry_diff <- function(s, b) {
  st <- mlvtoolbox:::.segment_table(b$drawing)
  st$contour <- b$contour_map[st$contour_index]
  st$segment <- mapply(function(ci, si) b$map[[ci]][si], st$contour_index, st$segment_index)
  key_s <- paste(s$contour, s$segment); key_r <- paste(st$contour, st$segment)
  i <- match(key_r, key_s)
  ref <- s[i, ]
  od <- abs(((st$orientation - ref$orientation) + 180) %% 360 - 180)
  # zero-length neighbours in the stored numbering (MLV-M015): MLV measures
  # curvature against them; R drops them before measuring.
  zl <- paste(s$contour, s$segment)[zero_seg(s)]
  zero_nb <- paste(st$contour, st$segment + 1) %in% zl | paste(st$contour, st$segment - 1) %in% zl
  cbad <- !is.finite(ref$curvature) | abs(st$curvature - ref$curvature) > tol_geom * pmax(1, abs(ref$curvature))
  data.frame(curv_bad_zero_neighbour = sum(cbad & zero_nb), curv_bad_subpixel_length = sum(cbad & !zero_nb & st$length < 1e-9),
    curv_bad_other = sum(cbad & !zero_nb & st$length >= 1e-9),
    ori_bad_subpixel_length = sum(od > tol_geom & st$length < 1e-9), segments = nrow(s), compared = nrow(st), dropped_zero = nrow(s) - nrow(st),
    ori_bad = sum(od > tol_geom), len_bad = sum(abs(st$length - ref$length) > tol_geom),
    curv_bad = sum(!is.finite(ref$curvature) | abs(st$curvature - ref$curvature) >
                     tol_geom * pmax(1, abs(ref$curvature))),
    ori_max = max(c(0, od)), len_max = max(c(0, abs(st$length - ref$length))),
    curv_max = max(c(0, abs(st$curvature - ref$curvature)[is.finite(ref$curvature)])))
}

run_dataset <- function(ds) {
  segs <- read_gz(file.path(export_root, ds, "segments.csv"))
  draws <- read_gz(file.path(export_root, ds, "drawings.csv"))
  stored_j <- read_gz(file.path(export_root, ds, "junctions.csv"))
  if (reference != "stored") {
    pj <- read_gz(file.path(pinned_root, ds, "junctions.csv"))
    sample <- utils::read.csv(file.path(pinned_root, ds, "sample.csv"))$drawing
  } else sample <- draws$drawing
  out <- parallel::mclapply(sample, function(k) {
    s <- segs[segs$drawing == k, ]
    w <- draws$width[draws$drawing == k]; h <- draws$height[draws$drawing == k]
    A <- with_minmax(stored_j[stored_j$drawing == k, ])
    if (reference == "pinned-vs-stored") {
      B <- with_minmax(pj[pj$drawing == k, ])
      return(list(cmp = compare_one(A, B, k), geom = NULL, nA = nrow(A), nB = nrow(B)))
    }
    if (reference == "pinned") A <- with_minmax(pj[pj$drawing == k, ])
    b <- build_drawing(s, w, h)
    B <- with_minmax(r_junctions(b))
    list(cmp = compare_one(A, B, k), geom = cbind(drawing = k, geometry_diff(s, b)),
         nA = nrow(A), nB = nrow(B))
  }, mc.cores = cores)
  bad <- vapply(out, inherits, NA, "try-error")
  if (any(bad)) stop(ds, ": ", paste(unique(vapply(out[bad], as.character, "")), collapse = "\n"))
  pick <- function(f) do.call(rbind, lapply(out, function(o) o$cmp[[f]]))
  res <- list(matched = pick("matched"), ref_only = pick("ref_only"), cand_only = pick("cand_only"),
              geom = do.call(rbind, lapply(out, `[[`, "geom")),
              nA = sum(vapply(out, `[[`, 0, "nA")), nB = sum(vapply(out, `[[`, 0, "nB")))
  saveRDS(res, file.path(detail_root, paste0(ds, ".rds")))
  res
}

summarise <- function(ds, r) {
  m <- r$matched
  exact <- if (is.null(m)) 0L else sum(m$pos_ok & m$type_a == m$type_b & m$contours_ok & m$angles_ok & m$minmax_ok)
  data.frame(dataset = ds, ref_junctions = r$nA, cand_junctions = r$nB,
    matched_2px = NROW(m), exact = exact,
    pos_ok = sum(m$pos_ok), type_ok = sum(m$type_a == m$type_b),
    contours_ok = sum(m$contours_ok), pairs_ok = sum(m$pairs_ok), angles_ok = sum(m$angles_ok),
    ref_only = NROW(r$ref_only), cand_only = NROW(r$cand_only),
    seg_total = if (is.null(r$geom)) NA else sum(r$geom$segments),
    seg_dropped_zero = if (is.null(r$geom)) NA else sum(r$geom$dropped_zero),
    ori_bad = if (is.null(r$geom)) NA else sum(r$geom$ori_bad),
    len_bad = if (is.null(r$geom)) NA else sum(r$geom$len_bad),
    curv_bad = if (is.null(r$geom)) NA else sum(r$geom$curv_bad),
    curv_bad_zero_neighbour = if (is.null(r$geom)) NA else sum(r$geom$curv_bad_zero_neighbour),
    curv_bad_subpixel_length = if (is.null(r$geom)) NA else sum(r$geom$curv_bad_subpixel_length),
    curv_bad_other = if (is.null(r$geom)) NA else sum(r$geom$curv_bad_other),
    ori_bad_subpixel_length = if (is.null(r$geom)) NA else sum(r$geom$ori_bad_subpixel_length))
}

sets <- complete_sets(export_root)
if (reference != "stored") sets <- intersect(sets, complete_sets(pinned_root))
summary <- do.call(rbind, lapply(sets, function(ds) {
  t0 <- Sys.time(); r <- run_dataset(ds); s <- summarise(ds, r)
  message(sprintf("%s: %d ref, %d exact (%.0fs)", ds, s$ref_junctions, s$exact,
                  as.numeric(Sys.time() - t0, units = "secs")))
  s
}))
tag <- if (length(args) >= 4) args[4] else ""
out_file <- file.path(receipt_root, sprintf("summary-%s%s.csv", reference, tag))
utils::write.csv(summary, out_file, row.names = FALSE)
tot <- colSums(summary[, -1], na.rm = TRUE)
print(tot)
