# Synthetic junction suite with analytical expectations (track A, P2).
# Each case states the correct answer derived by hand from elementary
# geometry, independent of both implementations. The script runs mlvtoolbox
# and the pinned MLV source under Octave and reports agreement of each with
# the analytical answer. Usage:
#   Rscript tools/mlvref/synthetic_junctions.R [lib]
# Writes docs/evidence/mlvref/junctions/synthetic-suite.csv.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else ".agent-work/track-a/lib"
suppressPackageStartupMessages(library(mlvtoolbox, lib.loc = lib))
work <- ".agent-work/track-a/synth"; dir.create(work, recursive = TRUE, showWarnings = FALSE)
receipt <- "docs/evidence/mlvref/junctions/synthetic-suite.csv"

seg <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)
# Ray from (cx,cy) at angle theta (degrees counterclockwise, image y down).
ray <- function(cx, cy, theta, len = 8) seg(cx, cy, cx + len * cospi(theta / 180), cy - len * sinpi(theta / 180))
J <- function(type, x, y, gaps) list(list(type = type, x = x, y = y, gaps = sort(gaps)))
none <- list()
s30 <- sinpi(1 / 6); c30 <- cospi(1 / 6)
poly <- function(v) { v <- rbind(v, v[1, ]); cbind(v[-nrow(v), , drop = FALSE], v[-1, , drop = FALSE]) }
# Closed polygon arms at (10,10): up (90) and along (10,11)->(20,20), i.e.
# 360 - atan2(9,10); bar arms 0 and 180.
pb <- atan2(9, 10) * 180 / pi; closed_gaps <- c(pb, 90, 90, 180 - pb)

cases <- list(
  L_corner = list(list(seg(0, 10, 10, 10), seg(10, 10, 10, 0)), none,
    "two arms are a corner, not a junction"),
  T = list(list(seg(0, 10, 20, 10), seg(10, 10, 10, 20)), J("T", 10, 10, c(90, 90, 180)),
    "stem ends on the bar interior"),
  T_diagonal = list(list(seg(0, 0, 10, 10), seg(5, 5, 10, 0)), J("T", 5, 5, c(90, 90, 180)),
    "perpendicular stem on a diagonal bar (MLV-M007)"),
  X = list(list(seg(0, 10, 20, 10), seg(10, 0, 10, 20)), J("X", 10, 10, c(90, 90, 90, 90)),
    "two crossing lines"),
  X_oblique = list(list(seg(0, 10, 20, 10), seg(10 - 10 * c30, 10 + 10 * s30, 10 + 10 * c30, 10 - 10 * s30)),
    J("X", 10, 10, c(30, 30, 150, 150)), "lines crossing at 30 degrees"),
  Y = list(list(ray(10, 10, 90), ray(10, 10, 210), ray(10, 10, 330)), J("Y", 10, 10, c(120, 120, 120)),
    "three contours share a start point"),
  Y_reversed_ends = list(lapply(list(ray(10, 10, 90), ray(10, 10, 210), ray(10, 10, 330)), function(z) z[, c(3, 4, 1, 2), drop = FALSE]),
    J("Y", 10, 10, c(120, 120, 120)), "three contours share an end point (reversal invariance)"),
  Arrow = list(list(ray(10, 10, 0), ray(10, 10, 30), ray(10, 10, 60)), J("Arrow", 10, 10, c(30, 30, 300)),
    "three rays within 60 degrees"),
  Star5 = list(lapply(seq(0, 288, 72), function(a) ray(10, 10, a)), J("Star", 10, 10, rep(72, 5)),
    "five rays"),
  X_from_ends = list(lapply(c(0, 90, 180, 270), function(a) ray(10, 10, a)), J("X", 10, 10, rep(90, 4)),
    "four contours share a start point"),
  psi = list(list(seg(10, 0, 10, 20), ray(10, 10, 45), ray(10, 10, 135)), J("X", 10, 10, c(45, 45, 135, 135)),
    "trident: a through line plus two branches (four arms)"),
  gap_0.9 = list(list(seg(0, 10, 20, 10), seg(10, 10.9, 10, 20)), J("T", 10, 10, c(90, 90, 180)),
    "stem 0.9 px short of the bar"),
  gap_1.1 = list(list(seg(0, 10, 20, 10), seg(10, 11.1, 10, 20)), none, "stem 1.1 px short of the bar"),
  gap_diag_1.2 = list(list(seg(0, 10, 20, 10), seg(10 + 1.2 / sqrt(2), 10 + 1.2 / sqrt(2), 18, 18)), none,
    "45-degree stem 1.2 px (Euclidean) short of the bar: beyond the default 1 px Euclidean gap; MLV's Chebyshev gap (one raster step) joins it (MLV-M018)"),
  gap_diag_0.9 = list(list(seg(0, 10, 20, 10), seg(10 + 0.9 / sqrt(2), 10 + 0.9 / sqrt(2), 18, 18)),
    J("T", 10, 10, c(45, 135, 180)), "45-degree stem 0.9 px short of the bar"),
  collinear_continuation = list(list(seg(0, 10, 10, 10), seg(10, 10, 20, 10)), none,
    "collinear end-to-end contours"),
  collinear_overlap = list(list(seg(0, 10, 12, 10), seg(8, 10, 20, 10)), none,
    "collinear overlap has no unique intersection"),
  T_at_vertex = list(list(seg(0, 10, 10, 10, 10, 10, 20, 10), seg(10, 10, 10, 20)), J("T", 10, 10, c(90, 90, 180)),
    "stem meets the bar at a bar vertex"),
  bent_bar = list(list(seg(0, 10, 10, 10, 10, 10, 20, 0), seg(10, 10, 10, 20)), J("Y", 10, 10, c(90, 135, 135)),
    "stem at the vertex of a bent bar"),
  merge_chain = list(list(seg(0, 10, 30, 10), seg(10, 0, 10, 20), seg(13, 0, 13, 20), seg(11.5, 0, 11.5, 20)),
    J("Star", 11.5, 10, c(0, 0, 0, 0, 90, 90, 90, 90)),
    "intersections at x=10,13,11.5 chain within 2 px (MLV-M006; permutation invariance)"),
  merge_chain_permuted = list(list(seg(0, 10, 30, 10), seg(10, 0, 10, 20), seg(11.5, 0, 11.5, 20), seg(13, 0, 13, 20)),
    J("Star", 11.5, 10, c(0, 0, 0, 0, 90, 90, 90, 90)), "same geometry, contours renumbered"),
  negative = list(list(seg(-6, -3, 0, -3), seg(-3, -6, -3, 0)), J("X", -3, -3, rep(90, 4)),
    "crossing at negative coordinates (MLV-M008)"),
  snap_long = list(list(seg(0, 5, 40000, 5), seg(3, 0, 3, 10)), J("X", 3, 5, rep(90, 4)),
    "crossing at parameter 7.5e-5 of a long segment, 3 px from its start (MLV-M019)"),
  zero_length_neighbour = list(list(seg(0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 10), seg(10, 10, 10, 20)),
    J("T", 10, 10, c(90, 90, 180)), "zero-length bar segment at the stem (MLV-M015/M020)"),
  offset_merge_rays = list(list(seg(20, 30, 50, 30), seg(20.6, 30.5, 20.6, 60), seg(20, 30, 0, 10)),
    J("Y", 20.4, 30.2, c(90, 135, 135)),
    "three straight rays whose pairwise intersections differ by < 1 px (MLV-R005)"),
  offset_merge_lines = list(list(seg(0, 30, 40, 30), seg(20, 10, 20, 50), seg(1, 50, 41, 10)),
    J("Star", 61 / 3, 91 / 3, c(45, 45, 45, 45, 90, 90)),
    "three straight lines crossing pairwise within 1.5 px (MLV-R005, MLV-M021)"),
  short_inside = list(list(seg(0, 10, 10, 10, 10, 10, 18, 4), seg(10, 8.5, 10, 11.5)), none,
    "3 px contour bisected by a bar vertex lies inside the 2 px neighbourhood: no arm leaves, so two arms remain (MLV-R006)"),
  short_inside_reversed = list(list(seg(0, 10, 10, 10, 10, 10, 18, 4), seg(10, 11.5, 10, 8.5)), none,
    "same, short contour reversed"),
  long_crossing_bent = list(list(seg(0, 10, 10, 10, 10, 10, 18, 4), seg(10, 7, 10, 13)),
    J("X", 10, 10, c(53.13010235415598, 90, 90, 126.86989764584402)), "6 px contour crossing a bar vertex"),
  chain_split = list(list(seg(0, 10, 10, 10), seg(10, 10, 11, 11), seg(11, 11, 20, 11)), none,
    "one polyline split into three contours with a 1.4 px middle piece (contour-partition invariance, MLV-M023)"),
  overshoot = list(list(seg(0, 10, 20, 10), seg(10, 9, 10, 10.5, 10, 10.5, 10, 20)), J("T", 10, 10, c(90, 90, 180)),
    "stem overshooting the bar by 1 px, inside the neighbourhood (MLV-M024: duplicate arm)"),
  stub_stem = list(list(seg(0, 10, 20, 10), seg(10, 10.2, 10, 11.2, 10, 11.2, 10, 20)), J("T", 10, 10, c(90, 90, 180)),
    "stem 0.2 px short of the bar whose 1 px first segment lies inside the 2 px neighbourhood"),
  stub_stem_reversed = list(list(seg(0, 10, 20, 10), seg(10, 20, 10, 11.2, 10, 11.2, 10, 10.2)), J("T", 10, 10, c(90, 90, 180)),
    "same stem drawn toward the bar: MLV takes the start point as the attachment (MLV-M022)"),
  X_vertices_near_crossing = list(list(seg(0, 10, 9, 10, 9, 10, 11, 10, 11, 10, 20, 10), seg(10, 0, 10, 20)),
    J("X", 10, 10, rep(90, 4)), "straight line with vertices 1 px either side of the crossing: one pass, two arms (MLV-M024)"),
  merged_chain_stem = list(list(seg(0, 10, 30, 10), seg(10, 0, 10, 20), seg(11.9, 0, 11.9, 20), seg(13.8, 0, 13.8, 20),
    seg(15.7, 10, 15.7, 20)), J("Star", 12.85, 10, c(0, 0, 0, 0, 0, 90, 90, 90, 90)),
    "merged chain whose last member is a T stem starting 2.85 px from the centroid: one arm, not two (MLV-R008, MLV-M025)"),
  closed_start_1 = list(list(seg(0, 10, 20, 10), poly(rbind(c(10, 11), c(10, 0), c(30, 0), c(20, 20)))), J("X", 10, 10, closed_gaps),
    "closed polygon passing 1 px below the bar crossing; polygon starts at the vertex inside the neighbourhood (MLV-R007)"),
  closed_start_2 = list(list(seg(0, 10, 20, 10), poly(rbind(c(10, 0), c(30, 0), c(20, 20), c(10, 11)))), J("X", 10, 10, closed_gaps),
    "same polygon, second start vertex"),
  closed_start_3 = list(list(seg(0, 10, 20, 10), poly(rbind(c(30, 0), c(20, 20), c(10, 11), c(10, 0)))), J("X", 10, 10, closed_gaps),
    "same polygon, third start vertex"),
  closed_start_4 = list(list(seg(0, 10, 20, 10), poly(rbind(c(20, 20), c(10, 11), c(10, 0), c(30, 0)))), J("X", 10, 10, closed_gaps),
    "same polygon, fourth start vertex"),
  merge_order = list(list(seg(9, 11.2, 11, 11.2), seg(10, 11, 10, 12.5), seg(11.5, 9, 11.5, 10.5), seg(8.5, 9, 8.5, 10.5),
    seg(5, 10, 9.5, 10, 9.5, 10, 10.5, 10, 10.5, 10, 20, 10)), none,
    "three crossings merged by a chain whose partners MLV visits in ascending index order; every small contour lies inside the 2 px neighbourhood, so only the bar's two arms leave (mlv2023 visit order)"),
  T_reversed_stem = list(list(seg(0, 10, 20, 10), seg(10, 20, 10, 10)), J("T", 10, 10, c(90, 90, 180)),
    "stem drawn toward the bar")
)

# An mlv_drawing cannot hold zero-length segments, so for cases that contain
# them Octave is also run on the cleaned geometry; that is the reference for
# rule = "mlv2023" (the uncleaned run documents MLV-M020).
zero_len <- function(z) z[, 1] == z[, 3] & z[, 2] == z[, 4]
cleaned <- lapply(cases, function(cs) lapply(cs[[1]], function(z) z[!zero_len(z), , drop = FALSE]))
dirty <- which(vapply(cases, function(cs) any(unlist(lapply(cs[[1]], zero_len))), NA))
octave_inputs <- c(lapply(cases, `[[`, 1), cleaned[dirty])
flat <- do.call(rbind, lapply(seq_along(octave_inputs), function(i)
  do.call(rbind, lapply(seq_along(octave_inputs[[i]]), function(c) cbind(i, c, octave_inputs[[i]][[c]])))))
colnames(flat) <- c("case", "contour", "x1", "y1", "x2", "y2")
cases_csv <- file.path(work, "cases.csv"); octave_csv <- file.path(work, "octave.csv")
utils::write.csv(format(as.data.frame(flat), digits = 17), cases_csv, row.names = FALSE, quote = FALSE)
status <- system2("octave-cli", c("--quiet", "--no-history", "tools/mlvref/synthetic_junctions.m", cases_csv, octave_csv),
                  stdout = TRUE, stderr = TRUE)
oct <- utils::read.csv(octave_csv, stringsAsFactors = FALSE)
mlv_ref_id <- seq_along(cases); mlv_ref_id[dirty] <- length(cases) + seq_along(dirty)
# Committed fixture for test-junction-oracle.R: cleaned geometry and the pinned
# MLV output on it (synthetic data only).
fx <- "inst/fixtures/octave"
fc <- as.data.frame(flat); fc <- fc[fc$case %in% mlv_ref_id, ]
fc$case <- names(cases)[match(fc$case, mlv_ref_id)]
utils::write.csv(format(fc, digits = 17), file.path(fx, "junction_cases.csv"), row.names = FALSE, quote = FALSE)
fo <- oct[oct$case %in% mlv_ref_id, ]; fo$case <- names(cases)[match(fo$case, mlv_ref_id)]
utils::write.csv(fo[, c("case", "junction", "x", "y", "type", "contours", "segments", "angles", "min_angle", "max_angle")],
                 file.path(fx, "junction_octave.csv"), row.names = FALSE)

gaps_of <- function(s) sort(as.numeric(strsplit(s, ";", fixed = TRUE)[[1]]))
seq_of <- function(s) as.numeric(strsplit(s, ";", fixed = TRUE)[[1]])
agrees <- function(found, expected, tol = 1e-9) {
  if (length(found) != length(expected)) return(FALSE)
  all(mapply(function(f, e) f$type == e$type && abs(f$x - e$x) <= tol && abs(f$y - e$y) <= tol &&
    length(f$gaps) == length(e$gaps) && all(abs(f$gaps - e$gaps) <= tol), found, expected))
}
# Exact MLV reproduction: same junction set; per junction type, position,
# ordered angle sequence and ordered contour/segment lists.
mlv_exact <- function(j, o, tol = 1e-9) {
  if (nrow(j) != nrow(o)) return(FALSE)
  if (!nrow(o)) return(TRUE)
  all(vapply(seq_len(nrow(o)), function(k) {
    i <- which(abs(j$x - o$x[k]) <= tol & abs(j$y - o$y[k]) <= tol)
    if (length(i) != 1L) return(FALSE)
    a <- j$angles[[i]]; b <- seq_of(o$angles[k]); m <- j$members[[i]]
    j$type[i] == o$type[k] && length(a) == length(b) && all(abs(a - b) <= tol) &&
      identical(as.numeric(m$contour_index), seq_of(o$contours[k])) &&
      identical(as.numeric(m$segment_index), seq_of(o$segments[k]))
  }, NA))
}
show <- function(found) if (!length(found)) "none" else paste(vapply(found, function(f)
  sprintf("%s@(%.6g,%.6g)[%s]", f$type, f$x, f$y, paste(signif(f$gaps, 6), collapse = " ")), ""), collapse = " + ")

rows <- lapply(seq_along(cases), function(i) {
  cs <- cases[[i]]
  d <- mlv_drawing(cs[[1]], 100000, 100000, zero_length = "drop")
  j <- detect_junctions(d)
  r <- lapply(seq_len(nrow(j)), function(k) list(type = j$type[k], x = j$x[k], y = j$y[k], gaps = sort(j$angles[[k]])))
  o <- oct[oct$case == i & oct$type != "none", ]
  oe <- if (any(o$type == "error")) o$error[1] else ""
  o <- o[o$type != "error", ]
  m <- lapply(seq_len(nrow(o)), function(k) list(type = o$type[k], x = o$x[k], y = o$y[k], gaps = gaps_of(o$angles[k])))
  oref <- oct[oct$case == mlv_ref_id[i] & !oct$type %in% c("none", "error"), ]
  data.frame(case = names(cases)[i], description = cs[[3]], expected = show(cs[[2]]),
    r = show(r), r_correct = agrees(r, cs[[2]]),
    mlv = if (nzchar(oe)) paste("error:", oe) else show(m), mlv_correct = !nzchar(oe) && agrees(m, cs[[2]]),
    r_mlv2023_reproduces_mlv = mlv_exact(detect_junctions(d, rule = "mlv2023"), oref),
    mlv_reference = if (i %in% dirty) "octave on cleaned geometry" else "octave")
})
res <- do.call(rbind, rows)
utils::write.csv(res, receipt, row.names = FALSE)
print(res[, c("case", "r_correct", "mlv_correct", "r_mlv2023_reproduces_mlv")], row.names = FALSE)
cat(sprintf("R correct %d/%d; MLV correct %d/%d; rule='mlv2023' reproduces MLV %d/%d\n", sum(res$r_correct), nrow(res),
            sum(res$mlv_correct), nrow(res), sum(res$r_mlv2023_reproduces_mlv), nrow(res)))
