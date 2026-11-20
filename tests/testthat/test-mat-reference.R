# Track-B reference qualification (docs/evidence/mlvref/mat/REPORT.md). All inputs are
# synthetic; expectations are analytical or invariance oracles, not stored MATLAB values.

test_that("width-one rasterization is independent of segment direction (MLV-R010)", {
 one <- function(z) which(rasterize(mlv_drawing(list(matrix(z, 1)), 10, 10)))
 # Exact line y = 1 + (x-1)/4 crosses y = 1.5 at x = 3: positive slope takes the larger row.
 expect_identical(one(c(1, 1, 5, 2)), c(1L, 11L, 22L, 32L, 42L))
 expect_identical(one(c(5, 2, 1, 1)), c(1L, 11L, 22L, 32L, 42L))
 # Negative slope takes the smaller row at the same half-pixel tie.
 expect_identical(one(c(1, 2, 5, 1)), c(2L, 12L, 21L, 31L, 41L))
 expect_identical(one(c(5, 1, 1, 2)), c(2L, 12L, 21L, 31L, 41L))
 # Steep segments: ties on the column, same rule.
 expect_identical(one(c(1, 1, 2, 5)), one(c(2, 5, 1, 1)))
 set.seed(11)
 z <- matrix(sample(1:30, 400, TRUE), ncol = 4); z <- z[(z[, 1] - z[, 3])^2 + (z[, 2] - z[, 4])^2 > 0, ]
 rows <- function(M) lapply(seq_len(nrow(M)), function(i) M[i, , drop = FALSE])
 a <- rasterize(mlv_drawing(rows(z), 30, 30))
 b <- rasterize(mlv_drawing(rows(z[, c(3, 4, 1, 2)]), 30, 30))
 expect_identical(a, b)
 # Every segment is 8-connected with exactly max(|dx|,|dy|)+1 pixels when unclipped.
 for (i in seq_len(nrow(z))) {
  p <- which(rasterize(mlv_drawing(list(z[i, , drop = FALSE]), 30, 30)))
  expect_equal(length(p), max(abs(z[i, 3] - z[i, 1]), abs(z[i, 4] - z[i, 2])) + 1)
 }
})

test_that("hmt8 thinning removes Pratt ribbons and preserves topology", {
 m <- matrix(FALSE, 14, 24); m[3:12, 3:22] <- TRUE      # even height 10
 two_by_two <- function(s) sum(s[-1, -1] & s[-nrow(s), -1] & s[-1, -ncol(s)] & s[-nrow(s), -ncol(s)])
 expect_gt(two_by_two(skeletonize(m, method = "pratt")), 0)   # documented Pratt ribbon
 h <- skeletonize(m, method = "hmt8")
 expect_equal(two_by_two(h), 0)
 expect_true(all(h <= m)); expect_true(any(h))
 comps <- function(s) { lab <- 0L; seen <- matrix(FALSE, nrow(s), ncol(s))
  for (i in which(s)) if (!seen[i]) { lab <- lab + 1L; q <- i; seen[i] <- TRUE
   while (length(q)) { p <- q[1]; q <- q[-1]; r <- (p - 1) %% nrow(s) + 1; c <- (p - 1) %/% nrow(s) + 1
    for (dr in -1:1) for (dc in -1:1) { rr <- r + dr; cc <- c + dc
     if (rr >= 1 && rr <= nrow(s) && cc >= 1 && cc <= ncol(s)) { j <- rr + nrow(s) * (cc - 1)
      if (s[j] && !seen[j]) { seen[j] <- TRUE; q <- c(q, j) } } } } }; lab }
 expect_equal(comps(h), 1L)
 ring <- matrix(FALSE, 20, 20); ring[4:17, 4:17] <- TRUE; ring[8:13, 8:13] <- FALSE
 hr <- skeletonize(ring, method = "hmt8")
 expect_equal(comps(hr), 1L); expect_false(hr[10, 10]); expect_equal(two_by_two(hr), 0)
 # odd-width bar: exact midline
 bar <- matrix(FALSE, 11, 30); bar[4:8, 3:28] <- TRUE
 expect_true(all(skeletonize(bar, method = "hmt8")[6, 8:23]))
 expect_identical(skeletonize(h, method = "hmt8"), h)
 expect_error(skeletonize(m, method = "nope"))
})

test_that("MLV AOF sampling is misregistered; the default centred tie-mean AOF finds the medial row (MLV-M027)", {
 m <- matrix(FALSE, 41, 41); m[10, 4:38] <- m[22, 4:38] <- TRUE   # exact medial row 16
 thr <- 2 / pi * sin(28 * pi / 180) * 60
 a <- compute_mat(m, trace_branches = FALSE, aof = "mlv2023")
 c <- compute_mat(m, trace_branches = FALSE)
 expect_identical(c$parameters$aof, "centered")
 expect_identical(which(a$AOF[, 20] > thr), c(16L, 17L))     # MLV: two-row ridge, shifted half a pixel
 expect_identical(which(c$AOF[, 20] > thr), 16L)             # analytical medial row
 expect_equal(c$AOF[15, 20], c$AOF[17, 20], tolerance = 1e-12)
 expect_true(all(c$skeleton[16, 12:30])); expect_false(any(c$skeleton[c(15, 17), 12:30]))
 # vertical version: medial column 16
 cv <- compute_mat(t(m), trace_branches = FALSE)
 expect_identical(which(cv$AOF[20, ] > thr), 16L)
})

test_that("the default AOF and its threshold mask are equivariant under all eight lattice symmetries", {
 n <- 31
 D4 <- list(identity = function(x) x, rot90 = function(x) t(x[nrow(x):1, ]), rot180 = function(x) x[nrow(x):1, ncol(x):1],
            rot270 = function(x) t(x)[ncol(x):1, ], flip_rows = function(x) x[nrow(x):1, ], flip_cols = function(x) x[, ncol(x):1],
            transpose = function(x) t(x), antitranspose = function(x) t(x[nrow(x):1, ncol(x):1]))
 inverse <- list(identity = "identity", rot90 = "rot270", rot180 = "rot180", rot270 = "rot90", flip_rows = "flip_rows",
                 flip_cols = "flip_cols", transpose = "transpose", antitranspose = "antitranspose")
 # Inputs are fixed points of the (order-dependent) hmt8 ink thinning under every
 # transform, so the comparison isolates the flux definition. Diagonal outlines,
 # crossings and dots create many exact equidistant ties.
 dia <- matrix(FALSE, n, n); k <- 0:10
 dia[cbind(16 - k, 6 + k)] <- dia[cbind(6 + k, 16 + k)] <- dia[cbind(16 + k, 26 - k)] <- dia[cbind(26 - k, 16 - k)] <- TRUE
 ex <- matrix(FALSE, n, n); ex[cbind(5:27, 5:27)] <- TRUE; ex[cbind(5:27, 27:5)] <- TRUE
 plus <- matrix(FALSE, n, n); plus[16, 5:27] <- TRUE; plus[5:27, 16] <- TRUE
 par <- matrix(FALSE, n, n); par[8, 4:28] <- par[20, 4:28] <- TRUE
 dots <- matrix(FALSE, n, n); dots[c(8, 24), c(8, 24)] <- TRUE; dots[16, 16] <- TRUE
 masks <- list(dia, ex, plus, par, dots)
 for (mask in masks) for (g in D4) expect_identical(skeletonize(g(mask)), g(mask))
 ii <- 3:(n - 2)
 for (mask in masks) {
  A <- compute_mat(mask, trace_branches = FALSE, min_component = 1)$AOF
  M <- compute_mat(mask, trace_branches = FALSE, aof = "mlv2023", min_component = 1)$AOF
  worst_mlv <- 0
  for (g in names(D4)) {
   back <- D4[[inverse[[g]]]]
   B <- back(compute_mat(D4[[g]](mask), trace_branches = FALSE, min_component = 1)$AOF)
   expect_lt(max(abs(A[ii, ii] - B[ii, ii])), 1e-11)
   expect_identical(A[ii, ii] > 18, B[ii, ii] > 18)
   if (g != "identity") worst_mlv <- max(worst_mlv,
     max(abs(M[ii, ii] - back(compute_mat(D4[[g]](mask), trace_branches = FALSE, aof = "mlv2023", min_component = 1)$AOF)[ii, ii])))
  }
  expect_gt(worst_mlv, 1)                                       # mlv2023 is not equivariant
 }
 expect_error(compute_mat(par, aof = "other"))
 expect_error(compute_mat(par, samples = 2), "samples")
})

test_that("MAT methods are recorded, specification-controlled and versioned", {
 m <- matrix(FALSE, 30, 30); m[8, 5:25] <- m[20, 5:25] <- TRUE
 x <- compute_mat(m, min_branch_points = 2)
 expect_identical(x$parameters[c("aof", "morphology", "aof_definition", "morphology_definition")],
   list(aof = "centered", morphology = "minimal", aof_definition = "centered-tie-mean-1",
        morphology_definition = "hmt8-then-sequential-simple-pixel-removal-1"))
 p <- compute_mat(m, aof = "mlv2023", morphology = "pratt")$parameters
 expect_identical(c(p$aof_definition, p$morphology_definition), c("mlv2023-computeAOF-lowest-index", "octave_image_pratt_2.20.0"))
 d <- mlv_drawing(list(matrix(c(5, 8, 25, 8), 1), matrix(c(5, 20, 25, 20), 1)), 30, 30)
 s <- feature_spec("parallelism")
 expect_identical(s$settings$mat[c("aof", "morphology")], list(aof = "centered", morphology = "minimal"))
 expect_identical(unname(s$definitions[c("aof", "morphology")]), c("centered-tie-mean-1", "hmt8-then-sequential-simple-pixel-removal-1"))
 expect_identical(compute_mat(m, morphology = "hmt8")$parameters$morphology_definition, "sequential_hmt8_matlab_skel_order-1")
 expect_identical(unname(feature_spec("parallelism", mat = list(branch_rule = "mlv2023"))$definitions[["tracing"]]),
                  "kovesi_edgelink_cleanedgelist")
 expect_identical(unname(s$definitions[["tracing"]]), "pixel_graph_no_corner_cut")
 expect_identical(s$settings$mat$branch_rule, "invariant")
 expect_identical(unname(s$definitions[["branch_rule"]]), "direction-symmetric-tiemax-graph-1")
 expect_identical(x$parameters[c("branch_rule", "branch_rule_definition", "tracing")],
   list(branch_rule = "invariant", branch_rule_definition = "direction-symmetric-tiemax-graph-1", tracing = "pixel_graph_no_corner_cut"))
 compat <- feature_spec("parallelism", mat = list(aof = "mlv2023", branch_rule = "mlv2023"))
 mc <- compute_mat(d, spec = compat)
 expect_identical(mc$parameters[c("branch_rule", "branch_rule_definition", "tracing")],
   list(branch_rule = "mlv2023", branch_rule_definition = "mlv2023-edgelink-lastwriter-1", tracing = "kovesi_edgelink_cleanedgelist"))
 expect_error(compute_mat_properties(mc, spec = feature_spec("parallelism", mat = list(aof = "mlv2023"))), "mismatch")
 expect_error(compute_mat(d, spec = s, branch_rule = "mlv2023"), "not both")
 expect_error(feature_spec("parallelism", mat = list(branch_rule = "edgelink")), "mat.branch_rule")
 stale <- s; stale$definitions <- stale$definitions[setdiff(names(stale$definitions), "branch_rule")]
 expect_error(validate_feature_spec(stale), "noncanonical")
 expect_error(compute_mat(d, spec = s, morphology = "pratt"), "not both")
 expect_error(compute_mat(d, spec = s, aof = "mlv2023"), "not both")
 ms <- compute_mat(d, spec = s)
 expect_identical(ms$parameters$aof, "centered")
 expect_s3_class(compute_mat_properties(ms), "mlv_mat")
 legacy <- feature_spec("parallelism", mat = list(aof = "mlv2023", morphology = "pratt"))
 ml <- compute_mat(d, spec = legacy)
 expect_identical(c(ml$parameters$aof, ml$parameters$morphology), c("mlv2023", "pratt"))
 expect_false(identical(ml$AOF, ms$AOF))
 # A MAT built with one method set cannot be scored under a specification declaring another.
 expect_error(compute_mat_properties(ml, spec = s), "mismatch")
 # A cached 0.2.x specification (Pratt definition under the new defaults) is rejected.
 old <- s; old$definitions[["morphology"]] <- "octave_image_pratt_2.20.0"
 expect_error(validate_feature_spec(old), "noncanonical")
 old2 <- s; old2$settings$mat$aof <- NULL; old2$settings$mat$morphology <- NULL; old2$definitions <- old2$definitions[setdiff(names(old2$definitions), "aof")]
 expect_error(validate_feature_spec(old2))
 expect_error(feature_spec("parallelism", mat = list(aof = "corner")), "mat.aof")
 expect_error(feature_spec("parallelism", mat = list(morphology = "zhang")), "mat.morphology")
 # extraction thinning follows the same morphology
 img <- matrix(1, 20, 30); img[6:13, 4:26] <- 0
 expect_identical(trace_contours(img)$source$morphology, "hmt8-then-sequential-simple-pixel-removal-1")
 expect_identical(trace_contours(img, spec = legacy)$source$morphology, "octave_image_pratt_2.20.0")
})
test_that("coincident consecutive branch points map without NaN subscripts (MLV-R011)", {
 d <- mlv_drawing(list(matrix(c(5, 10, 35, 10), 1), matrix(c(5, 20, 35, 20), 1)), 40, 30)
 m <- compute_mat(d, min_branch_points = 2)
 m <- compute_mat_properties(m)
 b <- m$branches[[1]]; k <- ceiling(nrow(b) / 2)
 dup <- b[c(seq_len(k), k, (k + 1):nrow(b)), ]
 m2 <- m; m2$branches[[1]] <- dup
 for (p in c("parallelism", "separation", "mirror")) m2$branches[[1]][[p]] <- score_branch(dup, p)
 expect_no_error(y <- map_mat_properties(d, m2))
 expect_true(any(is.finite(y$mat_images$parallelism)))
})

test_that("mlv2023 branch scores depend on traversal direction; arclength scores do not (MLV-M028)", {
 x <- 1:25; b <- data.frame(x = x, y = round(3 * sin(x / 4)), radius = 2 + .1 * x + .3 * sin(x / 2))
 r <- b[nrow(b):1, ]
 for (p in c("parallelism", "separation", "mirror")) {
  expect_gt(max(abs(score_branch(b, p) - rev(score_branch(r, p)))), 1e-3)
  expect_lt(max(abs(score_branch(b, p, estimator = "arclength") - rev(score_branch(r, p, estimator = "arclength"))),
                na.rm = TRUE), 1e-12)
 }
})

test_that("skeleton pixels on ink (signed distance <= 0) are never branch points", {
 d <- mlv_drawing(list(matrix(c(5, 10, 35, 10), 1), matrix(c(5, 20, 35, 20), 1)), 40, 30)
 m <- compute_mat(d, min_branch_points = 2)
 expect_identical(m$parameters$skeleton_ink_excluded, 0L)
 ink <- which(m$distance_map <= 0)[1:3]
 m$skeleton[ink] <- TRUE
 for (rule in c("invariant", "mlv2023")) {
  b <- trace_skeleton(m, 2L, rule = rule)
  expect_identical(attr(b, "skeleton_ink_excluded"), 3L)
  idx <- unlist(lapply(b, function(z) z$row + nrow(m$skeleton) * (z$column - 1)))
  expect_false(any(idx %in% ink))
  expect_true(all(unlist(lapply(b, `[[`, "radius")) > 0))
 }
})

test_that("the edgelink/cleanedgelist port reproduces pinned MLV tracing exactly", {
 dir <- system.file("fixtures", "edgelink", package = "mlvtoolbox")
 for (nm in c("grid", "fan", "kite")) {
  m <- unname(as.matrix(utils::read.csv(file.path(dir, paste0(nm, ".mask.csv")), header = FALSE)) > 0)
  o <- utils::read.csv(file.path(dir, paste0(nm, ".octave.branches.csv")), header = FALSE)
  expected <- unname(lapply(split(o[, 2:3], o[, 1]), function(z) { z <- unname(as.matrix(z)); storage.mode(z) <- "integer"; z }))
  expect_identical(lapply(mlvtoolbox:::.trace_edgelink(m, 10), unname), expected)
  raw <- mlvtoolbox:::cpp_edgelink_trace(m)
  expect_identical(mlvtoolbox:::cpp_cleanedgelist(raw, 10), mlvtoolbox:::.mlv_cleanedgelist(raw, 10))
  # traced through the public API on a MAT carrying this skeleton (all medial)
  mat <- structure(list(skeleton = m, distance_map = matrix(2, nrow(m), ncol(m)), AOF = matrix(0, nrow(m), ncol(m)),
                        parameters = list(branch_rule = "mlv2023")), class = "mlv_mat")
  b <- trace_skeleton(mat, 10L)
  expect_identical(lapply(b, function(z) unname(cbind(z$row, z$column))), expected)
 }
 # duplicated merge points are part of the published chain (MLV-M030)
 m <- unname(as.matrix(utils::read.csv(file.path(dir, "fan.mask.csv"), header = FALSE)) > 0)
 expect_gt(sum(vapply(mlvtoolbox:::.trace_edgelink(m, 10), function(z) sum(diff(z[, 1]) == 0 & diff(z[, 2]) == 0), 1)), 0)
})

test_that("published-compatible mirror reproduces MLV fitLineSegments including MLV-M016", {
 rd <- function(name) unname(as.matrix(utils::read.csv(system.file("fixtures", "parity", name, package = "mlvtoolbox"), header = FALSE)))
 for (name in c("constant", "varying", "curved", "terminal_bend", "short")) {
  z <- rd(paste0(name, ".branch.csv")); ref <- rd(paste0(name, ".octave.scores.csv"))
  b <- data.frame(row = z[, 1], column = z[, 2], radius = z[, 3])   # MATLAB X = row, Y = column
  expect_lt(max(abs(mlvtoolbox:::.mirror_mlv(b, 10) - ref[, 3])), 1e-12)
 }
})

test_that("default branch rule: per-pixel scores are invariant to branch reversal and renumbering", {
 d <- mlv_drawing(list(matrix(c(8, 8, 52, 8, 52, 8, 52, 52, 52, 52, 8, 52, 8, 52, 8, 8), 4, byrow = TRUE),
                       matrix(c(30, 8, 30, 52), 1), matrix(c(5, 30, 55, 24), 1)), 60, 60)
 props <- c("parallelism", "separation", "mirror", "taper")
 m <- compute_mat(d, min_branch_points = 3)
 expect_identical(m$parameters$branch_rule, "invariant")
 base <- map_mat_properties(d, compute_mat_properties(m, props))
 same <- function(y) for (p in props) {
  expect_identical(is.na(y$mat_images[[p]]), is.na(base$mat_images[[p]]))
  expect_lt(max(abs(y$mat_images[[p]] - base$mat_images[[p]]), na.rm = TRUE), 1e-12)
 }
 rev_m <- m; rev_m$branches <- lapply(m$branches, function(b) b[nrow(b):1, ])
 same(map_mat_properties(d, compute_mat_properties(rev_m, props)))
 set.seed(3); perm <- m; perm$branches <- m$branches[sample(length(m$branches))]
 same(map_mat_properties(d, compute_mat_properties(perm, props)))
 o <- c(3L, 1L, 2L); dp <- mlv_drawing(d$contours[o], d$width, d$height)
 yp <- map_mat_properties(dp, compute_mat_properties(compute_mat(dp, min_branch_points = 3), props))
 same(yp)
 expect_identical(yp$mat_properties$parallelism, base$mat_properties$parallelism[o])
 # the published-compatible rule is direction dependent (MLV-M028)
 mm <- compute_mat_properties(compute_mat(d, min_branch_points = 3, branch_rule = "mlv2023"), props)
 mr <- mm; mr$branches <- lapply(mm$branches, function(b) b[nrow(b):1, ])
 a <- map_mat_properties(d, mm)$mat_images$parallelism
 r <- map_mat_properties(d, compute_mat_properties(mr, props))$mat_images$parallelism
 expect_gt(max(abs(a - r), na.rm = TRUE), 1e-3)
})

test_that("default branch rule is equivariant under the eight lattice symmetries", {
 n <- 41
 tf <- list(id = function(x, y) cbind(x, y), rot90 = function(x, y) cbind(n + 1 - y, x),
            rot180 = function(x, y) cbind(n + 1 - x, n + 1 - y), rot270 = function(x, y) cbind(y, n + 1 - x),
            flipx = function(x, y) cbind(n + 1 - x, y), flipy = function(x, y) cbind(x, n + 1 - y),
            transpose = function(x, y) cbind(y, x), anti = function(x, y) cbind(n + 1 - y, n + 1 - x))
 tfd <- function(d, g) mlv_drawing(lapply(d$contours, function(z) cbind(g(z[, 1], z[, 2]), g(z[, 3], z[, 4]))), n, n)
 tfi <- function(M, g) { out <- M; rc <- arrayInd(seq_along(M), dim(M)); p <- g(rc[, 2], rc[, 1]); out[cbind(p[, 2], p[, 1])] <- M; out }
 drawings <- list(
  mlv_drawing(list(matrix(c(5, 14, 37, 14), 1), matrix(c(5, 28, 37, 28), 1)), n, n),
  mlv_drawing(list(matrix(c(8, 8, 34, 8, 34, 8, 34, 34, 34, 34, 8, 34, 8, 34, 8, 8), 4, byrow = TRUE)), n, n),
  mlv_drawing(list(matrix(c(5, 21, 37, 21), 1), matrix(c(21, 5, 21, 37), 1)), n, n),
  mlv_drawing(list(matrix(c(6, 6, 36, 36), 1), matrix(c(6, 36, 36, 6), 1)), n, n),
  mlv_drawing(list(matrix(c(8, 8, 34, 8, 34, 8, 34, 34), 2, byrow = TRUE)), n, n))
 diamond <- mlv_drawing(list(matrix(c(21, 6, 36, 21, 36, 21, 21, 36, 21, 36, 6, 21, 6, 21, 21, 6), 4, byrow = TRUE)), n, n)
 check <- function(A, B) for (p in names(A)) {
  a <- tfi(A[[p]], g); expect_identical(is.na(a), is.na(B[[p]])); expect_lt(max(abs(a - B[[p]]), na.rm = TRUE), 1e-12) }
 # Branch rule alone: transform the MAT fields (hmt8 thinning itself is order dependent).
 for (d in drawings) {
  m <- compute_mat(d, min_branch_points = 3, trace_branches = FALSE)
  A <- map_mat_properties(d, compute_mat_properties(m))$mat_images
  for (g in tf) {
   mt <- m; for (f in c("distance_map", "AOF", "skeleton", "mask")) mt[[f]] <- tfi(m[[f]], g)
   check(A, map_mat_properties(tfd(d, g), compute_mat_properties(mt))$mat_images)
  }
 }
 # Whole default pipeline on a drawing whose thinned flux is itself equivariant.
 a <- analyse_mat(diamond, min_branch_points = 3)
 for (g in tf) {
  b <- analyse_mat(tfd(diamond, g), min_branch_points = 3)
  expect_identical(tfi(a$mat$skeleton, g), b$mat$skeleton)
  check(a$drawing$mat_images, b$drawing$mat_images)
 }
})

test_that("isolated closed loops are scored independently of start pixel and direction", {
 n <- 41
 ring <- matrix(FALSE, n, n); th <- seq(0, 2 * pi, length.out = 4000)
 ring[cbind(round(21 + 12 * cos(th)), round(21 + 9 * sin(th)))] <- TRUE
 ring <- skeletonize(ring)
 d <- mlv_drawing(list(matrix(c(3, 3, 39, 3), 1)), n, n)      # drawing only supplies the ink raster
 m <- compute_mat(d, min_branch_points = 3, trace_branches = FALSE)
 m$skeleton <- ring; m$distance_map[] <- 3 + (row(ring) %% 3) / 2; m$branches <- trace_skeleton(m, 3L)
 loop <- Filter(function(b) mlvtoolbox:::.is_isolated_loop(b, ring), m$branches)
 expect_length(loop, 1L); expect_length(m$branches, 1L)
 props <- c("parallelism", "separation", "mirror", "taper")
 base <- map_mat_properties(d, compute_mat_properties(m, props))$mat_images
 b <- m$branches[[1]]; u <- b[-nrow(b), ]; k <- nrow(u)
 for (start in c(2L, k %/% 3L, k - 1L)) for (reverse in c(FALSE, TRUE)) {
  ord <- c(start:k, seq_len(start - 1L)); if (reverse) ord <- rev(ord)
  mr <- m; mr$branches <- list(u[c(ord, ord[1]), ])
  y <- map_mat_properties(d, compute_mat_properties(mr, props))$mat_images
  for (p in props) {
   expect_identical(is.na(y[[p]]), is.na(base[[p]]))
   expect_lt(max(abs(y[[p]] - base[[p]]), na.rm = TRUE), 1e-12)
  }
 }
 # the published rule keeps MLV's start dependence
 mc <- m; mc$parameters$branch_rule <- "mlv2023"
 a <- compute_mat_properties(mc, "parallelism")$score_images$parallelism
 mc$branches <- list(u[c(5:k, 1:4, 5), ])
 r <- compute_mat_properties(mc, "parallelism")$score_images$parallelism
 expect_gt(max(abs(a - r), na.rm = TRUE), 1e-6)
})

test_that("the default thinning leaves no simple non-end pixel, is homotopic and idempotent (MLV-R013)", {
 two_by_two <- function(s) sum(s[-1, -1] & s[-nrow(s), -1] & s[-1, -ncol(s)] & s[-nrow(s), -ncol(s)])
 redundant <- function(s) mlvtoolbox:::cpp_redundant_pixels(s)
 comps <- function(s, conn8 = TRUE) { lab <- 0L; seen <- matrix(FALSE, nrow(s), ncol(s))
  off <- if (conn8) as.matrix(expand.grid(-1:1, -1:1)) else rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  for (i in which(s)) if (!seen[i]) { lab <- lab + 1L; q <- i; seen[i] <- TRUE
   while (length(q)) { p <- q[1]; q <- q[-1]; r <- (p - 1) %% nrow(s) + 1; c <- (p - 1) %/% nrow(s) + 1
    for (k in seq_len(nrow(off))) { rr <- r + off[k, 1]; cc <- c + off[k, 2]
     if (rr >= 1 && rr <= nrow(s) && cc >= 1 && cc <= ncol(s)) { j <- rr + nrow(s) * (cc - 1)
      if (s[j] && !seen[j]) { seen[j] <- TRUE; q <- c(q, j) } } } } }; lab }
 holes <- function(s) { p <- matrix(FALSE, nrow(s) + 2, ncol(s) + 2); p[-c(1, nrow(p)), -c(1, ncol(p))] <- s; comps(!p, FALSE) - 1L }
 sq <- matrix(FALSE, 8, 8); sq[3:6, 3:6] <- TRUE
 expect_gt(two_by_two(skeletonize(sq, method = "hmt8")), 0)          # reviewer case: hmt8 leaves a 2x2 block
 expect_identical(two_by_two(skeletonize(sq)), 0L)
 # 2x2 blocks can remain, but only where all four pixels are non-simple: the
 # redundant-pixel count (simple pixels with >= 2 neighbours) is zero there too.
 set.seed(7); kept <- 0L
 for (k in 1:400) { h <- sample(3:25, 1); w <- sample(3:25, 1); m <- matrix(runif(h * w) < runif(1, .3, .9), h, w)
  s <- skeletonize(m); expect_identical(redundant(s), 0L)
  if (two_by_two(s) > 0) { kept <- kept + 1L; expect_identical(comps(s), comps(m)); expect_identical(holes(s), holes(m)) } }
 expect_gt(kept, 0L)
 set.seed(19)
 for (k in 1:40) {
  m <- matrix(FALSE, 40, 50)
  for (j in 1:sample(2:5, 1)) { r <- sample(3:30, 1); c <- sample(3:40, 1); m[r:(r + sample(1:8, 1)), c:(c + sample(1:9, 1))] <- TRUE }
  if (k %% 3 == 0) m[sample(length(m), 150)] <- TRUE
  s <- skeletonize(m)
  expect_identical(redundant(s), 0L)
  expect_true(all(s <= m))
  expect_identical(comps(s), comps(m)); expect_identical(holes(s), holes(m))
  expect_identical(skeletonize(s), s)
 }
})

test_that("a junction-anchored loop (lollipop) keeps the maximum of its doubled pixel", {
 poly <- function(xy) { xy <- rbind(xy, xy[1, ]); cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE]) }
 circ <- function(cx, cy, r, k = 24) { t <- seq(0, 2 * pi, length = k + 1)[-1]; cbind(cx + r * cos(t), cy + r * sin(t)) }
 t <- seq(0.35, 2 * pi - 0.35, length = 40); outer <- rbind(cbind(40 + 18 * cos(t), 35 + 18 * sin(t)), c(40 + 34, 35))
 d <- mlv_drawing(list(poly(circ(40, 35, 7)), poly(outer)), 90, 70)
 m <- compute_mat(d, min_branch_points = 4)
 medial <- m$skeleton & m$distance_map > 0
 lol <- which(vapply(m$branches, function(b) { n <- nrow(b); b$row[1] == b$row[n] && b$column[1] == b$column[n] &&
   !mlvtoolbox:::.is_isolated_loop(b, medial) }, TRUE))
 expect_gt(length(lol), 0)
 P <- c("parallelism", "separation", "mirror", "taper")
 A <- compute_mat_properties(m, P); m2 <- m; m2$branches[lol] <- lapply(m2$branches[lol], function(b) b[nrow(b):1, ])
 B <- compute_mat_properties(m2, P)
 for (p in P) expect_lt(max(abs(A$score_images[[p]] - B$score_images[[p]]), na.rm = TRUE), 1e-12)
 b <- A$branches[[lol[1]]]; j <- c(b$row[1], b$column[1])
 expect_equal(A$score_images$parallelism[j[1], j[2]], max(b$parallelism[c(1, nrow(b))]))
})

test_that("loop detection uses the medial pixels that tracing used", {
 poly <- function(xy) { xy <- rbind(xy, xy[1, ]); cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE]) }
 circ <- function(cx, cy, r, k = 24) { t <- seq(0, 2 * pi, length = k + 1)[-1]; cbind(cx + r * cos(t), cy + r * sin(t)) }
 d <- mlv_drawing(list(poly(circ(30, 30, 6)), poly(circ(30, 30, 16, 40))), 60, 60)
 m <- compute_mat(d, min_branch_points = 4)
 medial <- m$skeleton & m$distance_map > 0
 L <- which(vapply(m$branches, mlvtoolbox:::.is_isolated_loop, TRUE, skeleton = medial)); expect_length(L, 1L)
 b <- m$branches[[L]]; r <- b$row[1]; c <- b$column[1]
 run <- function(mm) map_mat_properties(d, compute_mat_properties(mm, c("parallelism", "mirror")))$mat_images
 for (dr in -1:1) for (dc in -1:1) { rr <- r + dr; cc <- c + dc; if (m$skeleton[rr, cc]) next
  m2 <- m; m2$skeleton[rr, cc] <- TRUE; m2$distance_map[rr, cc] <- 0      # an ink pixel touching the start
  m2$branches <- trace_skeleton(m2, 4)
  j <- which(vapply(m2$branches, function(z) z$row[1] == r && z$column[1] == c && nrow(z) == nrow(b), TRUE))
  if (!length(j)) next
  A <- run(m2); m3 <- m2; u <- m3$branches[[j]]; u <- u[-nrow(u), ]; u <- u[c(7:nrow(u), 1:6), ]; m3$branches[[j]] <- rbind(u, u[1, ])
  B <- run(m3)
  for (p in names(A)) expect_lt(max(abs(A[[p]] - B[[p]]), na.rm = TRUE), 1e-12)
 }
})

test_that("centred AOF requires samples divisible by 4 and treats all borders alike (MLV-M027)", {
 m <- matrix(FALSE, 12, 20); m[c(1, 3), 3:18] <- TRUE
 a <- compute_mat(m, min_component = 1, trace_branches = FALSE)
 b <- compute_mat(m[12:1, ], min_component = 1, trace_branches = FALSE)
 expect_lt(max(abs(a$AOF - b$AOF[12:1, ])), 1e-12)                 # every row, including borders
 expect_lt(max(abs(a$AOF - compute_mat(m[, 20:1], min_component = 1, trace_branches = FALSE)$AOF[, 20:1])), 1e-12)
 sq <- matrix(FALSE, 25, 25); sq[c(1, 25), ] <- TRUE; sq[, c(2, 24)] <- TRUE
 A <- compute_mat(sq, min_component = 1, trace_branches = FALSE)$AOF
 expect_lt(max(abs(A - t(compute_mat(t(sq), min_component = 1, trace_branches = FALSE)$AOF))), 1e-12)
 expect_error(compute_mat(m, samples = 30), "multiple of 4")
 expect_error(feature_spec("parallelism", mat = list(samples = 50)), "multiple of 4")
 expect_s3_class(compute_mat(m, samples = 30, aof = "mlv2023", trace_branches = FALSE), "mlv_mat")
 expect_s3_class(feature_spec("parallelism", mat = list(samples = 50, aof = "mlv2023")), "mlv_feature_spec")
})

test_that("isolated-loop scoring is linear in loop length and matches the start average for mirror", {
 n <- 120; ring <- matrix(FALSE, n, n); th <- seq(0, 2 * pi, length.out = 20000)
 ring[cbind(round(60 + 50 * cos(th)), round(60 + 40 * sin(th)))] <- TRUE
 ring <- skeletonize(ring)
 d <- mlv_drawing(list(matrix(c(3, 3, 110, 3), 1)), n, n)
 m <- compute_mat(d, min_branch_points = 3, trace_branches = FALSE)
 m$skeleton <- ring; m$distance_map[] <- 3; m$branches <- trace_skeleton(m, 3L)
 expect_length(m$branches, 1L); expect_gt(nrow(m$branches[[1]]), 250)
 elapsed <- system.time(x <- compute_mat_properties(m, c("parallelism", "separation", "mirror", "taper")))[["elapsed"]]
 expect_lt(elapsed, 5)
 b <- m$branches[[1]]; u <- b[-nrow(b), ]; k <- nrow(u)
 brute <- numeric(k)
 for (s in seq_len(k)) for (rv in c(FALSE, TRUE)) {
  ord <- c(s:k, seq_len(s - 1L)); if (rv) ord <- rev(ord)
  z <- score_branch(u[c(ord, ord[1]), ], "mirror"); v <- z[seq_len(k)]; v[1] <- (z[1] + z[k + 1]) / 2
  brute[ord] <- brute[ord] + v
 }
 expect_lt(max(abs(x$branches[[1]]$mirror[seq_len(k)] - brute / (2 * k))), 1e-12)
})

test_that("long isolated loops keep start and direction invariance within 1e-12", {
 n <- 600; ring <- matrix(FALSE, n, n); th <- seq(0, 2 * pi, length.out = 60000)
 ring[cbind(round(300 + 290 * cos(th)), round(300 + 270 * sin(th)))] <- TRUE
 ring <- skeletonize(ring)
 m <- structure(list(skeleton = ring, distance_map = 3 + (row(ring) %% 3) / 2 + (col(ring) %% 5) / 7,
                     AOF = matrix(0, n, n), mask = matrix(FALSE, n, n), branches = list(),
                     parameters = list(branch_rule = "invariant", min_branch_points = 3L)), class = "mlv_mat")
 m$branches <- trace_skeleton(m, 3L)
 expect_length(m$branches, 1L); expect_gt(nrow(m$branches[[1]]), 1500)
 # Local estimators (circular windows); the mirror start average is checked
 # against brute force on a shorter loop above (it costs O(m) fits per loop).
 P <- c("parallelism", "separation", "taper")
 A <- compute_mat_properties(m, P)$score_images
 b <- m$branches[[1]]; u <- b[-nrow(b), ]; k <- nrow(u)
 for (ord in list(c(97:k, 1:96), c(800:k, 1:799), k:1)) {
  m2 <- m; m2$branches <- list(u[c(ord, ord[1]), ])
  B <- compute_mat_properties(m2, P)$score_images
  for (p in P) expect_lt(max(abs(A[[p]] - B[[p]]), na.rm = TRUE), 1e-12)
 }
})
