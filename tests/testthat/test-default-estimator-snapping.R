circle_drawing <- function(n, spacing = 1, r = 10) {
 w <- rep(spacing, length.out = n); t <- c(0, cumsum(w)) / sum(w) * 2 * pi
 xy <- cbind(20 + r * cos(t), 20 + r * sin(t)); xy[n + 1, ] <- xy[1, ]
 mlv_drawing(list(cbind(xy[-(n + 1), ], xy[-1, ])), 40, 40)
}

test_that("vertex turning is the default curvature estimator everywhere", {
 roof <- mlv_drawing(list(rbind(c(5, 20, 15, 10), c(15, 10, 25, 20))), 30, 30)
 expect_identical(feature_spec()$settings$estimator, "turning")
 integrated <- function(...) with(as.data.frame(roof, ...), sum(curvature * length))
 expect_equal(integrated(), 90)
 expect_equal(integrated(estimator = "mlv2023"), 180)
 expect_equal(sum(with(compute_contour_properties(roof)$segment_properties, curvature * length)), 90)
 # The specification-free path equals the turning specification.
 spec <- feature_spec(c("orientation", "length", "curvature", "horizontal_vertical"), bins = 5,
                      normalize = FALSE, length_range = c(2, 60))
 free <- contour_statistics(roof, spec$settings$properties, bins = 5, length_range = c(2, 60))
 specified <- contour_statistics(roof, spec = spec)
 expect_equal(free$value, specified$value, tolerance = 1e-12)
 expect_identical(feature_schema(feature_table(roof))$estimator, "vertex-dual-1")
 expect_identical(feature_schema(feature_table(roof, estimator = "mlv2023"))$estimator, "mlv2023-port-0.1.1")
 expect_error(contour_statistics(roof, spec = spec, estimator = "mlv2023"), "not both")
 expect_error(feature_table(roof, spec = spec, estimator = "turning"), "not both")
 expect_error(compute_contour_properties(roof, spec = spec, estimator = "turning"), "not both")
 expect_error(contour_statistics(roof, estimator = "arclength"), "should be one of")
 # Models trained under one estimator are not silently applied under the other.
 model <- bind_feature_model(function(f) rep(0, nrow(f)), feature_schema(feature_table(roof, estimator = "mlv2023")))
 expect_error(predict_contours(roof, model), "schema|estimator|incompatible")
 expect_equal(predict_contours(roof, model, estimator = "mlv2023"), 0)
})

test_that("turning histograms are stable under resampling; mlv2023 histograms are not", {
 bin_of <- function(d, e) which(contour_statistics(d, "curvature", estimator = e)$value > 0)
 for (n in c(16, 32, 64, 128)) {
  uniform <- circle_drawing(n); uneven <- circle_drawing(n, c(1, 2))
  expect_identical(bin_of(uniform, "turning"), 4L)
  expect_identical(bin_of(uneven, "turning"), 4L)
  expect_equal(with(as.data.frame(uneven), sum(curvature * length)), 360, tolerance = 1e-10)
  # The published estimator splits an unevenly sampled circle into two bins.
  expect_identical(bin_of(uniform, "mlv2023"), 4L)
  expect_identical(bin_of(uneven, "mlv2023"), c(3L, 5L))
 }
 # On a closed irregular polygon the published estimator skips the closing
 # vertex and counts the last one twice; turning sums the exterior angles.
 tri <- mlv_drawing(list(rbind(c(1, 1, 11, 1), c(11, 1, 1, 11), c(1, 11, 1, 1))), 12, 12)
 expect_equal(with(as.data.frame(tri), sum(curvature * length)), 360)
 expect_equal(with(as.data.frame(tri, estimator = "mlv2023"), sum(curvature * length)), 405)
})

test_that("the default turning estimator accepts exact retracing", {
 spur <- mlv_drawing(list(rbind(c(1, 1, 11, 1), c(11, 1, 6, 1))), 12, 12)
 p <- as.data.frame(spur)
 expect_equal(sum(p$curvature * p$length), 180)
 expect_true(all(is.na(p$signed_curvature[p$curvature > 0])))
 expect_equal(sum(contour_statistics(spur, "curvature")$value), 15)
 expect_error(curvature_profile(spur, "arclength"), "reversal")
 expect_error(contour_turning(spur), "reversal")
})

test_that("bin-edge snapping moves values by at most the tolerance", {
 snap <- mlvtoolbox:::.snap_to_edges; index <- mlvtoolbox:::.bin_index
 tol <- mlvtoolbox:::.bin_edge_tolerance
 edges <- seq(0, 120, length.out = 9)
 set.seed(31)
 near <- rep(edges, each = 20) + rep(edges, each = 20) * sample(-40:40, 180, TRUE) * .Machine$double.eps +
  sample(c(-1, 1), 180, TRUE) * 1e-12
 far <- c(runif(200, 0, 120), edges[-1] - 1e-6, edges + 1e-6, (head(edges, -1) + tail(edges, -1)) / 2)
 values <- c(near, far, -5, 200, NA)
 snapped <- snap(values, edges)
 moved <- which(!is.na(values) & snapped != values)
 expect_true(all(snapped[moved] %in% edges))
 expect_true(all(abs(snapped[moved] - values[moved]) <= tol(snapped[moved])))
 # Values clearly inside a bin are unaffected.
 expect_identical(snap(far, edges), far)
 expect_identical(index(far, edges), findInterval(far, edges, all.inside = TRUE))
 expect_identical(index(near, edges), findInterval(rep(edges, each = 20), edges, all.inside = TRUE))
 expect_true(is.na(index(NA_real_, edges)))
 # Infinite edges are never snapping targets.
 expect_identical(snap(c(0, 1e300), c(-Inf, Inf)), c(0, 1e300))
})

test_that("mathematically equal junction angles share a bin", {
 index <- mlvtoolbox:::.bin_index
 edges <- seq(0, 120, length.out = 9)
 for (angle in c(60, 90, 120)) {
  computed <- c(acos(cos(angle * pi / 180)) * 180 / pi, atan2(sin(angle * pi / 180), cos(angle * pi / 180)) * 180 / pi,
                angle * (1 - 4 * .Machine$double.eps), angle * (1 + 4 * .Machine$double.eps))
  expect_identical(index(computed, edges), rep(index(angle, edges), 4))
 }
 expect_identical(index(59.999999999999993, edges), index(60, edges))
 expect_identical(index(60, edges), 5L)
 # End to end: junctions drawn with trigonometric coordinates have smallest
 # gaps a few ulps below 60 or 120 degrees (e.g. 59.999999999999943); they
 # enter the same bin as an exact gap (60 enters [60, 75), the upper bin).
 star <- function(dirs) lapply(dirs, function(a) rbind(c(20, 20, 20 + 12 * cos(a * pi / 180), 20 - 12 * sin(a * pi / 180))))
 for (dirs in list(c(0, 60, 120), seq(0, 300, 60), c(0, 60, 180, 240), c(90, 210, 330), c(17, 137, 257))) {
  d <- mlv_drawing(star(dirs), 40, 40)
  j <- detect_junctions(d); exact <- if (min(diff(sort(dirs))) == 60) 60 else 120
  expect_equal(j$min_angle, exact, tolerance = 1e-12)
  expect_identical(contour_statistics(d, "junction_angle")$value, as.numeric(seq_len(8) == index(exact, edges)), info = paste(dirs, collapse = ","))
 }
 t_junction <- mlv_drawing(list(rbind(c(5, 20, 35, 20)), rbind(c(20, 20, 20, 35))), 40, 40)
 expect_identical(contour_statistics(t_junction, "junction_angle")$value, as.numeric(seq_len(8) == 7))
})

test_that("orientation snapping agrees between the R and fused native paths", {
 oi <- mlvtoolbox:::.orientation_index
 expect_identical(oi(c(44.999999999999993, 45, 45.000000000000007), 2), c(2, 2, 2))
 expect_identical(oi(c(44.999999999999993, 45), 8), c(3, 3))   # a bin centre on 8 bins
 expect_identical(oi(c(33.749999999999993, 33.75), 8), c(3, 3))
 expect_identical(oi(c(179.99999999999997, 0, 89.999999999998991, 90.000000000001037), 5), c(1, 1, 4, 4))
 expect_identical(oi(c(11.25 - 1e-6, 11.25 + 1e-6), 8), c(1, 2))
 # Segments at rounding distance from orientation, length and horizontal/vertical
 # edges; the fused kernel and the R binning of the published table agree.
 angles <- c(45, 33.75, 90, 22.5, 135, 67.5, 0, 157.5)
 eps <- c(-4, -1, 0, 1, 4) * .Machine$double.eps
 segs <- do.call(rbind, lapply(seq_along(angles), function(i) do.call(rbind, lapply(eps, function(e) {
  a <- angles[i] * (1 + e) * pi / 180; len <- 5 + i
  c(50, 50, 50 + len * cos(a), 50 - len * sin(a))
 }))))
 d <- mlv_drawing(lapply(seq_len(nrow(segs)), function(i) segs[i, , drop = FALSE]), 100, 100)
 st <- mlvtoolbox:::.segment_table(d)
 for (bins in c(2, 4, 5, 8, 13)) {
  fused <- contour_statistics(d, c("orientation", "horizontal_vertical", "curvature"), bins = bins,
                              by = "contour", estimator = "mlv2023")
  for (prop in c("orientation", "horizontal_vertical")) {
   values <- if (prop == "orientation") st$orientation else
    pmax(-1, pmin(1, abs(cos(st$orientation * pi / 180)) - abs(sin(st$orientation * pi / 180))))
   idx <- if (prop == "orientation") oi(values, bins) else mlvtoolbox:::.bin_index(values, seq(-1, 1, length.out = bins + 1L))
   r <- mlvtoolbox:::.grouped_bins(values, st$length, st$contour_index, seq(0, 1, length.out = bins + 1L), length(d$ids), idx)
   expect_identical(matrix(fused$value[fused$property == prop], ncol = bins, byrow = TRUE), r, info = paste(prop, bins))
  }
  # The R orientation path (turning specification) agrees with the fused one.
  spec <- feature_spec("orientation", bins = bins, normalize = FALSE)
  expect_equal(contour_statistics(d, spec = spec, by = "contour")$value,
               fused$value[fused$property == "orientation"], tolerance = 1e-12, info = bins)
 }
 # Mathematically equal orientations share a bin on every bin count.
 groups <- split(seq_len(nrow(st)), rep(seq_along(angles), each = length(eps)))
 for (bins in c(2, 4, 5, 8, 13)) for (g in groups) expect_length(unique(oi(st$orientation[g], bins)), 1)
})

random_drawing <- function(seed) {
 set.seed(seed)
 contours <- lapply(seq_len(sample(2:6, 1)), function(k) {
  n <- sample(1:12, 1)
  xy <- cbind(cumsum(c(sample(5:40, 1), sample(-6:6, n, TRUE))), cumsum(c(sample(5:40, 1), sample(-6:6, n, TRUE))))
  if (runif(1) < .3) xy <- xy + matrix(runif(length(xy), -.5, .5), ncol = 2)
  if (n > 2 && runif(1) < .3) xy[n + 1, ] <- xy[1, ]                       # closed
  if (runif(1) < .3) xy <- rbind(xy[1, ], (xy[1, ] + xy[2, ]) / 2, xy[-1, ]) # collinear vertex
  if (runif(1) < .2) xy <- rbind(xy, xy[nrow(xy) - 1, ])                    # exact retracing
  keep <- c(TRUE, rowSums(abs(diff(xy))) > 0); xy <- xy[keep, , drop = FALSE]
  if (nrow(xy) < 2) xy <- rbind(xy, xy[1, ] + c(3, 4))
  cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE])
 })
 mlv_drawing(contours, 120, 120)
}
edge_drawings <- function() list(
 single = mlv_drawing(list(rbind(c(1, 1, 11, 1))), 20, 20),
 collinear = mlv_drawing(list(rbind(c(1, 1, 5, 1), c(5, 1, 11, 1))), 20, 20),
 spur = mlv_drawing(list(rbind(c(1, 1, 11, 1), c(11, 1, 6, 1))), 20, 20),
 there_and_back = mlv_drawing(list(rbind(c(1, 1, 11, 1), c(11, 1, 1, 1))), 20, 20),
 square_mid_edge = mlv_drawing(list(rbind(c(5, 1, 9, 1), c(9, 1, 9, 9), c(9, 9, 1, 9), c(1, 9, 1, 1), c(1, 1, 5, 1))), 20, 20),
 triangle = mlv_drawing(list(rbind(c(1, 1, 11, 1), c(11, 1, 1, 11), c(1, 11, 1, 1))), 20, 20),
 uneven_circle = circle_drawing(32, c(1, 2)),
 empty = mlv_drawing(list(), 20, 20))

test_that("the native turning histogram equals the exact integral of the R profile", {
 r_histogram <- function(d, bins, range, by) {
  edges <- 10^seq(log10(range[1] + 1), log10(range[2] + 1), length.out = bins + 1L) - 1
  ng <- if (by == "drawing") 1L else length(d$ids)
  p <- curvature_profile(d, "turning")
  groups <- if (by == "drawing") rep.int(1L, nrow(p)) else p$contour_index
  mlvtoolbox:::.linear_histogram(p$curvature_start, p$curvature_end, p$length, groups, edges, ng)
 }
 drawings <- c(lapply(1:40, random_drawing), edge_drawings())
 for (k in seq_along(drawings)) for (bins in c(3, 8)) for (by in c("drawing", "contour")) {
  d <- drawings[[k]]; range <- if (bins == 3) c(0, 90) else c(.5, 40)
  native <- contour_statistics(d, "curvature", bins = bins, by = by, curvature_range = range)
  expected <- r_histogram(d, bins, range, by)
  # Same values and bins (identical arithmetic and snapping); weights are summed
  # per whole dual cell natively but per half-cell in the profile, so sums can
  # differ in the last bits.
  expect_equal(matrix(native$value, ncol = bins, byrow = TRUE), unname(expected),
               tolerance = 1e-12, info = paste(names(drawings)[k], k, bins, by))
 }
 # Exact bin agreement on values that sit on bin edges up to rounding.
 tri <- edge_drawings()$triangle
 h <- contour_statistics(tri, "curvature", bins = 8)
 expect_equal(sum(h$value), sum(as.data.frame(tri)$length))
})

test_that("segment tables have one row per segment with exact mean curvature", {
 independent <- function(d, estimator, scale) {
  out <- numeric(); total <- numeric()
  for (i in seq_along(d$contours)) {
   z <- unname(d$contours[[i]]); len <- sqrt((z[, 3] - z[, 1])^2 + (z[, 4] - z[, 2])^2)
   cp <- mlvtoolbox:::.contour_profile(d, i, estimator, scale, 1); p <- cp$profile; L <- cp$total
   s_in <- c(0, cumsum(len)) - c(0, cumsum(len))[cp$curve$kept[1]]
   shifts <- if (cp$curve$closed) c(-L, 0, L) else 0
   for (k in seq_along(len)) {
    acc <- 0
    for (j in seq_len(nrow(p))) for (sh in shifts) {
     a <- p$s_start[j] + sh; b <- p$s_end[j] + sh
     lo <- max(a, s_in[k]); hi <- min(b, s_in[k + 1])
     if (hi > lo) {
      at <- function(t) p$curvature_start[j] + (p$curvature_end[j] - p$curvature_start[j]) * (t - a) / (b - a)
      acc <- acc + (hi - lo) * (at(lo) + at(hi)) / 2
     }
    }
    out <- c(out, acc / len[k])
   }
   total <- c(total, sum(p$curvature * p$length))
  }
  list(mean = out, total = total)
 }
 drawings <- c(lapply(41:70, random_drawing), edge_drawings()[-8])
 for (k in seq_along(drawings)) {
  d <- drawings[[k]]; n <- sum(vapply(d$contours, nrow, 1L))
  for (e in c("turning", "arclength")) {
   if (e == "arclength" && inherits(try(curvature_profile(d, e), silent = TRUE), "try-error")) next
   s <- if (e == "turning") as.data.frame(d) else as.data.frame(d, spec = feature_spec(estimator = e, curvature_scale = 1.5))
   ref <- independent(d, e, 1.5)
   expect_identical(nrow(s), n)
   expect_identical(s[, c("x1", "y1", "x2", "y2", "orientation")],
                    as.data.frame(d, estimator = "mlv2023")[, c("x1", "y1", "x2", "y2", "orientation")])
   expect_equal(s$curvature, ref$mean, tolerance = 1e-10, info = paste(k, e))
   per_contour <- tapply(s$curvature * s$length, factor(s$contour_index, seq_along(d$ids)), sum)
   expect_equal(unname(as.vector(per_contour)), ref$total, tolerance = 1e-10, info = paste(k, e))
  }
 }
 # Closed contour starting in the middle of an edge, and exact retracing.
 sq <- as.data.frame(edge_drawings()$square_mid_edge)
 expect_equal(sq$curvature * sq$length, c(45, 90, 90, 90, 45))
 spur <- as.data.frame(edge_drawings()$spur)
 expect_equal(spur$curvature * spur$length, c(120, 60))
 expect_true(all(is.na(spur$signed_curvature)))
 roof <- as.data.frame(mlv_drawing(list(rbind(c(5, 20, 15, 10), c(15, 10, 25, 20))), 30, 30))
 expect_equal(roof$signed_curvature, -roof$curvature)
})

test_that("dual cells narrower than 64 ulp of the contour length keep their turn", {
 seg <- function(xy) cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE])
 zig <- cbind(rep(c(100, 2), 500), seq(60, 90, length.out = 1000))
 xy <- rbind(c(1, 1), c(1, 1 + 1e-11), c(1 + 1e-11, 1 + 1e-11), c(1 + 1e-11, 50), zig)
 d <- mlv_drawing(list(seg(xy)), 200, 200)
 exact <- sum(mlvtoolbox:::.canonical_curve(d$contours[[1]], TRUE)$absolute)
 expect_equal(sum(with(as.data.frame(d), curvature * length)), exact, tolerance = 1e-12)
 expect_equal(sum(with(curvature_profile(d, "turning"), curvature * length)), exact, tolerance = 1e-12)
 native <- contour_statistics(d, "curvature", bins = 4, curvature_range = c(0, 90), length_range = c(2, 1e6))
 specified <- contour_statistics(d, spec = feature_spec("curvature", bins = 4, normalize = FALSE, length_range = c(2, 1e6)))
 expect_equal(native$value, specified$value, tolerance = 1e-12)
 # the 90 deg corner's 1e-11 px dual cell, measured at arc position ~1 (ulp 2e-16)
 expect_equal(native$value[4], 1e-11, tolerance = 1e-4)
})

test_that("compute_contour_properties(estimator =) is the drawing's default thereafter", {
 a <- rbind(c(0, 0, 10, 0), c(10, 0, 10, 10))                              # 90 deg: turning 90, mlv2023 180
 b <- rbind(c(0, 30, 20 / 3, 30), c(20 / 3, 30, 20 / 3 + 20 / 3 * cos(50 * pi / 180), 30 - 20 / 3 * sin(50 * pi / 180)))
 b <- rbind(b, c(b[2, 3:4], b[2, 3] + 20 / 3, b[2, 4]))                      # turns 50, 50: turning 100, mlv2023 150
 plain <- mlv_drawing(list(a + 20, b + 20), 80, 80, ids = c("a", "b"))
 old <- compute_contour_properties(plain, estimator = "mlv2023")
 new <- compute_contour_properties(plain)
 expect_identical(old$curvature_estimator, "mlv2023"); expect_identical(new$curvature_estimator, "turning")
 published <- function(f, ...) f(plain, ..., estimator = "mlv2023")
 expect_identical(as.data.frame(old), published(as.data.frame))
 expect_identical(old$segment_properties, published(as.data.frame))
 expect_identical(contour_statistics(old), published(contour_statistics))
 expect_identical(feature_table(old), published(feature_table))
 expect_identical(average_property(old, "curvature"), published(average_property, "curvature"))
 model <- bind_feature_model(function(f) f$curvature_1, published(feature_table))
 expect_equal(predict_contours(old, model), published(predict_contours, model))
 expect_error(predict_contours(new, model), "schema|estimator|incompatible")
 # Explicit arguments still win.
 expect_identical(contour_statistics(old, estimator = "turning"), contour_statistics(plain))
 expect_identical(as.data.frame(old, estimator = "turning"), as.data.frame(plain))
 # Densities, plots and contour splits use the recorded estimator.
 expect_equal(sum(feature_density(old, "curvature")), 180 + 150, tolerance = 1e-9)
 expect_equal(sum(feature_density(new, "curvature")), 90 + 100, tolerance = 1e-9)
 expect_identical(split_contours(old, "integrated_curvature")$scores, c(a = 2, b = 1))
 expect_identical(split_contours(new, "integrated_curvature")$scores, c(a = 1, b = 2))
 w <- c(curvature_8 = 1, curvature_5 = 1, curvature_6 = 1, curvature_7 = 1)
 h <- published(feature_table, c("orientation", "length", "curvature"), normalize = FALSE)
 expect_equal(split_by_histogram(old, w)$scores, setNames(as.numeric(as.matrix(h[, names(w)]) %*% w), c("a", "b")))
 pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
 expect_identical(plot(old, property = "curvature"), old)
 expect_identical(mlvtoolbox:::.measurement_table(old), published(as.data.frame))
 # A specification governs instead; geometry operations keep the recipe, not derived data.
 expect_null(compute_contour_properties(old, spec = feature_spec())$curvature_estimator)
 moved <- rotate_drawing(old, 30)
 expect_identical(moved$curvature_estimator, "mlv2023"); expect_null(moved$segment_properties)
 expect_identical(subset_contours(old, "a")$curvature_estimator, "mlv2023")
})
