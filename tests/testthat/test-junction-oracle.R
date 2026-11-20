# Junction regressions with analytical answers, derived by hand from the
# construction (no implementation output is used as a reference). Angles are
# counterclockwise degrees with image y down. See
# docs/evidence/mlvref/junctions/REPORT.md for the natural-data adjudication.
jd <- function(...) mlv_drawing(list(...), 100000, 100000, zero_length = "drop")
sg <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)
gaps <- function(j, i = 1L) sort(j$angles[[i]])

test_that("arm directions follow the arms, not the merged junction position (MLV-R005)", {
  # Three straight rays whose pairwise intersections are (20.6,30), (20,30)
  # and (20.6,30.6); they merge at the centroid (20.4,30.2). The rays point at
  # 0, 270 and 135 degrees, so the arm gaps are exactly 90, 135 and 135.
  j <- detect_junctions(jd(sg(20, 30, 50, 30), sg(20.6, 30.5, 20.6, 60), sg(20, 30, 0, 10)))
  expect_equal(nrow(j), 1L)
  expect_equal(c(j$x, j$y), c(20.4, 30.2), tolerance = 1e-12)
  expect_equal(j$type, "Y")
  expect_equal(gaps(j), c(90, 135, 135), tolerance = 1e-12)
  # Three lines crossing pairwise at (20,30), (21,30) and (20,31): six arms at
  # 0, 45, 90, 180, 225 and 270 degrees, whatever the merged position.
  j <- detect_junctions(jd(sg(0, 30, 40, 30), sg(20, 10, 20, 50), sg(1, 50, 41, 10)))
  expect_equal(c(j$x, j$y), c(61, 91) / 3, tolerance = 1e-12)
  expect_equal(j$type, "Star")
  expect_equal(gaps(j), c(45, 45, 45, 45, 90, 90), tolerance = 1e-12)
})

test_that("merged parallel crossings keep exact arm directions and are order invariant", {
  # Crossings at x = 10, 13 and 11.5 chain within the 2 px merge radius into
  # one junction at x = 11.5 with arms 0, 180 and three times 90 and 270.
  bar <- sg(0, 10, 30, 10)
  stems <- list(sg(10, 0, 10, 20), sg(13, 0, 13, 20), sg(11.5, 0, 11.5, 20))
  for (order in list(1:3, c(1, 3, 2), c(3, 2, 1))) {
    j <- detect_junctions(do.call(jd, c(list(bar), stems[order])))
    expect_equal(nrow(j), 1L)
    expect_equal(c(j$x, j$y), c(11.5, 10), tolerance = 1e-12)
    expect_equal(j$type, "Star")
    expect_equal(gaps(j), c(0, 0, 0, 0, 90, 90, 90, 90), tolerance = 1e-12)
  }
})

test_that("arms are the directions in which contours leave the neighbourhood (MLV-R006)", {
  bar <- sg(0, 10, 10, 10, 10, 10, 18, 4)
  bend <- atan2(6, 8) * 180 / pi # bar arms at bend and 180 degrees
  # A 3 px contour centred on the bar vertex never leaves the 2 px disc, so
  # only the two bar arms remain, in either direction: no junction.
  expect_equal(nrow(detect_junctions(jd(bar, sg(10, 8.5, 10, 11.5)))), 0L)
  expect_equal(nrow(detect_junctions(jd(bar, sg(10, 11.5, 10, 8.5)))), 0L)
  # A 6 px crossing leaves it twice: X with arms at bend, 90, 180 and 270.
  j <- detect_junctions(jd(bar, sg(10, 7, 10, 13)))
  expect_equal(j$type, "X")
  expect_equal(gaps(j), c(90 - bend, 90, 90, 90 + bend), tolerance = 1e-12)
  # One polyline split into three contours with a short middle piece is a
  # continuation, not a junction (two arms leave the merged neighbourhood).
  expect_equal(nrow(detect_junctions(jd(sg(0, 10, 10, 10), sg(10, 10, 11, 11), sg(11, 11, 20, 11)))), 0L)
  # A 1 px overshoot stays inside the neighbourhood: a T with one stem arm.
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 9, 10, 10.5, 10, 10.5, 10, 20)))
  expect_equal(j$type, "T")
  expect_equal(gaps(j), c(90, 90, 180))
  # A straight line with vertices 1 px either side of a crossing passes once:
  # two arms, not one per segment inside the neighbourhood.
  j <- detect_junctions(jd(sg(0, 10, 9, 10, 9, 10, 11, 10, 11, 10, 20, 10), sg(10, 0, 10, 20)))
  expect_equal(j$type, "X")
  expect_equal(gaps(j), rep(90, 4))
  # With a zero radius only exact contour ends terminate an arm.
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)), endpoint_radius = 0)
  expect_equal(j$type, "T")
  expect_equal(gaps(j), c(90, 90, 180))
})

test_that("reversing or renumbering contours leaves junctions unchanged", {
  rays <- lapply(c(90, 210, 330), function(a) sg(10, 10, 10 + 8 * cospi(a / 180), 10 - 8 * sinpi(a / 180)))
  a <- detect_junctions(do.call(jd, rays))
  b <- detect_junctions(do.call(jd, lapply(rev(rays), function(z) z[, c(3, 4, 1, 2), drop = FALSE])))
  expect_equal(a$type, "Y"); expect_equal(b$type, "Y")
  expect_equal(gaps(a), c(120, 120, 120), tolerance = 1e-12)
  expect_equal(gaps(b), gaps(a), tolerance = 1e-12)
  t1 <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)))
  t2 <- detect_junctions(jd(sg(10, 20, 10, 10), sg(20, 10, 0, 10)))
  expect_equal(c(t1$type, t2$type), c("T", "T"))
  expect_equal(gaps(t1), c(90, 90, 180)); expect_equal(gaps(t2), c(90, 90, 180))
})

test_that("the attachment end is the nearer endpoint, not the start point (MLV-M022)", {
  # A stem 0.2 px short of the bar, with a 1 px first piece inside the 2 px
  # neighbourhood. Drawn either way it is a T stem pointing away (270).
  bar <- sg(0, 10, 20, 10)
  forward <- sg(10, 10.2, 10, 11.2, 10, 11.2, 10, 20)
  backward <- sg(10, 20, 10, 11.2, 10, 11.2, 10, 10.2)
  for (stem in list(forward, backward)) {
    j <- detect_junctions(jd(bar, stem))
    expect_equal(c(j$x, j$y), c(10, 10))
    expect_equal(j$type, "T")
    expect_equal(gaps(j), c(90, 90, 180))
  }
})

test_that("intersections near the end of a long segment are exact (MLV-M019)", {
  # Parameter 3/40000 < 1e-4 along the bar: MLV snaps to (0,5), 3 px away.
  j <- detect_junctions(jd(sg(0, 5, 40000, 5), sg(3, 0, 3, 10)))
  expect_equal(c(j$x, j$y), c(3, 5), tolerance = 1e-9)
  expect_equal(j$type, "X")
  expect_equal(gaps(j), rep(90, 4), tolerance = 1e-9)
})

test_that("zero-length segments add no arm (MLV-M020)", {
  d <- jd(sg(0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 10), sg(10, 10, 10, 20))
  j <- detect_junctions(d)
  expect_equal(j$type, "T")
  expect_equal(gaps(j), c(90, 90, 180))
})

test_that("corners, collinear contacts and overlaps are not junctions", {
  expect_equal(nrow(detect_junctions(jd(sg(0, 10, 10, 10), sg(10, 10, 10, 0)))), 0L)
  expect_equal(nrow(detect_junctions(jd(sg(0, 10, 10, 10), sg(10, 10, 20, 10)))), 0L)
  expect_equal(nrow(detect_junctions(jd(sg(0, 10, 12, 10), sg(8, 10, 20, 10)))), 0L)
  psi <- detect_junctions(jd(sg(10, 0, 10, 20), sg(10, 10, 10 + 8 * cospi(.25), 10 - 8 * sinpi(.25)),
                             sg(10, 10, 10 + 8 * cospi(.75), 10 - 8 * sinpi(.75))))
  expect_equal(psi$type, "X")
  expect_equal(gaps(psi), c(45, 45, 135, 135), tolerance = 1e-12)
})

test_that("the default Euclidean gap is rotation invariant; Chebyshev is optional (MLV-M018)", {
  # A 45-degree stem whose end lies 1.2 px (Euclidean) from the bar: within
  # one raster step (Chebyshev 0.85 px), beyond a 1 px Euclidean gap.
  bar <- sg(0, 10, 20, 10)
  stem <- function(g) sg(10 + g / sqrt(2), 10 + g / sqrt(2), 18, 18)
  expect_equal(nrow(detect_junctions(jd(bar, stem(1.2)))), 0L) # default: Euclidean
  expect_equal(nrow(detect_junctions(jd(bar, stem(1.2)), gap_metric = "euclidean")), 0L)
  expect_equal(detect_junctions(jd(bar, stem(1.2)), gap_metric = "chebyshev")$type, "T")
  j <- detect_junctions(jd(bar, stem(.9)))
  expect_equal(c(j$x, j$y), c(10, 10), tolerance = 1e-12)
  expect_equal(gaps(j), c(45, 135, 180), tolerance = 1e-12)
  # Rotating a T with a 0.95 px (found) or 1.05 px (not found) gap.
  for (g in c(.95, 1.05)) for (angle in c(0, 17, 45, 71)) {
    d <- rotate_drawing(mlv_drawing(list(sg(-10, 0, 10, 0), sg(0, g, 0, 10)), 100, 100), angle, center = c(0, 0))
    expect_equal(nrow(detect_junctions(d)), as.integer(g < 1))
  }
  expect_error(detect_junctions(jd(bar, stem(1)), gap_metric = "manhattan"))
})

test_that("a segment beyond the junction has no arm pointing back across it (MLV-R008)", {
  # Stem 0.8 px short of the bar, endpoint radius 0.5 or 0: the stem stays
  # outside the neighbourhood, wholly below the junction, so it has one arm
  # (270) and the junction is a T, not an X.
  for (r in c(.5, 0)) {
    j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10.8, 10, 20)), endpoint_radius = r)
    expect_equal(c(j$x, j$y), c(10, 10))
    expect_equal(j$type, "T")
    expect_equal(gaps(j), c(90, 90, 180))
  }
  # Merged chain: crossings at x = 10, 11.9, 13.8 and a T stem starting at
  # x = 15.7, 2.85 px from the centroid 12.85: 2 + 3 * 2 + 1 = 9 arms.
  j <- detect_junctions(jd(sg(0, 10, 30, 10), sg(10, 0, 10, 20), sg(11.9, 0, 11.9, 20),
                           sg(13.8, 0, 13.8, 20), sg(15.7, 10, 15.7, 20)))
  expect_equal(c(j$x, j$y), c(12.85, 10), tolerance = 1e-12)
  expect_equal(gaps(j), c(rep(0, 5), rep(90, 4)), tolerance = 1e-12)
  # The reversed stem (drawn toward the bar) gives the same junction.
  k <- detect_junctions(jd(sg(0, 10, 30, 10), sg(10, 0, 10, 20), sg(11.9, 0, 11.9, 20),
                           sg(13.8, 0, 13.8, 20), sg(15.7, 20, 15.7, 10)))
  expect_equal(gaps(k), gaps(j), tolerance = 1e-12)
})

test_that("closed contours give the same junction for every start vertex (MLV-R007)", {
  # The polygon passes 1 px below the bar crossing at (10,10): arms up (90)
  # and along (10,11) -> (20,20), i.e. 360 - atan2(9,10) degrees, plus the
  # bar arms 0 and 180.
  v <- rbind(c(10, 11), c(10, 0), c(30, 0), c(20, 20))
  pb <- atan2(9, 10) * 180 / pi
  for (k in 0:3) for (rev in c(FALSE, TRUE)) {
    w <- v[(seq_len(4) + k - 1L) %% 4L + 1L, ]
    if (rev) w <- w[4:1, ]
    w <- rbind(w, w[1, ])
    j <- detect_junctions(jd(sg(0, 10, 20, 10), cbind(w[-5, ], w[-1, ])))
    expect_equal(j$type, "X")
    expect_equal(gaps(j), c(pb, 90, 90, 180 - pb), tolerance = 1e-12)
  }
  # A closed contour lying wholly inside the neighbourhood has no arm.
  tiny <- rbind(c(10, 9), c(11, 10), c(10, 11), c(9, 10)); tiny <- rbind(tiny, tiny[1, ])
  expect_equal(nrow(detect_junctions(jd(sg(0, 10, 20, 10), cbind(tiny[-5, ], tiny[-1, ])))), 0L)
})

test_that("sweep and exhaustive candidate search agree for both gap metrics", {
  set.seed(2611)
  for (i in seq_len(60)) {
    d <- mlv_drawing(lapply(seq_len(sample(5:40, 1)), function(k) {
      p <- runif(2, 0, 30); m <- sample(1:3, 1)
      xy <- rbind(p, p + apply(matrix(rnorm(2 * m, 0, 4), m, 2), 2, cumsum))
      cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE])
    }), 40, 40)
    for (g in c("chebyshev", "euclidean")) {
      a <- detect_junctions(d, gap_metric = g, method = "sweep")
      b <- detect_junctions(d, gap_metric = g, method = "exhaustive")
      expect_equal(a[, 1:6], b[, 1:6], tolerance = 1e-12)
      expect_equal(unclass(a$angles), unclass(b$angles), tolerance = 1e-12)
    }
  }
})

test_that("rule = 'mlv2023' reproduces pinned MLV exactly on the synthetic suite", {
  # inst/fixtures/octave/junction_*.csv: geometry (zero-length segments
  # removed, as an mlv_drawing requires) and the output of the pinned MATLAB
  # source (f12418a) under Octave 11.1; tools/mlvref/synthetic_junctions.R.
  geom <- utils::read.csv(system.file("fixtures", "octave", "junction_cases.csv", package = "mlvtoolbox"))
  ref <- utils::read.csv(system.file("fixtures", "octave", "junction_octave.csv", package = "mlvtoolbox"),
                         colClasses = c(contours = "character", segments = "character", angles = "character"))
  num <- function(s) as.numeric(strsplit(s, ";", fixed = TRUE)[[1]])
  expect_equal(length(unique(geom$case)), 41L)
  for (cs in unique(geom$case)) {
    g <- geom[geom$case == cs, ]
    d <- mlv_drawing(lapply(split(g, g$contour), function(z) as.matrix(z[, c("x1", "y1", "x2", "y2")])), 100000, 100000)
    j <- detect_junctions(d, rule = "mlv2023")
    o <- ref[ref$case == cs & !ref$type %in% c("none", "error"), ]
    expect_equal(nrow(j), nrow(o), info = cs)
    for (k in seq_len(nrow(o))) {
      i <- which(abs(j$x - o$x[k]) <= 1e-9 & abs(j$y - o$y[k]) <= 1e-9)
      expect_length(i, 1L)
      expect_equal(j$type[i], o$type[k], info = cs)
      expect_equal(j$angles[[i]], num(o$angles[k]), tolerance = 1e-12, info = cs)
      expect_equal(as.numeric(j$members[[i]]$contour_index), num(o$contours[k]), info = cs)
      expect_equal(as.numeric(j$members[[i]]$segment_index), num(o$segments[k]), info = cs)
    }
  }
  expect_error(detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)), rule = "mlv2023", gap_metric = "euclidean"),
               "chebyshev")
  expect_equal(detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)), rule = "mlv2023", gap_metric = "chebyshev")$type, "T")
})

test_that("near-zero-length segments are zero-length under the connectivity tolerance (MLV-R009)", {
  # A segment of length about 2.5e-13 at coordinate 871 (2 and 1 ulp): its
  # endpoints differ only by rounding residue (<= 64 ulp of 871).
  z <- rbind(c(879, 951, 871, 946), c(871, 946, 871 - 2^-42, 946 - 2^-43))
  expect_error(mlv_drawing(list(z), 1000, 1000), "zero-length")
  d <- mlv_drawing(list(z), 1000, 1000, zero_length = "drop")
  expect_equal(nrow(d$contours[[1]]), 1L)
  expect_equal(as.data.frame(d)$curvature, 0)
  # Genuine short segments are kept, at any scale.
  expect_equal(nrow(mlv_drawing(list(rbind(c(0, 0, 1e-6, 0))), 10, 10)$contours[[1]]), 1L)
  expect_equal(nrow(mlv_drawing(list(rbind(c(0, 0, 1e-15, 0))), 10, 10)$contours[[1]]), 1L)
  expect_equal(nrow(mlv_drawing(list(rbind(c(871, 946, 871 - 1e-10, 946))), 1000, 1000)$contours[[1]]), 1L)
})

test_that("rule = 'mlv2023' visits merge partners in MATLAB order and handles long chains", {
  # Octave (pinned f12418a): contours 1;2;3;5;4;5;5, segments 1;1;1;3;1;1;2,
  # angles 0;90;0;90;0;90;90 (also in the synthetic fixture as merge_order).
  d <- jd(sg(9, 11.2, 11, 11.2), sg(10, 11, 10, 12.5), sg(11.5, 9, 11.5, 10.5), sg(8.5, 9, 8.5, 10.5),
          sg(5, 10, 9.5, 10, 9.5, 10, 10.5, 10, 10.5, 10, 20, 10))
  j <- detect_junctions(d, rule = "mlv2023")
  expect_equal(j$members[[1]]$contour_index, c(1, 2, 3, 5, 4, 5, 5))
  expect_equal(j$members[[1]]$segment_index, c(1, 1, 1, 3, 1, 1, 2))
  expect_equal(j$angles[[1]], c(0, 90, 0, 90, 0, 90, 90))
  # Every small contour lies within 2 px of a merged crossing; only the bar
  # leaves the neighbourhood (two arms), so the default rule finds nothing.
  expect_equal(nrow(detect_junctions(d)), 0L)
  # 3,000 crossings 1.5 px apart form one merge chain: one junction whose
  # arms are the bar's two plus two per stem, under both rules.
  xs <- 10 + 1.5 * (0:2999)
  chain <- do.call(jd, c(list(sg(0, 30, max(xs) + 10, 30)), lapply(xs, function(x) sg(x, 0, x, 60))))
  for (rule in c("mlv2023", "exit")) {
    j <- detect_junctions(chain, rule = rule)
    expect_equal(nrow(j), 1L)
    expect_length(j$angles[[1]], 6002L)
  }
})

test_that("merged junctions keep the arms of every merged crossing (MLV-R014)", {
  # Five crossings 2 px apart merge into one junction at x = 14; the outer
  # stems lie 4 px from that position. Each crossing alone is an X (bar 0/180,
  # stem 90/270); merged, the arms are their union: bar twice, stems 5 x 2.
  bar <- sg(0, 10, 30, 10)
  stems <- lapply(seq(10, 18, by = 2), function(x) sg(x, 0, x, 20))
  j <- detect_junctions(do.call(jd, c(list(bar), stems)))
  expect_equal(nrow(j), 1L)
  expect_equal(c(j$x, j$y), c(14, 10))
  expect_equal(sort(j$angles[[1]]), c(rep(0, 8), rep(90, 4)))
  separate <- detect_junctions(do.call(jd, c(list(bar), stems)), merge_radius = 1)
  expect_equal(nrow(separate), 5L)
  expect_true(all(separate$type == "X"))
})

test_that("an extension exactly at the limit is accepted in either direction (MLV-R015)", {
  # Stem length 3.5 ends 1 px from the bar: limit min(1, 0.3 * 3.5) = 1.
  for (stem in list(sg(3, 6, 3, 2.5), sg(3, 2.5, 3, 6))) {
    j <- detect_junctions(jd(sg(0, 1.5, 6, 1.5), stem))
    expect_equal(j$type, "T"); expect_equal(c(j$x, j$y), c(3, 1.5))
  }
  # Gap equal to the limit found by sweep and exhaustive search alike.
  d <- mlv_drawing(list(sg(0, 0, 0, 6), sg(1, 3, 4.5, 3)), 10, 10)
  expect_equal(nrow(detect_junctions(d, method = "exhaustive")), 1L)
  expect_identical(detect_junctions(d), detect_junctions(d, method = "exhaustive"))
})

test_that("grid drawings with exact ties: sweep == exhaustive bitwise, reversal invariant", {
  set.seed(4401)
  grid_drawing <- function() lapply(seq_len(sample(4:25, 1)), function(k) {
    p <- sample(0:40, 2) / 2; m <- sample(1:3, 1)
    xy <- rbind(p, p + apply(matrix(sample(-8:8, 2 * m, TRUE) / 2, m, 2), 2, cumsum))
    xy <- xy[c(TRUE, rowSums(abs(diff(xy))) > 0), , drop = FALSE]
    if (nrow(xy) < 2) xy <- rbind(p, p + c(1, 0))
    cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE])
  })
  # Reversal permutes segment indices, so merged positions may differ in the
  # last place (summation order); compare at 1e-9.
  key <- function(j) if (!nrow(j)) character() else
    sort(paste(round(j$x, 9), round(j$y, 9), j$type,
               vapply(j$angles, function(a) paste(round(sort(a), 9), collapse = "/"), "")))
  for (i in seq_len(150)) {
    cs <- grid_drawing(); d <- mlv_drawing(cs, 50, 50)
    for (g in c("euclidean", "chebyshev"))
      expect_identical(detect_junctions(d, gap_metric = g), detect_junctions(d, gap_metric = g, method = "exhaustive"))
    rev <- mlv_drawing(lapply(cs, function(z) z[nrow(z):1, c(3, 4, 1, 2), drop = FALSE]), 50, 50)
    expect_identical(key(detect_junctions(rev)), key(detect_junctions(d)))
  }
})

test_that("dropping residue segments keeps contours connected and closed (MLV-R016)", {
  e <- 1e-14; line <- sg(0, 15, 15, 15)
  starts <- list(sg(10, 10, 20, 10, 20, 10, 20, 20, 20, 20, 10, 20, 10, 20, 10, 10 + e, 10, 10 + e, 10, 10),
                 sg(20, 10, 20, 20, 20, 20, 10, 20, 10, 20, 10, 10 + e, 10, 10 + e, 10, 10, 10, 10, 20, 10))
  for (poly in starts) {
    d <- mlv_drawing(list(poly, line), 40, 40, zero_length = "drop")
    z <- d$contours[[1]]
    expect_equal(nrow(z), 4L)
    expect_identical(unname(z[1, 1:2]), unname(z[4, 3:4])) # exactly closed
    expect_identical(unname(z[-1, 1:2]), unname(z[-4, 3:4])) # exactly connected
    j <- detect_junctions(d) # the line crosses the left edge at (10,15)
    expect_equal(j$type, "X"); expect_equal(sort(j$angles[[1]]), rep(90, 4))
  }
  # Consecutive residue segments accumulate beyond 64 ulp but still clean.
  z <- mlv_drawing(list(sg(0, 1, 1, 1, 1, 1, 1 + 1e-14, 1, 1 + 1e-14, 1, 1 + 2e-14, 1, 1 + 2e-14, 1, 5, 1)),
                   10, 10, zero_length = "drop")$contours[[1]]
  expect_equal(unname(z), sg(0, 1, 1, 1, 1, 1, 5, 1))
})

test_that("residue closing gaps are exact without any dropped segment (MLV-R016)", {
  # The polygon closes up to 1e-14 at its start vertex and has no residue
  # segment to drop; both start vertices are the same closed square, which the
  # line y = 15 crosses at its left edge, so the junction is an X either way.
  e <- 1e-14; line <- sg(0, 15, 15, 15)
  starts <- list(sg(10, 10, 20, 10, 20, 10, 20, 20, 20, 20, 10, 20, 10, 20, 10, 10 + e),
                 sg(20, 10, 20, 20, 20, 20, 10, 20, 10, 20, 10, 10 + e, 10, 10, 20, 10))
  for (poly in starts) {
    d <- mlv_drawing(list(poly, line), 40, 40)
    z <- d$contours[[1]]; n <- nrow(z)
    expect_identical(unname(z[1, 1:2]), unname(z[n, 3:4]))
    expect_identical(unname(z[-1, 1:2]), unname(z[-n, 3:4]))
    j <- detect_junctions(d)
    expect_equal(j$type, "X"); expect_equal(sort(j$angles[[1]]), rep(90, 4))
  }
  # A genuinely open contour is left open.
  open <- mlv_drawing(list(sg(0, 0, 5, 0, 5, 0, 5, 5)), 10, 10)$contours[[1]]
  expect_false(identical(unname(open[1, 1:2]), unname(open[2, 3:4])))
})
