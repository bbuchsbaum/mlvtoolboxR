# Undefined branch scores are NA under the default branch rule
# (definition "direction-symmetric-tiemax-graph-defined-3"), not the published zeros.

P <- c("parallelism", "separation", "mirror", "taper")
line_mat <- function(cols, row = 10L, radius = 4, rule = "invariant", h = 20L, w = 40L) {
  sk <- matrix(FALSE, h, w); sk[row, cols] <- TRUE
  m <- structure(list(skeleton = sk, distance_map = matrix(radius, h, w), AOF = matrix(0, h, w),
                      mask = matrix(FALSE, h, w), branches = list(),
                      parameters = list(branch_rule = rule, min_branch_points = 10L)), class = "mlv_mat")
  m$branches <- trace_skeleton(m, 10L, rule = "invariant"); m   # same pixels for both rules
}

test_that("a straight branch of constant radius scores its end points like its interior", {
  # Hand derivation: every interior parallelism/taper ratio is 1 (no radius
  # change), so the end points, whose own windows are empty, take the mean of
  # their one defined neighbour: 1. The published estimator scores the ends 0
  # before the width-3 mean: (0 + 1) / 2 and (0 + 1 + 1) / 3, to the power 10.
  m <- compute_mat_properties(line_mat(5:20), P)
  b <- m$branches[[1]]; expect_identical(nrow(b), 16L)
  for (p in c("parallelism", "taper")) expect_identical(b[[p]], rep(1, 16))
  expect_equal(b$separation, rep((1 - 1 / 4)^10, 16))
  expect_equal(b$mirror, rep((1 - 1 / 16)^10, 16))
  mc <- compute_mat_properties(line_mat(5:20, rule = "mlv2023"), P)$branches[[1]]
  expect_equal(mc$parallelism[c(1, 2, 3, 16)], c(.5^10, (2 / 3)^10, 1, .5^10))
  # the public function keeps the published values
  expect_equal(score_branch(b, "parallelism")[1:2], c(.5^10, (2 / 3)^10))
  expect_identical(score_branch(b, "parallelism"), mc$parallelism)
})

test_that("two- and three-point branches are NA where the estimator is undefined, not 0", {
  for (n in 2:3) {
    m <- compute_mat_properties(line_mat(5:(4 + n)), P)
    b <- m$branches[[1]]; expect_identical(nrow(b), n)
    expect_true(all(is.na(b$mirror)))                                 # mirror needs more than 3 points
    expect_equal(b$separation, rep((1 - 1 / 4)^10, n))
    if (n == 2L) { expect_true(all(is.na(b$parallelism))); expect_true(all(is.na(b$taper))) }
    else { expect_identical(b$parallelism, rep(1, 3)); expect_identical(b$taper, rep(1, 3)) }
    expect_true(all(is.na(m$score_images$mirror)))
    # the published rule scores the same branches 0
    mc <- compute_mat_properties(line_mat(5:(4 + n), rule = "mlv2023"), P)$branches[[1]]
    expect_identical(mc$mirror, rep(0, n))
    if (n == 2L) expect_identical(mc$parallelism, rep(0, 2))
  }
})

test_that("a short component kept by the tracer maps NA unless another branch covers the pixel", {
  # A two-pixel skeleton between two short lines: separation is defined (radius
  # 3), parallelism is not, so the contours get separation values and NA parallelism.
  d <- mlv_drawing(list(matrix(c(8, 7, 13, 7), 1), matrix(c(8, 13, 13, 13), 1)), 20, 20)
  sk <- matrix(FALSE, 20, 20); sk[10, 10:11] <- TRUE
  m <- structure(list(skeleton = sk, distance_map = matrix(3, 20, 20), AOF = matrix(0, 20, 20),
                      mask = rasterize(d), branches = list(),
                      parameters = list(branch_rule = "invariant", min_branch_points = 10L)), class = "mlv_mat")
  m$branches <- trace_skeleton(m, 10L); expect_length(m$branches, 1L)
  y <- map_mat_properties(d, compute_mat_properties(m, P))
  ink <- rasterize(d)
  expect_true(any(is.finite(y$mat_images$separation[ink])))
  for (p in c("parallelism", "mirror", "taper")) expect_true(all(is.na(y$mat_images[[p]][ink])))
  # Add a long branch whose tangent samples reach the same pixels: they now carry its values.
  sk2 <- sk; sk2[10, 2:18] <- TRUE; sk2[10, 10:11] <- TRUE
  m2 <- m; m2$skeleton <- sk2; m2$branches <- c(trace_skeleton(m2, 10L), m$branches)
  y2 <- map_mat_properties(d, compute_mat_properties(m2, P))
  expect_true(all(is.finite(y2$mat_images$parallelism[ink & is.finite(y$mat_images$separation)])))
  # The running example's two-pixel branch inside the V is NA for parallelism,
  # and the mapped contour values come from the branch above it.
  seg <- function(...) matrix(c(...), nrow = 1)
  ex <- mlv_drawing(list(seg(10, 18, 70, 18), seg(10, 32, 70, 30), rbind(c(20, 40, 40, 44), c(40, 44, 60, 40))), 80, 50)
  a <- analyse_mat(ex, properties = P)
  short <- Filter(function(b) nrow(b) == 2L, a$mat$branches); expect_length(short, 1L)
  expect_true(all(is.na(short[[1]]$parallelism)))
  expect_false(any(a$drawing$mat_images$parallelism[rasterize(ex)] == 0, na.rm = TRUE))
})

test_that("genuine zeros remain valid zeros", {
  # Lines on rows 20 and 23: the skeleton lies one pixel from the nearer line,
  # radius 1, so separation is (1 - 1/1)^10 = 0 at every point and mapped pixel.
  poly <- function(x, y) cbind(head(x, -1), head(y, -1), x[-1], y[-1])
  d <- mlv_drawing(list(poly(c(10, 50), c(20, 20)), poly(c(10, 50), c(23, 23))), 60, 40)
  a <- analyse_mat(d)
  v <- a$drawing$mat_images$separation[rasterize(d)]
  expect_true(all(is.finite(v))); expect_true(all(v == 0))
  expect_true(all(unlist(lapply(a$mat$branches, `[[`, "separation")) == 0))
})

test_that("NA-aware direction averaging and mapping keep reversal invariance", {
  # Radii below 1 at one end: the 2-point moving mean of the radius is
  # direction dependent, so a point can be defined in one direction only.
  m <- line_mat(5:20); r <- c(.6, .8, 1.2, rep(4, 13))
  m$branches[[1]]$radius <- r; m$distance_map[10, 5:20] <- r
  A <- compute_mat_properties(m, P)
  mr <- m; mr$branches <- list(m$branches[[1]][16:1, ])
  B <- compute_mat_properties(mr, P)
  for (p in P) {
    expect_identical(is.na(A$score_images[[p]]), is.na(B$score_images[[p]]))
    expect_lt(max(abs(A$score_images[[p]] - B$score_images[[p]]), na.rm = TRUE), 1e-12)
  }
  # defined in at least one direction everywhere; the published forward
  # scores leave the first three points NA
  expect_false(anyNA(A$branches[[1]]$separation))
  expect_true(all(is.na(score_branch(m$branches[[1]], "separation")[1:3])))
})

test_that("the corrected estimator sums windows directly, so an exact radius 1 gives separation 0", {
  sample_na <- function(R, p = 10) mlvtoolbox:::.score_branch_sample(seq_along(R), rep(1, length(R)), R, "separation", 5, p, "na")
  sample_pub <- function(R) mlvtoolbox:::.score_branch_sample(seq_along(R), rep(1, length(R)), R, "separation", 5, 10, "zero")
  # radius pair (1, 1): the cumulative-sum difference gives 1 - 1.1e-16, the direct mean exactly 1
  R <- c(sqrt(2), 1, 1)
  expect_lt(mlvtoolbox:::.moving_mean(R, 2L)[3], 1)
  expect_identical(mlvtoolbox:::.moving_mean_direct(R, 2L)[3], 1)
  s2 <- 1 - 1 / ((sqrt(2) + 1) / 2)
  expect_identical(sample_na(R, 1)[3], (s2 + 0) / 2)                 # point 3 is a defined 0 before smoothing
  # where three consecutive raw values are 0, the smoothed score is exactly 0
  R5 <- c(sqrt(2), 1, 1, 1, 1)
  expect_identical(sample_na(R5)[4:5], c(0, 0))
  expect_true(all(is.na(sample_pub(R5)[2:4])))                       # the published rounding loses them
})

test_that("arclength mirror over a window whose ends coincide is NA under the default rule", {
  sq <- data.frame(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), radius = 2)
  na <- mlvtoolbox:::.score_branch(sq, "mirror", 2, 1, "arclength", 1, 1, "pixels", NULL, "na")
  zero <- mlvtoolbox:::.score_branch(sq, "mirror", 2, 1, "arclength", 1, 1, "pixels", NULL, "zero")
  expect_true(is.na(na[3])); expect_identical(zero[3], 0)            # window [0, 4] spans the closed square
  expect_identical(na[-3], zero[-3])
  expect_identical(score_branch(sq, "mirror", 2, 1, estimator = "arclength"), zero)
})

test_that("split_contours leaves contours without a finite medial-axis criterion unassigned", {
  seg <- function(...) matrix(c(...), nrow = 1)
  # Two 2 px lines 3 px apart: the only branch has two points, so parallelism is
  # undefined everywhere and neither contour can be ranked (this used to error).
  d0 <- mlv_drawing(list(seg(8, 8, 10, 8), seg(8, 11, 10, 11)), 24, 24)
  x0 <- analyse_mat(d0, properties = c("parallelism", "mirror", "taper"))$drawing
  expect_false(any(is.finite(unlist(x0$mat_properties$parallelism))))
  s0 <- split_contours(x0, "parallelism")
  expect_identical(s0$unranked, x0$ids); expect_true(all(is.na(s0$scores)))
  expect_length(s0$top$ids, 0L); expect_length(s0$bottom$ids, 0L); expect_setequal(s0$unassigned$ids, x0$ids)
  # Mixed: four contours, the first two made unrankable.
  d <- mlv_drawing(list(seg(5, 10, 40, 10), seg(5, 16, 40, 16), seg(5, 30, 40, 30), seg(5, 38, 40, 38)), 50, 50)
  x <- analyse_mat(d, properties = c("parallelism", "mirror", "taper"))$drawing
  ref <- split_contours(x, "parallelism", fraction = 1)
  expect_length(ref$unranked, 0L)
  y <- x; y$mat_properties$parallelism[1:2] <- lapply(y$mat_properties$parallelism[1:2], function(z) z * NA)
  s <- split_contours(y, "parallelism", fraction = 1)
  expect_identical(s$unranked, x$ids[1:2])
  expect_true(all(is.na(s$scores[1:2])))
  expect_identical(unname(s$scores[3:4]), unname(rank(ref$scores[3:4])))   # ranks among rankable contours
  expect_setequal(s$top$ids, x$ids[3:4]); expect_setequal(s$bottom$ids, x$ids[3:4])
  expect_setequal(s$unassigned$ids, x$ids[1:2])
  # an unrankable criterion makes the contour unranked under combined criteria too
  expect_identical(split_contours(y, c("length", "parallelism"))$unranked, x$ids[1:2])
  # the length budget is still a fraction of the total length of all four
  # equal contours: half of it holds both rankable contours (a budget over the
  # rankable ones only would hold one)
  s5 <- split_contours(y, "parallelism", fraction = .5)
  expect_setequal(s5$top$ids, x$ids[3:4])
  # the designer optimises additive contour targets only, never medial-axis criteria
  spec <- feature_spec("length", normalize = FALSE)
  expect_error(design_subsets(x, spec, target = "parallelism", length_bounds = c(0, 100)), "should be one of")
})
