# Junction views of plot.mlv_drawing(): arms, gaps, reach, pixel grid and
# overlays; and colour_key().
sg <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)
# Render f() to a k-pixels-per-unit PNG with no margins over a w x h canvas
# and return the image array; image pixel (x, y) is at device column
# round((x - 0.5) * k), row round((y - 0.5) * k).
shoot <- function(f, w, h, k = 10) {
  skip_if_not_installed("png")
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = w * k, height = h * k)
  op <- graphics::par(mar = c(0, 0, 0, 0))
  f()
  graphics::par(op); grDevices::dev.off()
  on.exit(unlink(file))
  png::readPNG(file)
}
at <- function(img, x, y, k = 10) img[round((y - 0.5) * k), round((x - 0.5) * k), 1:3]
is_dark <- function(rgb) all(rgb < 0.3)
is_white <- function(rgb) all(rgb > 0.95)
with_png <- function(code) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 500, height = 400)
  on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
  force(code)
}

test_that("junction views render with every new option", {
  t_shape <- mlv_drawing(list(sg(6, 20, 34, 20), sg(20, 20, 20, 34)), 40, 40)
  diagonal <- mlv_drawing(list(sg(10, 10, 30, 30), sg(20, 20, 30, 10)), 40, 40)
  mlv <- compute_contour_properties(diagonal, rule = "mlv2023")
  with_png({
    expect_invisible(plot(t_shape, property = "junctions", arms = TRUE, gaps = TRUE))
    expect_silent(plot(mlv, property = "junctions", arms = TRUE, fan = 10,
                       reference_junctions = detect_junctions(diagonal)))
    expect_silent(plot(mlv, arms = TRUE, fan = 0, gaps = TRUE))
    expect_silent(plot(t_shape, pixel_grid = TRUE, reach = "chebyshev",
                       xlim = c(15, 25), ylim = c(25, 15)))
    expect_silent(plot(t_shape, add = TRUE, col = "grey80"))
  })
})

test_that("new options are validated", {
  d <- mlv_drawing(list(sg(6, 20, 34, 20), sg(20, 20, 20, 34)), 40, 40)
  with_png({
    expect_error(plot(d, arms = TRUE, fan = -1), "fan")
    expect_error(plot(d, gaps = NA), "TRUE or FALSE")
    expect_error(plot(d, add = "yes"), "TRUE or FALSE")
    expect_error(plot(d, reach = "manhattan"), "should be one of")
    expect_error(plot(d, reach = "euclidean", reach_epsilon = 1), "reach_epsilon")
    expect_error(plot(d, arms = TRUE, reference_junctions = data.frame(a = 1)),
                 "reference_junctions")
    expect_error(plot(d, arms = TRUE, arm_lwd = 0), "arm_lwd")
  })
})

test_that("arms are drawn beneath the ink, so they never hide a stroke", {
  d <- mlv_drawing(list(sg(2, 10, 38, 10), sg(20, 10, 20, 30)), 40, 32)
  shot <- shoot(function() plot(d, property = "junctions", arms = TRUE, arm_length = 12,
                                lwd = 4, axes = FALSE, frame = FALSE), 40, 32)
  # Halfway along the right and downward arms the ink is still black on top.
  expect_true(is_dark(at(shot, 26, 10)))
  expect_true(is_dark(at(shot, 20, 16)))
  # Beside the ink the translucent arm band shows (not white, not black).
  band <- at(shot, 26, 10.25)
  expect_false(is_white(band)); expect_false(is_dark(band))
})

test_that("reach draws the metric's ball: a square reaches further on the diagonal", {
  # A stem ending at (20, 20); the point 0.9 px right and 0.9 px down lies
  # inside the Chebyshev square of half-width 1 but outside the unit disc.
  d <- mlv_drawing(list(sg(20, 20, 32, 8)), 40, 40)
  draw <- function(metric) function()
    plot(d, reach = metric, axes = FALSE, frame = FALSE, lwd = 0.5)
  corner <- c(20 + 0.9, 20 + 0.9)
  expect_false(is_white(at(shoot(draw("chebyshev"), 40, 40), corner[1], corner[2])))
  expect_true(is_white(at(shoot(draw("euclidean"), 40, 40), corner[1], corner[2])))
  # Closed contours have no free ends and get no allowance drawn.
  ring <- mlv_drawing(list(rbind(sg(10, 10, 30, 10), sg(30, 10, 20, 30), sg(20, 30, 10, 10))), 40, 40)
  plain <- shoot(function() plot(ring, axes = FALSE, frame = FALSE), 40, 40)
  reach <- shoot(function() plot(ring, reach = "chebyshev", axes = FALSE, frame = FALSE), 40, 40)
  expect_identical(plain, reach)
})

test_that("the pixel grid lies on pixel boundaries, half-way between centres", {
  d <- mlv_drawing(list(sg(1, 1, 2, 1)), 10, 10)
  shot <- shoot(function() plot(d, col = NA, pixel_grid = TRUE, axes = FALSE, frame = FALSE),
                10, 10, k = 20)
  # The boundary y = 4.5 falls on device row (4.5 - 0.5) * 20 = 80; the
  # pixel centre (5, 5) on row 90 has no line through it.
  column <- round((5 - 0.5) * 20)
  # The grid is faint, so look for any departure from white.
  expect_true(any(vapply(79:81, function(r) any(shot[r, column, 1:3] < 0.99), TRUE)))
  expect_true(is_white(at(shot, 5, 5, k = 20)))
})

test_that("add = TRUE draws in the current coordinates without a new plot", {
  d <- mlv_drawing(list(sg(5, 5, 35, 25)), 40, 30)
  with_png({
    plot(d, property = "length")
    before <- graphics::par("usr")
    plot(rotate_drawing(d, 30), add = TRUE)
    expect_identical(graphics::par("usr"), before)
  })
})

test_that("arm drawing is replayed at the new device size", {
  d <- mlv_drawing(list(sg(6, 20, 34, 20), sg(20, 20, 20, 34)), 40, 40)
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 600, height = 400)
  grDevices::dev.control("enable")
  plot(d, property = "junctions", arms = TRUE, gaps = TRUE)
  layout(matrix(1:2, 2), heights = c(1, graphics::lcm(1.3)))
  plot(d, property = "orientation", legend = FALSE)
  colour_key(d, "orientation", glyphs = TRUE)
  recorded <- grDevices::recordPlot()
  grDevices::dev.off(); unlink(file)
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 250, height = 300)
  expect_silent(grDevices::replayPlot(recorded))
  grDevices::dev.off(); unlink(file)
})

test_that("colour_key draws keys and agrees with the plotted colours", {
  d <- mlv_drawing(list(sg(5, 30, 20, 10), sg(20, 10, 35, 30)), 40, 40)
  d <- compute_contour_properties(d)
  with_png({
    k <- colour_key(d, "orientation")
    expect_identical(k$limits, c(0, 180))
    expect_false(anyNA(k$colours))
    # The cyclic key wraps: its two ends are as close as neighbours elsewhere.
    rgb <- grDevices::col2rgb(k$colours)
    step <- sqrt(colSums((rgb[, -1] - rgb[, -ncol(rgb)])^2))
    wrap <- sqrt(sum((rgb[, 1] - rgb[, ncol(rgb)])^2))
    expect_lt(wrap, 3 * stats::median(step))
    expect_silent(colour_key(d, "curvature", horizontal = FALSE))
    expect_silent(colour_key(grDevices::hcl.colors(64, "Viridis"), c(0, 1), "score",
                             mark = 0.3, mark_label = "threshold"))
    expect_silent(colour_key(function(n) grDevices::gray.colors(n), c(-1, 1), NULL,
                             horizontal = FALSE, mark = 0))
    k <- colour_key(d, "length", zlim = c(0, 100))
    expect_identical(k$limits, c(0, 100))
    expect_error(colour_key(d, "junctions"), "not colour coded")
    expect_error(colour_key("red", c(0, 1)), "at least two colours")
    expect_error(colour_key(c("red", "blue"), c(1, 0)), "increasing")
    expect_error(colour_key(c("red", "blue"), c(0, 1), cyclic = NA), "TRUE or FALSE")
  })
})

test_that("palettes with a missing colour are repaired, never painted as no data", {
  cols <- grDevices::hcl.colors(101, "Blue-Red 3")
  fixed <- mlvtoolbox:::.no_na_colours(cols)
  expect_false(anyNA(fixed))
  expect_identical(fixed[!is.na(cols)], cols[!is.na(cols)])
  # A repaired colour lies between its neighbours.
  if (anyNA(cols)) {
    i <- which(is.na(cols))[1]
    rgb <- grDevices::col2rgb(fixed[(i - 1):(i + 1)])
    expect_true(all(rgb[, 2] >= pmin(rgb[, 1], rgb[, 3]) - 1 &
                    rgb[, 2] <= pmax(rgb[, 1], rgb[, 3]) + 1))
  }
  expect_false(anyNA(mlvtoolbox:::.no_na_colours(c(NA, "#000000", NA))))
  expect_error(mlvtoolbox:::.no_na_colours(c(NA, NA)), "no colours")
})
