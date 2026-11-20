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

test_that("reach draws the extension along the stem, longer on the diagonal under Chebyshev", {
  # A stem ending at (20, 20) and running up-right at 45 degrees; its end
  # extends down-left along its own line. 1.2 px along that line lies beyond
  # the 1 px Euclidean reach but within the sqrt(2) px Chebyshev reach.
  d <- mlv_drawing(list(sg(20, 20, 32, 8)), 40, 40)
  draw <- function(metric) function()
    plot(d, reach = metric, axes = FALSE, frame = FALSE, lwd = 0.5)
  probe <- c(20, 20) + 1.2 * c(-1, 1) / sqrt(2)
  blue <- function(rgb) rgb[3] > rgb[1] + 0.2
  expect_true(blue(at(shoot(draw("chebyshev"), 40, 40), probe[1], probe[2])))
  expect_false(blue(at(shoot(draw("euclidean"), 40, 40), probe[1], probe[2])))
  # Under both metrics the extension covers 0.6 px along the line.
  mid <- c(20, 20) + 0.6 * c(-1, 1) / sqrt(2)
  expect_true(blue(at(shoot(draw("euclidean"), 40, 40), mid[1], mid[2])))
  # The outline is not filled: a point inside the disc but off the line stays white.
  expect_true(is_white(at(shoot(draw("euclidean"), 40, 40), 19.6, 19.4)))
  # The reach length follows detect_junctions(): metric, then relative cap.
  len <- mlvtoolbox:::.reach_length
  expect_equal(len(c(1, 1) / sqrt(2), 20, "chebyshev", c(1, 0.3)), sqrt(2))
  expect_equal(len(c(1, 0), 20, "chebyshev", c(1, 0.3)), 1)
  expect_equal(len(c(1, 1) / sqrt(2), 2, "euclidean", c(1, 0.3)), 0.6)
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

test_that("gaps can be limited to the smallest or the largest gap", {
  arrow <- mlv_drawing(list(rbind(c(8, 26, 20, 20), c(20, 20, 32, 26)), sg(20, 20, 20, 34)), 40, 40)
  with_png({
    expect_silent(plot(arrow, arms = TRUE, gaps = "min"))
    expect_silent(plot(arrow, arms = TRUE, gaps = "max"))
    expect_error(plot(arrow, arms = TRUE, gaps = "all"), "TRUE or FALSE")
  })
  # Only the largest gap's sector is shaded under "max"; "min" shades none.
  shade_at <- function(gaps) {
    at(shoot(function() plot(arrow, property = "junctions", arms = TRUE, gaps = gaps,
                             arm_length = 10, axes = FALSE, frame = FALSE), 40, 40), 20, 18)
  }
  expect_false(is_white(shade_at("max")))
  expect_true(is_white(shade_at("min")))
})

test_that("coincident arms keep their direction and stack their heads", {
  # The published rule gives this overshoot two coincident downward arms.
  d <- mlv_drawing(list(sg(5, 10, 35, 10), rbind(c(20, 9, 20, 11), c(20, 11, 20, 30))), 40, 32)
  d <- compute_contour_properties(d, rule = "mlv2023")
  shot <- shoot(function() plot(d, property = "junctions", arms = TRUE, arm_length = 8,
                                axes = FALSE, frame = FALSE), 40, 32)
  # Both heads lie on the stem's line (x = 20), one behind the other; nothing
  # is drawn at the fanned positions either side.
  orange <- function(rgb) rgb[1] > 0.6 && rgb[3] < 0.3
  expect_true(orange(at(shot, 20.25, 17.6)))
  expect_false(orange(at(shot, 21.5, 17.6)))
  expect_false(orange(at(shot, 18.5, 17.6)))
})

test_that("ink distance measures to the nearest segment", {
  ink <- rbind(c(0, 0, 10, 0), c(0, 5, 10, 5))
  expect_equal(mlvtoolbox:::.ink_distance(c(5, 1), ink), 1)
  expect_equal(mlvtoolbox:::.ink_distance(c(13, 0), ink), 3)
  expect_identical(mlvtoolbox:::.ink_distance(c(1, 1), ink[0, , drop = FALSE]), Inf)
})

test_that("legend labels use the degree sign", {
  d <- mlv_drawing(list(sg(5, 10, 35, 10)), 40, 40)
  expect_match(mlvtoolbox:::.plot_label(d, "orientation"), "°, mod 180", fixed = TRUE)
  expect_match(mlvtoolbox:::.plot_label(d, "curvature"), "°/px", fixed = TRUE)
})

test_that("arrowhead size can be set and is validated", {
  d <- mlv_drawing(list(sg(6, 20, 34, 20), sg(20, 20, 20, 34)), 40, 40)
  with_png({
    expect_silent(plot(d, property = "junctions", arms = TRUE, arm_head = 0.04))
    expect_error(plot(d, arms = TRUE, arm_head = 0), "arm_head")
  })
  big <- shoot(function() plot(d, property = "junctions", arms = TRUE, arm_length = 10,
                               arm_head = 0.3, axes = FALSE, frame = FALSE), 40, 40)
  small <- shoot(function() plot(d, property = "junctions", arms = TRUE, arm_length = 10,
                                 arm_head = 0.03, axes = FALSE, frame = FALSE), 40, 40)
  orange <- function(img) sum(img[, , 1] > 0.6 & img[, , 3] < 0.3)
  expect_gt(orange(big), orange(small))
})

test_that("a mapped property's underlay takes col, and col = NA leaves it out", {
  # Shallow diagonal lines: the vector underlay shows as slivers beside the
  # stepped raster pixels.
  d <- mlv_drawing(list(sg(4, 5, 26, 9), sg(4, 15, 26, 19)), 30, 22)
  mapped <- analyse_mat(d, properties = "parallelism", min_branch_points = 2)$drawing
  draw <- function(col) function() plot(mapped, property = "parallelism", legend = FALSE,
                                        axes = FALSE, frame = FALSE, lwd = 3, col = col)
  red <- function(img) sum(img[, , 1] > 0.8 & img[, , 2] < 0.2 & img[, , 3] < 0.2)
  expect_gt(red(shoot(draw("red"), 30, 22)), 0)
  expect_equal(red(shoot(draw(NA), 30, 22)), 0)
})
