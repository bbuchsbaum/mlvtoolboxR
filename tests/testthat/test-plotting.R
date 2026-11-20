# Plotting: rendering, palettes, legends, junction arms and graphics state.
sg <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)
plot_drawing <- function() {
  th <- seq(0, 2 * pi, length.out = 13)
  ring <- cbind(30 + 8 * cos(th[-13]), 25 - 8 * sin(th[-13]),
                30 + 8 * cos(th[-1]), 25 - 8 * sin(th[-1]))
  mlv_drawing(list(sg(5, 10, 55, 10), sg(20, 10, 20, 45), ring), 60, 50)
}
with_png <- function(code) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 400, height = 400)
  on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
  force(code)
}

test_that("every plot property renders and returns the drawing invisibly", {
  d <- compute_contour_properties(plot_drawing())
  with_png({
    for (p in list(NULL, "orientation", "length", "curvature", "junctions")) {
      expect_invisible(out <- plot(d, property = p))
      expect_identical(out, d)
    }
    expect_silent(plot(d, property = "junctions", arms = TRUE, axes = FALSE, main = "arms"))
    expect_silent(plot(d, property = "orientation", arms = TRUE, legend = FALSE))
    expect_silent(plot(mlv_drawing(list(), 10, 10), property = "length"))
  })
})

test_that("the orientation palette is cyclic over 0 to 180 degrees", {
  pal <- mlvtoolbox:::.orientation_colours
  map <- function(v) mlvtoolbox:::.map_colours(v, c(0, 180), pal, cyclic = TRUE)
  expect_identical(map(0), map(180))
  expect_identical(map(c(10, 45)), map(c(190, 225) %% 180))
  expect_false(identical(map(0), map(90)))
  # Adjacent colours across the wrap are as close as adjacent colours elsewhere.
  rgbdist <- function(a, b) sqrt(sum((grDevices::col2rgb(a) - grDevices::col2rgb(b))^2))
  cols <- pal(256)
  steps <- vapply(seq_len(256), function(i) rgbdist(cols[i], cols[i %% 256 + 1]), 0)
  expect_lt(steps[256], 3 * stats::median(steps))
  # Directed orientations 10 and 190 degrees are the same axial orientation.
  d <- mlv_drawing(list(sg(0, 10, 10, 10 - 10 * tan(pi / 18)), sg(20, 20, 10, 20 + 10 * tan(pi / 18))), 30, 30)
  ori <- as.data.frame(d)$orientation
  expect_equal(ori, c(10, 190), tolerance = 1e-12)
  expect_identical(map(ori[1] %% 180), map(ori[2] %% 180))
})

test_that("sequential colours clamp to zlim and palettes can be supplied", {
  map <- function(v, ...) mlvtoolbox:::.map_colours(v, c(0, 10), function(n) grDevices::gray.colors(n), ...)
  expect_identical(map(-5), map(0))
  expect_identical(map(50), map(10))
  expect_true(is.na(map(NA_real_)))
  d <- compute_contour_properties(plot_drawing())
  with_png({
    expect_silent(plot(d, property = "length", zlim = c(0, 100), palette = c("white", "black")))
    expect_silent(plot(d, property = "curvature", palette = grDevices::heat.colors))
    expect_error(plot(d, property = "length", zlim = c(5, 1)), "zlim")
    expect_error(plot(d, property = "length", palette = "red"), "palette")
    expect_error(plot(d, legend = NA), "TRUE or FALSE")
    expect_error(plot(d, arms = TRUE, arm_length = -1), "arm_length")
    expect_error(plot(d, property = "absent"), "unknown")
  })
})

test_that("the colour bar takes room only when requested", {
  d <- plot_drawing()
  usr_right <- function(...) with_png({plot(d, ...); graphics::par("usr")[2]})
  plain <- usr_right()
  with_bar <- usr_right(property = "length")
  without_bar <- usr_right(property = "length", legend = FALSE)
  expect_equal(without_bar, plain)
  expect_gt(with_bar, plain + 0.15 * d$width)
  # A user xlim is respected.
  expect_lt(usr_right(property = "length", xlim = c(0, 60)), 64)
})

test_that("plotting restores graphics parameters", {
  d <- compute_contour_properties(plot_drawing())
  keep <- c("mar", "mfrow", "xpd", "col", "lwd", "cex", "lty", "las", "bg", "fg")
  with_png({
    graphics::par(mfrow = c(1, 2), mar = c(2, 3, 1, 1), las = 0)
    before <- graphics::par(keep)
    plot(d, property = "orientation", arms = TRUE)
    plot(d, property = "junctions", arms = TRUE, axes = FALSE)
    expect_identical(graphics::par(keep), before)
    expect_silent(plot_properties(d))
    expect_identical(graphics::par(keep), before)
  })
})

test_that("drawings and MAT stages share plot coordinates on the same canvas", {
  d <- plot_drawing()
  mat <- compute_mat(d)
  usr <- function(f) with_png({f(); graphics::par("usr")})
  for (shape in list(c(400, 400), c(600, 300))) {
    with_shape <- function(f) {
      file <- tempfile(fileext = ".png")
      grDevices::png(file, width = shape[1], height = shape[2])
      on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
      f(); graphics::par("usr")
    }
    expect_equal(with_shape(function() plot(d)),
                 with_shape(function() plot(mat, "ink", legend = FALSE)))
    expect_equal(with_shape(function() plot(d, property = "length")),
                 with_shape(function() plot(mat, "flux")))
  }
  # The canvas frame is the raster extent: pixel centres 1..width, edges at 0.5.
  # asp = 1 widens one axis; the other is exactly the raster extent.
  u <- usr(function() plot(d, legend = FALSE))
  x_tight <- isTRUE(all.equal(u[1:2], c(0.5, d$width + 0.5)))
  y_tight <- isTRUE(all.equal(u[3:4], c(d$height + 0.5, 0.5)))
  expect_true(x_tight || y_tight)
  expect_true(u[1] <= 0.5 && u[2] >= d$width + 0.5 && u[3] >= d$height + 0.5 && u[4] <= 0.5)
})

test_that("mapped MAT scores fill whole pixels on the raster extent", {
  skip_if_not_installed("png")
  d <- mlv_drawing(list(sg(4, 6, 26, 6), sg(4, 14, 26, 14)), 30, 20)
  mapped <- analyse_mat(d, properties = "parallelism", min_branch_points = 2)$drawing
  image <- mapped$mat_images$parallelism
  where <- which(is.finite(image))
  expect_gt(length(where), 10L)
  k <- 10L # device pixels per image pixel
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 30 * k, height = 20 * k)
  op <- graphics::par(mar = c(0, 0, 0, 0))
  plot(mapped, property = "parallelism", legend = FALSE, axes = FALSE, frame = FALSE)
  graphics::par(op); grDevices::dev.off()
  shot <- png::readPNG(file); unlink(file)
  rc <- arrayInd(where, dim(image))
  expected <- grDevices::col2rgb(mlvtoolbox:::.map_colours(image[where], c(0, 1),
    function(n) grDevices::hcl.colors(n, "Viridis"))) / 255
  # Sample each mapped pixel near its centre and near a corner of its square:
  # a filled unit square covers both, a point marker would not reach the corner.
  for (offset in list(c(0, 0), c(-0.35, -0.35))) {
    # Image pixel r spans device rows (r - 1) * k to r * k; its centre is at
    # (r - 0.5) * k.
    yy <- pmin(20 * k, pmax(1, round((rc[, 1] + offset[1] - 0.5) * k)))
    xx <- pmin(30 * k, pmax(1, round((rc[, 2] + offset[2] - 0.5) * k)))
    got <- t(sapply(seq_along(yy), function(i) shot[yy[i], xx[i], 1:3]))
    # Bitmap devices may colour-manage output (macOS shifts channels by up
    # to ~0.07); white background or a grey contour differs by more than 0.8.
    expect_true(all(rowSums(abs(got - t(expected))) < 0.25))
  }
})

test_that("range labels use an en dash only where the device can draw it", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  expect_identical(mlvtoolbox:::.range_dash(), "-")
  grDevices::dev.off(); unlink(file)
  with_png(expect_identical(mlvtoolbox:::.range_dash(), "–"))
})
