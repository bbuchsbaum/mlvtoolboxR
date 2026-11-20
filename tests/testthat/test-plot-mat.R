# plot.mlv_mat(): every stage renders on a raster device, invalid input fails
# with an explicit message, and no persistent graphics parameter changes.

plot_fixture <- function() {
  d <- mlv_drawing(list(matrix(c(5, 10, 45, 10), 1), matrix(c(5, 22, 45, 22), 1)), 50, 32)
  list(drawing = d, mat = compute_mat_properties(compute_mat(d), c("parallelism", "separation")))
}

with_png <- function(code) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 400, height = 300)
  on.exit({ grDevices::dev.off(); unlink(file) })
  force(code)
}

graphics_state <- function() par(c("mar", "mfrow", "mfcol", "xpd", "bg", "fg", "col", "cex", "las", "mgp", "oma"))

test_that("every stage and score renders and returns the MAT invisibly", {
  f <- plot_fixture()
  with_png({
    for (what in c("ink", "distance", "flux", "skeleton", "branches", "parallelism", "separation")) {
      expect_invisible(out <- plot(f$mat, what))
      expect_identical(out, f$mat)
      expect_no_error(plot(f$mat, what, drawing = f$drawing, legend = FALSE, axes = TRUE))
    }
    expect_no_error(plot(f$mat, "flux", zlim = c(-50, 50), col = grDevices::gray.colors(10)))
  })
})

test_that("plot coordinates are image coordinates: y down, pixel centres at integers", {
  f <- plot_fixture()
  with_png({
    plot(f$mat, "ink")
    usr <- par("usr")
    expect_gt(usr[3], usr[4])                          # reversed y axis
    expect_true(usr[3] >= 32.5 && usr[4] <= 0.5)       # the whole image is visible
    expect_lt(grconvertY(1, "user", "device"), grconvertY(32, "user", "device"))  # row 1 on top
    expect_lte(usr[1], 0.5)
    expect_gt(usr[2], 50.5)                            # legend room right of the image
    # equal aspect: one pixel spans the same device distance in x and y
    expect_equal(abs(diff(grconvertX(c(0, 10), "user", "device"))),
                 abs(diff(grconvertY(c(0, 10), "user", "device"))), tolerance = 1e-6)
  })
})

test_that("graphics state is unchanged after plotting", {
  f <- plot_fixture()
  with_png({
    par(mfrow = c(1, 2), mar = c(2, 3, 4, 1))
    before <- graphics_state()
    plot(f$mat, "flux", drawing = f$drawing)
    plot(f$mat, "skeleton")
    expect_identical(graphics_state(), before)
  })
})

test_that("invalid input fails with explicit messages", {
  f <- plot_fixture()
  with_png({
    expect_error(plot(f$mat, "taper"), "what must be one of")
    expect_error(plot(f$mat, c("ink", "flux")), "what must be one of")
    expect_error(plot(compute_mat(f$drawing), "parallelism"), "compute_mat_properties")
    expect_error(plot(f$mat, "flux", zlim = c(1, 0)), "zlim")
    expect_error(plot(f$mat, "flux", legend = NA), "legend")
    other <- mlv_drawing(list(matrix(c(1, 1, 5, 5), 1)), 20, 20)
    expect_error(plot(f$mat, "flux", drawing = other), "dimensions differ")
    expect_error(plot.mlv_mat(list(), "flux"), "mlv_mat")
  })
})

test_that("empty skeletons and empty ink render without error", {
  single <- mlv_drawing(list(matrix(c(10, 20, 50, 20), 1)), 60, 40)
  m <- compute_mat(single)
  expect_false(any(m$skeleton))
  blank <- compute_mat(matrix(FALSE, 12, 15))
  with_png({
    for (what in c("ink", "distance", "flux", "skeleton", "branches"))
      expect_no_error(plot(m, what, drawing = single))
    for (what in c("ink", "distance", "flux", "skeleton", "branches"))
      expect_no_error(plot(blank, what))
  })
})

test_that("the raster is drawn with row 1 at the top (not flipped)", {
  skip_if_not_installed("png")
  mask <- matrix(FALSE, 20, 20); mask[2:4, ] <- TRUE   # ink only near the top
  m <- compute_mat(mask)
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 200, height = 200)
  par(mar = c(0, 0, 0, 0))
  plot(m, "ink", legend = FALSE, main = "")
  ink_top <- grconvertY(3, "user", "device"); ink_bottom <- grconvertY(17, "user", "device")
  centre_x <- grconvertX(10, "user", "device")
  grDevices::dev.off(); on.exit(unlink(file))
  img <- png::readPNG(file)
  dark <- function(y) mean(img[round(y), round(centre_x), 1:3]) < 0.5
  expect_true(dark(ink_top))
  expect_false(dark(ink_bottom))
})

test_that("the branches legend has one swatch per branch, in the drawn colours", {
  key <- mlvtoolbox:::.mat_branch_key(2L)
  expect_identical(key$labels, c("ink", "branch 1", "branch 2"))
  expect_identical(key$fill[-1], mlvtoolbox:::.mat_branch_colours[1:2])
  expect_false(identical(key$fill[2], key$fill[3]))
  many <- mlvtoolbox:::.mat_branch_key(9L)
  expect_length(many$labels, 7L)
  expect_match(many$labels[7], "\\+ 4 more")
  expect_identical(mlvtoolbox:::.mat_branch_key(0L)$labels, "ink")
  f <- plot_fixture()
  with_png(expect_no_error(plot(f$mat, "branches", cex.main = 1.4, font.main = 2)))
})

test_that("every stage renders on a pdf device without encoding warnings", {
  f <- plot_fixture()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) })
  for (what in c("ink", "distance", "flux", "skeleton", "branches", "parallelism", "separation"))
    expect_no_warning(plot(f$mat, what, drawing = f$drawing))
})

test_that("default palettes contain no missing colours", {
  # hcl.colors(101, "Blue-Red 3")[51] is NA in R 4.5.1; an NA colour is drawn
  # transparent, so zero flux would read as "no data".
  for (what in c("distance", "flux", "parallelism")) {
    pal <- mlvtoolbox:::.mat_palette(what)
    expect_length(pal, 101L)
    expect_false(anyNA(pal))
  }
  expect_false(anyNA(mlvtoolbox:::.no_na_colours(grDevices::hcl.colors(101, "Blue-Red 3"))))
  # The flux midpoint is a neutral near-white, not missing.
  mid <- grDevices::col2rgb(mlvtoolbox:::.mat_palette("flux")[51])
  expect_true(all(mid > 230))
})

test_that("a zero-flux pixel is painted, and only missing values are transparent", {
  z <- matrix(c(-1, 0, 1, NA), 2)
  pal <- mlvtoolbox:::.mat_palette("flux")
  painted <- mlvtoolbox:::.mat_colour_matrix(z, pal, c(-1, 1))
  expect_identical(painted[2, 1], pal[51])
  expect_identical(painted[2, 2], "transparent")
  expect_false(anyNA(painted))
  skip_if_not_installed("png")
  # Render a MAT with a region of exactly zero flux on a magenta background:
  # the zero pixels must not show the background through.
  d <- mlv_drawing(list(matrix(c(5, 10, 45, 10), 1), matrix(c(5, 22, 45, 22), 1)), 50, 32)
  m <- compute_mat(d)
  zero <- which(m$AOF == 0, arr.ind = TRUE)
  expect_gt(nrow(zero), 0L)
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 500, height = 320, bg = "#FF00FF")
  par(mar = c(0, 0, 0, 0))
  plot(m, "flux", legend = FALSE, main = "", levels = numeric(0))
  px <- grconvertX(zero[1, 2], "user", "device"); py <- grconvertY(zero[1, 1], "user", "device")
  grDevices::dev.off(); on.exit(unlink(file))
  rgb <- png::readPNG(file)[round(py), round(px), 1:3]
  expect_true(all(rgb > 0.85))                 # near-white, not magenta
})

test_that("discs, isolines, zoom and grid render and are validated", {
  f <- plot_fixture()
  with_png({
    expect_no_error(plot(f$mat, "branches", discs = 5))
    expect_no_error(plot(f$mat, "parallelism", discs = 3, legend = FALSE))
    expect_no_error(plot(f$mat, "distance", levels = 1:8))
    expect_no_error(plot(f$mat, "flux", levels = numeric(0)))
    plot(f$mat, "ink", xlim = c(10, 20), ylim = c(15, 5), grid = TRUE, legend = FALSE)
    usr <- par("usr")
    expect_lte(usr[1], 10); expect_gte(usr[2], 20)
    expect_gt(usr[3], usr[4])                    # y still increases downward
    expect_error(plot(f$mat, "ink", discs = -1), "discs")
    expect_error(plot(f$mat, "ink", discs = c(1, 2)), "discs")
    expect_error(plot(f$mat, "distance", levels = "a"), "levels")
    expect_error(plot(f$mat, "ink", grid = NA), "grid")
    expect_error(plot(f$mat, "ink", xlim = c(1, 1)), "xlim")
    expect_error(plot(f$mat, "ink", ylim = c(1, NA)), "ylim")
  })
})

test_that("every stage reserves the same legend room, so stage figures align", {
  f <- plot_fixture()
  usr <- lapply(c("ink", "distance", "flux", "skeleton", "branches", "parallelism"),
                function(what) with_png({ par(mar = c(0, 0, 0, 0)); plot(f$mat, what); par("usr") }))
  for (u in usr[-1]) expect_equal(u, usr[[1]])
  with_png({
    par(mar = c(0, 0, 0, 0))
    plot(f$mat, "flux", legend = FALSE)
    u <- par("usr")
    expect_true(isTRUE(all.equal(u[2], 50.5)) || isTRUE(all.equal(u[3], 32.5)))
  })
})

test_that("unscored skeleton pixels get a legend entry", {
  short <- compute_mat_properties(compute_mat(mlv_drawing(list(
    matrix(c(5, 10, 45, 10), 1), matrix(c(5, 22, 45, 22), 1),
    rbind(c(15, 28, 25, 31), c(25, 31, 35, 28))), 50, 36)), "parallelism")
  undefined <- short$skeleton & !is.finite(short$score_images$parallelism)
  expect_true(any(undefined))
  skip_if_not_installed("png")
  swatch_pixels <- function(legend) {
    file <- tempfile(fileext = ".png")
    grDevices::png(file, width = 600, height = 400)
    par(mar = c(0, 0, 0, 0))
    plot(short, "parallelism", legend = legend, main = "")
    edge <- grconvertX(50.5, "user", "device")
    grDevices::dev.off(); on.exit(unlink(file))
    img <- png::readPNG(file)
    right <- img[, seq(ceiling(edge) + 2, dim(img)[2]), 1:3, drop = FALSE]
    grey75 <- grDevices::col2rgb("grey75")[, 1] / 255
    sum(abs(right[, , 1] - grey75[1]) < 0.02 & abs(right[, , 2] - grey75[2]) < 0.02 &
        abs(right[, , 3] - grey75[3]) < 0.02)
  }
  expect_gt(swatch_pixels(TRUE), 20)             # the "no score" swatch is drawn
  expect_identical(swatch_pixels(FALSE), 0L)
})

test_that("branch colours can be replaced and the key follows them", {
  key <- mlvtoolbox:::.mat_branch_key(3L, col = c("red", "blue"))
  expect_identical(key$fill, c("black", "red", "blue", "red"))
  f <- plot_fixture()
  with_png(expect_no_error(plot(f$mat, "branches", col = "grey40", discs = 4)))
})

test_that("colour_key() for a MAT matches the stage's colours, limits and threshold", {
  f <- plot_fixture()
  with_png({
    layout(matrix(1:2), heights = c(1, lcm(1.3)))
    plot(f$mat, "flux", legend = FALSE)
    key <- colour_key(f$mat, "flux")
    limit <- max(abs(f$mat$AOF))
    expect_equal(key$limits, c(-limit, limit))
    expect_no_error(colour_key(f$mat, "parallelism", horizontal = FALSE))
    expect_no_error(colour_key(f$mat, "distance", zlim = c(0, 5)))
    expect_error(colour_key(f$mat, "skeleton"), "what must be one of")
    expect_error(colour_key(f$mat, "flux", zlim = c(2, 1)), "zlim")
  })
  expect_false(anyNA(mlvtoolbox:::.mat_scale(f$mat, "flux")$colours))
  expect_equal(mlvtoolbox:::.mat_scale(f$mat, "flux")$mark,
               f$mat$parameters$threshold * f$mat$parameters$samples)
})

test_that("discs crossing the shown region are omitted and the flux stage draws ink black", {
  b <- data.frame(x = c(5, 20), y = c(10, 10), radius = c(8, 3))
  view <- list(x = c(0.5, 40.5), y = c(0.5, 20.5))
  with_png({
    plot.new(); plot.window(c(0, 40), c(0, 20))
    drawn <- mlvtoolbox:::.mat_draw_discs(list(b), 1, view)
    expect_identical(drawn[[1]], 2L)          # the disc at x = 5 (radius 8) crosses x = 0.5
    expect_identical(mlvtoolbox:::.mat_draw_discs(list(b), 1)[[1]], 1:2)
  })
  skip_if_not_installed("png")
  f <- plot_fixture()
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 500, height = 320)
  par(mar = c(0, 0, 0, 0))
  plot(f$mat, "flux", legend = FALSE, main = "", levels = numeric(0))
  ink <- which(f$mat$mask, arr.ind = TRUE)[1, ]
  px <- grconvertX(ink[2], "user", "device"); py <- grconvertY(ink[1], "user", "device")
  grDevices::dev.off(); on.exit(unlink(file))
  expect_true(all(png::readPNG(file)[round(py), round(px), 1:3] < 0.1))
})
