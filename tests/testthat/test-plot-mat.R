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
