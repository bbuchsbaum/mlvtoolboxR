test_that("MAT mapping and density use the raster width compute_mat thresholded", {
  d <- mlv_drawing(list(rbind(c(5, 5, 35, 5)), rbind(c(5, 15, 35, 15))), 40, 25)
  a <- analyse_mat(d, spec = feature_spec(properties = "parallelism", raster_line_width = 3))
  wide <- rasterize(d, line_width = 3)
  expect_identical(a$mat$mask, wide)
  mapped <- a$drawing$mat_mapping_masks$parallelism
  expect_false(any(mapped & !wide))
  expect_gt(sum(mapped & !rasterize(d)), 0)
  expect_equal(a$drawing$mat_parameters$raster_line_width, 3)
  pieces <- a$drawing$mat_properties$parallelism
  expect_equal(lengths(pieces), vapply(seq_along(d$ids), function(i)
    length(mlvtoolbox:::.raster_indices(d, i, line_width = 3)), integer(1)))
  density <- feature_density(a$drawing, "parallelism")
  expect_equal(density[mapped], a$drawing$mat_images$parallelism[mapped])
})

test_that("MAT mapping rejects a specification whose raster width differs from the MAT's", {
  d <- mlv_drawing(list(rbind(c(5, 5, 35, 5)), rbind(c(5, 15, 35, 15))), 40, 25)
  wide <- feature_spec(properties = "parallelism", raster_line_width = 3)
  mat <- compute_mat_properties(compute_mat(d), spec = wide)
  expect_equal(mat$parameters$raster_line_width, 1)
  expect_error(map_mat_properties(d, mat, spec = wide), "raster width mismatch")
  # without a specification the recorded width still applies
  m3 <- compute_mat_properties(compute_mat(d, spec = wide), spec = wide)
  x <- map_mat_properties(d, m3, spec = NULL)
  expect_false(any(x$mat_mapping_masks$parallelism & !m3$mask))
  expect_equal(x$mat_parameters$raster_line_width, 3)
})

test_that("segment_at_junctions accepts double member indices", {
  d <- mlv_drawing(list(rbind(c(1, 5, 10, 5)), rbind(c(5, 1, 5, 10))), 12, 12)
  j <- detect_junctions(d)
  j$members[[1]]$contour_index <- as.double(j$members[[1]]$contour_index)
  j$members[[1]]$segment_index <- as.double(j$members[[1]]$segment_index)
  expect_equal(length(segment_at_junctions(d, j)$contours), length(segment_at_junctions(d)$contours))
})

test_that("split_by_histogram works on drawings carrying a specification", {
  d <- mlv_drawing(list(rbind(c(2, 2, 20, 2)), rbind(c(2, 8, 10, 8), c(10, 8, 10, 18))), 30, 30)
  e <- compute_contour_properties(d, spec = feature_spec())
  plain <- split_by_histogram(d, c(orientation_1 = 1))
  expect_identical(split_by_histogram(e, c(orientation_1 = 1))$scores, plain$scores)
  expect_type(split_by_histogram(e, c(orientation_1 = 1), spec = e$feature_spec)$scores, "double")
  expect_error(split_by_histogram(e, c(orientation_1 = 1), bins = 4, spec = e$feature_spec), "not both")
  expect_error(split_by_histogram(e, c(contour_id = 1)), "distinct feature-table columns")
  expect_error(split_by_histogram(e, c(orientation_1 = 1, orientation_1 = 2)), "distinct")
  expect_error(split_by_histogram(e, c(1)), "distinct")
})

test_that("zero-length cleaning reaches a fixed point and rejects underflowing lengths", {
  loop <- rbind(c(1, 1, 1 + 1e-14, 1), c(1 + 1e-14, 1, 1 + 2e-14, 1), c(1 + 2e-14, 1, 1, 1))
  d <- mlv_drawing(list(loop, rbind(c(1, 1, 5, 5))), 10, 10, zero_length = "drop")
  expect_equal(d$ids, "2")
  expect_error(mlv_drawing(list(rbind(c(0, 0, 1e-200, 0))), 10, 10), "zero-length")
  expect_length(mlv_drawing(list(rbind(c(0, 0, 1e-200, 0))), 10, 10, zero_length = "drop")$contours, 0)
  expect_no_error(contour_statistics(d))
  # small but representable geometry is still kept
  expect_length(mlv_drawing(list(rbind(c(0, 0, 1e-100, 0))), 10, 10)$contours, 1)
})

test_that("piece drawings keep parent ids aligned when a piece is dropped", {
  d <- mlv_drawing(list(rbind(c(1, 5, 10, 5)), rbind(c(5, 1, 5, 10))), 12, 12)
  j <- detect_junctions(d)
  j$x[1] <- 1 + 9e-15
  s <- segment_at_junctions(d, j)
  expect_equal(length(s$parent_ids), length(s$ids))
  expect_equal(s$parent_ids, c("1", "2", "2"))
})

test_that("feature_density smooths one-pixel-wide images along their long axis", {
  tall <- mlv_drawing(list(rbind(c(1, 2, 1, 8))), 1, 10)
  wide <- mlv_drawing(list(rbind(c(2, 1, 8, 1))), 10, 1)
  a <- feature_density(tall, sigma = 1)
  b <- feature_density(wide, sigma = 1)
  expect_equal(dim(a), c(10L, 1L))
  expect_equal(dim(b), c(1L, 10L))
  expect_equal(as.vector(a), as.vector(b))
  expect_gt(a[1, 1], 0)
})

test_that("colour mapping clamps extreme values to the end colours", {
  pal <- function(n) paste0("c", seq_len(n))
  expect_warning(out <- mlvtoolbox:::.map_colours(c(-1e12, .5, 1e12), c(0, 1), pal), NA)
  expect_equal(out, c("c1", "c129", "c256"))
})

test_that("count arguments are bounded and validated on every path", {
  img <- matrix(0, 20, 20); img[10, 3:17] <- 1
  expect_error(trace_contours(img, min_points = 1e12), "min_points")
  expect_error(compute_mat(matrix(c(TRUE, FALSE), 4, 4), min_branch_points = -5, trace_branches = FALSE),
               "min_branch_points")
  m <- compute_mat(matrix(c(TRUE, FALSE), 4, 4), trace_branches = FALSE)
  expect_error(trace_skeleton(m, min_points = 1e12), "min_points")
})

test_that("native entry points guard their inputs", {
  far <- mlvtoolbox:::cpp_raster_indices(matrix(c(-1e15, 1, 1e15, 1), 1), 5L, 5L, 3)
  expect_equal(sort(far), c(1L, 2L, 6L, 7L, 11L, 12L, 16L, 17L, 21L, 22L))
  expect_length(mlvtoolbox:::cpp_raster_indices(matrix(c(1e15, 1e15, 2e15, 2e15), 1), 5L, 5L, 3), 0)
  expect_error(mlvtoolbox:::cpp_kuwahara(numeric(10), 3L, 3L, 1L, 1, 1, 8, .01), "buffer")
  expect_error(mlvtoolbox:::cpp_line_flow(numeric(0), 0L, 3L, 1L, 1), "positive")
})

test_that("argument contracts reject ambiguous or malformed input", {
  d <- mlv_drawing(list(rbind(c(2, 2, 20, 2)), rbind(c(2, 8, 10, 8))), 30, 30)
  expect_error(render_drawing(d, "length", values = c(1, 2)), "not both")
  expect_error(batch_statistics(list(d), image_ids = ""), "image_ids")
  expect_error(as_mlv_drawing(list(contours = list(rbind(c(1, 1, 2, 2))), imsize = c(5, 5), numContours = NA)),
               "numContours")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_error(colour_key(c("red", "blue"), c(0, 1), label = character(0), horizontal = FALSE), "single string")
})

test_that("SVG import follows the specification for edge cases", {
  skip_if_not_installed("xml2")
  svg <- function(body) {
    p <- tempfile(fileext = ".svg")
    writeLines(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="40" height="40" viewBox="0 0 40 40">%s</svg>', body), p)
    p
  }
  expect_length(read_mlv_svg(svg('<rect x="1" y="1" width="0" height="10"/>'))$contours, 0)
  expect_length(read_mlv_svg(svg('<ellipse cx="5" cy="5" rx="0" ry="10"/>'))$contours, 0)
  expect_length(read_mlv_svg(svg('<circle cx="5" cy="5" r="0"/>'))$contours, 0)
  compact <- read_mlv_svg(svg('<path d="M 10 10 A 10 10 0 0120 20"/>'))
  spaced <- read_mlv_svg(svg('<path d="M 10 10 A 10 10 0 0 1 20 20"/>'))
  expect_equal(compact$contours, spaced$contours)
  expect_error(read_mlv_svg(svg('<path d="M 10 10 L 20 # 20"/>')), "invalid character")
  expect_error(read_mlv_svg(svg('<line x1="px" y1="1" x2="5" y2="5"/>')), "SVG length")
  expect_error(read_mlv_svg(svg('<line x1="  " y1="1" x2="5" y2="5"/>')), "SVG length")
  expect_length(read_mlv_svg(svg('<line x1="2px" y1=" 1 " x2="5" y2="5"/>'))$contours, 1)
})

test_that("MAT axis drawing survives a 180-degree reversal", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  graphics::plot.new(); graphics::plot.window(c(0, 5), c(0, 5))
  b <- data.frame(x = c(1, 2, 1), y = c(1, 1, 1))
  expect_no_error(mlvtoolbox:::.mat_draw_axis(list(b)))
})

test_that("design_subsets rejects invalid MAT diagnostics before searching", {
  z <- designer_v1_fixtures()$small_set; a <- z$arguments
  a$spec <- feature_spec(c("length", "curvature", "orientation"), normalize = FALSE,
    units = "visual_degrees", pixels_per_degree = 1, estimator = "arclength", curvature_scale = .1)
  solved <- FALSE
  local_mocked_bindings(.subset_solve = function(...) { solved <<- TRUE; stop("searched") })
  expect_error(do.call(design_subsets, c(list(drawings = z$drawings, diagnostics = "parallelism"), a)),
               "branch_estimator")
  expect_false(solved)
})
