#' Plot a stage of a medial-axis representation
#'
#' Draws one stage of the pipeline stored in an `mlv_mat`, in image
#' coordinates: x increases to the right, y increases downward, pixel centres
#' lie at integers and the aspect ratio is one. Quantities are colour-coded
#' with a legend that names the quantity and its unit:
#'
#' * `"ink"`: the rasterized ink mask (black on white).
#' * `"distance"`: the signed distance map in pixels (sequential palette).
#'   Ink pixels have nonpositive values.
#' * `"flux"`: the average outward flux, stored as the sum of `samples`
#'   dot products (diverging palette with limits symmetric about zero). The
#'   legend marks the skeleton threshold, `threshold * samples`.
#' * `"skeleton"`: the pixels whose flux exceeds the threshold (light) and the
#'   thinned skeleton (dark), with the ink in grey.
#' * `"branches"`: traced branches, one colour per branch, each labelled with
#'   its index in `x$branches`.
#' * a branch score name present in `x$score_images` (for example
#'   `"parallelism"`; see [compute_mat_properties()]): the score at each
#'   skeleton pixel on a 0-1 sequential palette. Pixels without a score are
#'   not drawn.
#'
#' The function changes no persistent graphics parameter: it reserves room for
#' the legend inside the plotting region, so later calls such as
#' [graphics::points()] use the same image coordinates.
#'
#' @param x An `mlv_mat` from [compute_mat()].
#' @param what Stage to draw; see Description.
#' @param drawing Optional `mlv_drawing` whose segments are drawn on top, for
#'   example the drawing the MAT was computed from. Its dimensions must match.
#' @param legend Draw a colour bar or key to the right of the image.
#' @param zlim Optional colour limits for `"distance"`, `"flux"` and scores.
#'   The default flux limits are symmetric about zero.
#' @param col Optional vector of colours replacing the default palette for a
#'   continuous stage.
#' @param main Panel title; the default names the stage. Titles use
#'   `cex.main = 1` and `font.main = 1` unless set in `...`.
#' @param axes Draw light pixel axes.
#' @param overlay_col,overlay_lwd Colour and width of the `drawing` overlay.
#' @param ... Further graphical parameters passed to [graphics::plot.window()]
#'   and [graphics::title()], such as `cex.main`.
#' @return `x`, invisibly.
#' @examples
#' d <- mlv_drawing(list(matrix(c(5, 10, 45, 10), 1),
#'                       matrix(c(5, 22, 45, 22), 1)), 50, 32)
#' m <- compute_mat(d)
#' plot(m, "flux", drawing = d)
#' plot(m, "branches")
#' @export
plot.mlv_mat <- function(x, what = "flux", drawing = NULL, legend = TRUE,
                         zlim = NULL, col = NULL, main = NULL, axes = FALSE,
                         overlay_col = "black", overlay_lwd = 0.6, ...) {
  if (!inherits(x, "mlv_mat")) stop("x must be an mlv_mat", call. = FALSE)
  stages <- c("ink", "distance", "flux", "skeleton", "branches")
  scores <- names(x$score_images)
  if (!is.character(what) || length(what) != 1L || is.na(what) ||
      !what %in% c(stages, scores))
    stop("what must be one of ", paste(c(stages, scores), collapse = ", "),
         if (!length(scores)) " (branch scores need compute_mat_properties())",
         call. = FALSE)
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend))
    stop("legend must be TRUE or FALSE", call. = FALSE)
  if (!is.null(zlim) && (!is.numeric(zlim) || length(zlim) != 2L ||
      any(!is.finite(zlim)) || zlim[2] <= zlim[1]))
    stop("zlim must be two increasing finite numbers", call. = FALSE)
  h <- nrow(x$mask); w <- ncol(x$mask)
  if (!is.null(drawing)) {
    .check_drawing(drawing)
    if (!identical(c(as.integer(drawing$height), as.integer(drawing$width)), c(h, w)))
      stop("drawing and MAT image dimensions differ", call. = FALSE)
  }
  if (is.null(main)) main <- switch(what, ink = "Ink", distance = "Signed distance map",
    flux = "Average outward flux", skeleton = "Thresholded flux and skeleton",
    branches = "Traced branches", paste("Branch score:", what))

  # Room for the legend is reserved inside the plotting region (no par() change).
  pad <- if (legend) 0.3 * max(w, h) else 0
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, w + 0.5 + pad), ylim = c(h + 0.5, 0.5),
                        asp = 1, xaxs = "i", yaxs = "i", ...)
  # rasterImage() places the first matrix row at ytop, which is the image top
  # because the y axis is reversed; image(useRaster = TRUE) would flip it.
  draw <- function(colours) graphics::rasterImage(grDevices::as.raster(colours),
    0.5, h + 0.5, w + 0.5, 0.5, interpolate = FALSE)
  paint <- function(z, colours, limits) {
    k <- 1L + floor((pmin(pmax(z, limits[1]), limits[2]) - limits[1]) /
                      diff(limits) * (length(colours) - 1L) + 0.5)
    draw(matrix(ifelse(is.na(k), "transparent", colours[k]), h, w))
  }
  layer <- function(m, colour) if (any(m))
    draw(matrix(ifelse(m, colour, "transparent"), h, w))
  range_of <- function(z) {
    if (!any(is.finite(z))) return(NULL)
    r <- range(z[is.finite(z)])
    if (r[1] == r[2]) r <- r + c(-0.5, 0.5)
    r
  }
  bar <- NULL; key <- NULL; note <- NULL
  medial <- x$skeleton & is.finite(x$distance_map) & x$distance_map > 0
  if (what == "ink") {
    layer(x$mask, "black")
  } else if (what %in% c("distance", "flux") || what %in% scores) {
    z <- switch(what, distance = x$distance_map, flux = x$AOF, x$score_images[[what]])
    limits <- if (!is.null(zlim)) zlim else if (what %in% scores) c(0, 1) else range_of(z)
    if (what == "flux" && is.null(zlim) && !is.null(limits)) limits <- c(-1, 1) * max(abs(limits))
    colours <- if (!is.null(col)) col else switch(what,
      distance = grDevices::hcl.colors(101, "Blues 3", rev = TRUE),
      flux = grDevices::hcl.colors(101, "Blue-Red 3"),
      grDevices::hcl.colors(101, "Viridis"))
    if (what %in% scores) {
      layer(x$mask, "grey80"); layer(x$skeleton & !is.finite(z), "grey60")
    }
    if (is.null(limits)) note <- "no finite values" else {
      z[!is.finite(z)] <- NA
      paint(z, colours, limits)
      bar <- list(colours = colours, limits = limits, label = switch(what,
        distance = "signed distance (px)",
        flux = paste0("flux (sum of ", x$parameters$samples, " samples)"),
        paste0(what, " (0", .range_dash(), "1)")))
      if (what == "flux" && !is.null(x$parameters$threshold))
        bar$mark <- x$parameters$threshold * x$parameters$samples
    }
  } else if (what == "skeleton") {
    above <- x$AOF > x$parameters$threshold * x$parameters$samples
    layer(above & !x$mask, "#9ecae1"); layer(x$mask, "grey55"); layer(x$skeleton, "#08306b")
    key <- list(labels = c("ink", "flux > threshold", "skeleton"),
                fill = c("grey55", "#9ecae1", "#08306b"))
  } else {
    layer(x$mask, "grey80")
    layer(x$skeleton & !medial, "grey50")
    okabe <- .mat_branch_colours
    n <- length(x$branches)
    if (!n) note <- "no branches" else for (i in seq_len(n)) {
      b <- x$branches[[i]]
      graphics::rect(b$x - 0.5, b$y - 0.5, b$x + 0.5, b$y + 0.5,
                     col = okabe[(i - 1L) %% length(okabe) + 1L], border = NA)
    }
    if (n) for (i in seq_len(n)) {
      b <- x$branches[[i]]; k <- ceiling(nrow(b) / 2)
      graphics::text(b$x[k], b$y[k], i, cex = 0.7, font = 2, pos = 3, offset = 0.4)
    }
    key <- .mat_branch_key(n)
  }
  graphics::rect(0.5, h + 0.5, w + 0.5, 0.5, border = "grey70", lwd = 0.8)
  if (!is.null(drawing) && length(drawing$contours)) {
    s <- do.call(rbind, drawing$contours)
    graphics::segments(s[, 1], s[, 2], s[, 3], s[, 4], col = overlay_col, lwd = overlay_lwd)
  }
  if (!is.null(note)) graphics::text((w + 1) / 2, (h + 1) / 2, note, col = "grey30")
  if (axes) {
    graphics::axis(1, col = "grey70", col.axis = "grey40", cex.axis = 0.7)
    graphics::axis(2, col = "grey70", col.axis = "grey40", cex.axis = 0.7, las = 1)
  }
  do.call(graphics::title, utils::modifyList(
    list(main = main, cex.main = 1, font.main = 1), list(...)))
  if (legend && !is.null(bar)) .mat_colour_bar(bar, w, h)
  if (legend && !is.null(key)) graphics::legend(w + 0.5 + 0.02 * max(w, h), 0.5,
    legend = key$labels, fill = key$fill, border = NA, bty = "n", cex = 0.7,
    xpd = NA, y.intersp = 1.1)
  invisible(x)
}

# Okabe-Ito colours (without black), recycled over branches in index order.
.mat_branch_colours <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                         "#CC79A7", "#F0E442")

# Legend for the branches stage: one swatch per branch in its own colour, up
# to six; beyond that the first five and a count of the rest.
.mat_branch_key <- function(n) {
  if (!n) return(list(labels = "ink", fill = "grey80"))
  shown <- if (n <= 6L) n else 5L
  fill <- .mat_branch_colours[(seq_len(shown) - 1L) %% length(.mat_branch_colours) + 1L]
  labels <- paste("branch", seq_len(shown))
  if (shown < n) {
    labels <- c(labels, sprintf("+ %d more (colours repeat)", n - shown))
    fill <- c(fill, NA)
  }
  list(labels = c("ink", labels), fill = c("grey80", fill))
}

# A vertical colour bar drawn in user coordinates to the right of the image,
# with tick labels and a rotated title naming the quantity and unit.
.mat_colour_bar <- function(bar, w, h) {
  s <- max(w, h)
  x0 <- w + 0.5 + 0.04 * s; x1 <- x0 + 0.035 * s
  y_bottom <- h + 0.5 - 0.06 * h; y_top <- 0.5 + 0.06 * h
  n <- length(bar$colours)
  edges <- seq(y_bottom, y_top, length.out = n + 1L)
  graphics::rect(x0, edges[-(n + 1L)], x1, edges[-1L], col = bar$colours, border = NA)
  graphics::rect(x0, y_bottom, x1, y_top, border = "grey40", lwd = 0.6)
  pos <- function(v) y_bottom + (v - bar$limits[1]) / diff(bar$limits) * (y_top - y_bottom)
  at <- pretty(bar$limits, 4); at <- at[at >= bar$limits[1] & at <= bar$limits[2]]
  labels <- format(at, trim = TRUE)
  if (!is.null(bar$mark) && bar$mark >= bar$limits[1] && bar$mark <= bar$limits[2]) {
    keep <- abs(at - bar$mark) > 0.09 * diff(bar$limits)
    at <- c(at[keep], bar$mark)
    labels <- c(labels[keep], sprintf("%.1f threshold", bar$mark))
    graphics::segments(x0, pos(bar$mark), x1, pos(bar$mark), lwd = 1.2)
  }
  graphics::segments(x1, pos(at), x1 + 0.01 * s, pos(at), col = "grey40", lwd = 0.6)
  graphics::text(x1 + 0.015 * s, pos(at), labels, adj = c(0, 0.5), cex = 0.65, xpd = NA)
  # With asp = 1, |strheight()| (user y units; negative on the reversed y axis)
  # is also the horizontal extent of the rotated title in user x units.
  x_title <- x1 + 0.015 * s + max(graphics::strwidth(labels, cex = 0.65)) +
    0.9 * abs(graphics::strheight("M", cex = 0.7))
  graphics::text(x_title, (y_bottom + y_top) / 2, bar$label, srt = 90,
                 adj = c(0.5, 0.5), cex = 0.7, xpd = NA)
}
