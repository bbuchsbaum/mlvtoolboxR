#' Draw a colour key in its own panel
#'
#' Draws a colour bar that fills the current figure region, for use in a
#' [graphics::layout()] row or column beside panels drawn with
#' `legend = FALSE`. One key can then serve several panels that share a
#' scale, and the key never takes room from the data region of a panel.
#' The key is drawn in fixed physical sizes (text at `cex`, a bar about 0.12
#' inch thick), so it reads the same in a wide figure and in a narrow phone
#' copy of it; give it a strip about 0.5 inch tall (horizontal) or 1 inch wide
#' (vertical), for example `layout(..., heights = c(1, lcm(1.3)))`.
#'
#' The default method takes the colours and limits directly. The method for
#' drawings derives them exactly as [plot.mlv_drawing()] does for the same
#' `property`, `palette` and `zlim`, so a panel drawn with `legend = FALSE`
#' and its key always agree.
#'
#' The function opens a new plot (it calls [graphics::plot.new()]) and
#' changes no persistent graphics parameter.
#' @param x A character vector of at least two colours, or a palette function
#'   of the number of colours (default method); a drawing (drawing method).
#' @param limits Two increasing finite numbers: the values of the first and
#'   last colour.
#' @param label The quantity and its unit, for example `"length (px)"`.
#' @param horizontal Draw the bar horizontally (default) or vertically.
#' @param ticks Tick positions; defaults to [pretty()] values within `limits`,
#'   or five equally spaced values for a cyclic key.
#' @param mark Optional value marked with a heavier tick across the bar, such
#'   as a threshold, and `mark_label` its label.
#' @param mark_label Label for `mark`.
#' @param cyclic Treat the scale as cyclic, as for axial orientation: the
#'   colours wrap so that both ends share the first colour.
#' @param glyphs Draw a short line at each tick at the angle given by the tick
#'   value in degrees (counterclockwise, as drawn on screen). Useful for
#'   orientation keys.
#' @param cex Text size.
#' @param ... Passed to methods.
#' @return The colours, limits and label used, invisibly.
#' @examples
#' roof <- rbind(c(5, 30, 20, 10), c(20, 10, 35, 30))
#' d <- mlv_drawing(list(roof, rbind(c(2, 35, 38, 35))), 40, 40)
#' op <- par(no.readonly = TRUE)
#' layout(matrix(1:2), heights = c(1, lcm(1.3)))
#' plot(d, property = "orientation", legend = FALSE, lwd = 2)
#' colour_key(d, "orientation")
#' par(op)
#' colour_key(grDevices::hcl.colors(64, "Viridis"), c(0, 1), "score (0-1)")
#' @export
colour_key <- function(x, ...) UseMethod("colour_key")

#' @rdname colour_key
#' @export
colour_key.default <- function(x, limits, label = NULL, horizontal = TRUE, ticks = NULL,
                               mark = NULL, mark_label = NULL, cyclic = FALSE,
                               glyphs = FALSE, cex = 0.8, ...) {
  if (is.function(x)) palette <- x else {
    if (!is.character(x) || length(x) < 2L || anyNA(x))
      stop("x must be a palette function or at least two colours", call. = FALSE)
    palette <- grDevices::colorRampPalette(x)
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[2] <= limits[1])
    stop("limits must be two increasing finite numbers", call. = FALSE)
  for (flag in list(horizontal, cyclic, glyphs))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("horizontal, cyclic and glyphs must be TRUE or FALSE", call. = FALSE)
  n <- 128L
  cols <- if (cyclic) .map_colours((seq_len(n) - 0.5) / n * diff(limits) + limits[1],
                                   limits, palette, cyclic = TRUE) else palette(n)
  cols <- .no_na_colours(cols)
  if (is.null(ticks)) ticks <- if (cyclic) seq(limits[1], limits[2], length.out = 5L) else {
    z <- pretty(limits, 4L); z[z >= limits[1] - 1e-12 & z <= limits[2] + 1e-12]}
  ticks <- ticks[is.finite(ticks) & ticks >= limits[1] - 1e-12 & ticks <= limits[2] + 1e-12]

  # The key fills the figure region whatever the margins; mar is restored.
  op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  .replayable(.draw_key, cols = cols, limits = limits, label = label, horizontal = horizontal,
              ticks = ticks, mark = mark, mark_label = mark_label, glyphs = glyphs, cex = cex)
  invisible(list(colours = cols, limits = limits, label = label))
}

#' @rdname colour_key
#' @param property For drawings: `"orientation"`, `"length"`, `"curvature"`
#'   or a mapped medial-axis score, as in [plot.mlv_drawing()].
#' @param palette,zlim As in [plot.mlv_drawing()].
#' @export
colour_key.mlv_drawing <- function(x, property, palette = NULL, zlim = NULL, ...) {
  .check_drawing(x)
  bar <- .drawing_bar(x, property, palette, zlim)
  if (is.null(bar)) stop("property is not colour coded; use orientation, length, ",
                         "curvature or a mapped medial-axis score", call. = FALSE)
  # A constant property is drawn in the middle colour (as by plot()); centre
  # it on a unit-wide key.
  limits <- if (diff(bar$range) > 0) bar$range else bar$range + c(-0.5, 0.5)
  colour_key.default(bar$palette, limits, sub("\n", " ", bar$label, fixed = TRUE),
                     cyclic = bar$cyclic, ...)
}

# Replace missing colours (some hcl.colors() palettes return NA at an odd
# midpoint on some R versions) by interpolating their neighbours, so a palette
# never paints a value as "no data".
.no_na_colours <- function(cols) {
  bad <- is.na(cols)
  if (!any(bad)) return(cols)
  if (all(bad)) stop("palette returned no colours", call. = FALSE)
  if (sum(!bad) == 1L) return(rep(cols[!bad], length(cols)))
  rgb <- grDevices::col2rgb(cols[!bad], alpha = TRUE)
  idx <- which(!bad)
  fill <- vapply(seq_len(4L), function(k)
    stats::approx(idx, rgb[k, ], xout = seq_along(cols), rule = 2)$y, numeric(length(cols)))
  if (is.null(dim(fill))) fill <- matrix(fill, nrow = 1L)
  cols[bad] <- grDevices::rgb(fill[bad, 1], fill[bad, 2], fill[bad, 3], fill[bad, 4],
                              maxColorValue = 255)
  cols
}

# The key itself, in inches from the figure region's corners so that it keeps
# its physical size; replayed with the display list (see .replayable()).
.draw_key <- function(cols, limits, label, horizontal, ticks, mark, mark_label, glyphs, cex) {
  # Work in inches from the figure region's corners, so the key keeps its
  # physical size whatever the panel size and margins.
  fx <- function(inch) graphics::grconvertX(inch, "inches", "user")
  fy <- function(inch) graphics::grconvertY(inch, "inches", "user")
  fig <- graphics::par("fig"); din <- graphics::par("din")
  left <- fig[1] * din[1]; right <- fig[2] * din[1]
  bottom <- fig[3] * din[2]; top <- fig[4] * din[2]
  line_h <- graphics::strheight("Mg", units = "inches", cex = cex)
  thick <- 0.12; tick <- 0.05
  fmt <- function(v) format(v, trim = TRUE, drop0trailing = TRUE)
  if (horizontal) {
    width <- right - left
    x0 <- left + 0.12 * width; x1 <- right - 0.12 * width
    g <- if (glyphs) 0.07 else 0
    # Stack, top to bottom: label, bar, ticks and their values, then any
    # glyphs. The stack without glyphs is centred in the strip and glyphs hang
    # below it, so keys with and without glyphs put their bars at one height.
    head <- if (is.null(label)) 0 else 1.6 * line_h
    stack <- head + thick + tick + 1.15 * line_h
    y1 <- (top + bottom) / 2 + stack / 2 - head
    y0 <- y1 - thick
    at <- function(v) x0 + (v - limits[1]) / diff(limits) * (x1 - x0)
    graphics::rasterImage(grDevices::as.raster(matrix(cols, nrow = 1L)), fx(x0), fy(y0),
                          fx(x1), fy(y1), interpolate = FALSE, xpd = NA)
    graphics::rect(fx(x0), fy(y0), fx(x1), fy(y1), border = "grey40", lwd = 0.5, xpd = NA)
    graphics::segments(fx(at(ticks)), fy(y0), fx(at(ticks)), fy(y0 - tick),
                       col = "grey40", lwd = 0.5, xpd = NA)
    graphics::text(fx(at(ticks)), fy(y0 - tick - 0.15 * line_h), fmt(ticks), adj = c(0.5, 1),
                   cex = cex, xpd = NA)
    if (glyphs) {
      gy <- y0 - tick - 1.15 * line_h - 0.06 - g
      a <- ticks * pi / 180
      graphics::segments(fx(at(ticks) - g * cos(a)), fy(gy - g * sin(a)),
                         fx(at(ticks) + g * cos(a)), fy(gy + g * sin(a)),
                         lwd = 1.2, col = "grey20", xpd = NA)
    }
    if (!is.null(mark)) {
      graphics::segments(fx(at(mark)), fy(y0 - 0.03), fx(at(mark)), fy(y1 + 0.03),
                         lwd = 1.6, xpd = NA)
      if (!is.null(mark_label))
        graphics::text(fx(at(mark)), fy(y1 + 0.06), mark_label, adj = c(0.5, 0),
                       cex = 0.9 * cex, xpd = NA)
    }
    if (!is.null(label))
      graphics::text(fx(x0), fy(y1 + 0.55 * line_h), label, adj = c(0, 0), cex = cex, xpd = NA)
  } else {
    height <- top - bottom
    x0 <- left + 0.1 * (right - left); x1 <- x0 + thick
    y0 <- bottom + 0.08 * height
    y1 <- top - 0.08 * height - (if (is.null(label)) 0 else
      (length(strsplit(label, "\n", fixed = TRUE)[[1]]) + 0.6) * line_h * 1.2)
    at <- function(v) y0 + (v - limits[1]) / diff(limits) * (y1 - y0)
    graphics::rasterImage(grDevices::as.raster(matrix(rev(cols), ncol = 1L)), fx(x0), fy(y0),
                          fx(x1), fy(y1), interpolate = FALSE, xpd = NA)
    graphics::rect(fx(x0), fy(y0), fx(x1), fy(y1), border = "grey40", lwd = 0.5, xpd = NA)
    graphics::segments(fx(x1), fy(at(ticks)), fx(x1 + tick), fy(at(ticks)),
                       col = "grey40", lwd = 0.5, xpd = NA)
    graphics::text(fx(x1 + tick + 0.04), fy(at(ticks)), fmt(ticks), adj = c(0, 0.5),
                   cex = cex, xpd = NA)
    if (!is.null(mark)) {
      graphics::segments(fx(x0 - 0.03), fy(at(mark)), fx(x1 + 0.03), fy(at(mark)),
                         lwd = 1.6, xpd = NA)
      if (!is.null(mark_label))
        graphics::text(fx(x0 - 0.05), fy(at(mark)), mark_label, adj = c(1, 0.5),
                       cex = 0.9 * cex, xpd = NA)
    }
    if (!is.null(label))
      graphics::text(fx(x0), fy(y1 + 0.5 * line_h), label, adj = c(0, 0), cex = cex, xpd = NA)
  }
  invisible(NULL)
}
