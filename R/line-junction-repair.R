# Find short, supported extensions from free endpoints to a crossing contour.
# Geometry proposes candidates on a simplified drawing. The original image
# supplies the evidence, and the original skeleton supplies the target pixels.
.repair_line_junctions <- function(mask, lightness, max_gap) {
  if (max_gap <= 1.5 || !any(mask)) return(mask)
  skeleton <- skeletonize(mask)
  drawing <- trace_contours(skeleton, thin = FALSE)
  if (length(drawing$contours) < 2L) return(mask)
  fitted <- simplify_contours(drawing, tolerance = 1)
  sizes <- vapply(fitted$contours, nrow, integer(1))
  segments <- do.call(rbind, fitted$contours)
  contour <- rep.int(seq_along(sizes), sizes)
  hits <- cpp_intersections(segments, contour, max_gap, 1, TRUE, TRUE)
  if (!nrow(hits)) return(mask)

  last <- cumsum(sizes)
  first <- c(1L, last[-length(last)] + 1L)
  a <- segments[, 1:2, drop = FALSE]
  b <- segments[, 3:4, drop = FALSE]
  tips <- rbind(a[first, , drop = FALSE], b[last, , drop = FALSE])
  tip_pixels <- tips[, 2L] + (tips[, 1L] - 1) * nrow(mask)
  free <- !(duplicated(tip_pixels) | duplicated(tip_pixels, fromLast = TRUE))
  closed <- tip_pixels[seq_along(sizes)] == tip_pixels[length(sizes) + seq_along(sizes)]

  # Examine both roles for each hit. Only the stem may be extrapolated.
  stem <- c(hits[, 1L], hits[, 2L])
  cross <- c(hits[, 2L], hits[, 1L])
  point <- rbind(hits[, 3:4, drop = FALSE], hits[, 3:4, drop = FALSE])
  len <- sqrt(rowSums((b - a)^2))
  direction <- (b - a) / len
  along_stem <- rowSums((point - a[stem, , drop = FALSE]) * direction[stem, , drop = FALSE])
  along_cross <- rowSums((point - a[cross, , drop = FALSE]) * direction[cross, , drop = FALSE])
  at_start <- along_stem < -1.5 & stem == first[contour[stem]] & free[contour[stem]]
  at_end <- along_stem > len[stem] + 1.5 & stem == last[contour[stem]] &
    free[length(sizes) + contour[stem]]
  # The target must continue through the hit, with at least one pixel to
  # spare when its fitted segment ends at an actual contour endpoint.
  eps <- 1e-8
  through <- along_cross >= -eps & along_cross <= len[cross] + eps &
    (closed[contour[cross]] | ((cross != first[contour[cross]] | along_cross > 1) &
                                (cross != last[contour[cross]] | len[cross] - along_cross > 1)))
  transverse <- abs(rowSums(direction[stem, , drop = FALSE] * direction[cross, , drop = FALSE])) <= cos(pi / 6)
  candidates <- which((at_start | at_end) & len[stem] >= 5 & through & transverse)
  if (!length(candidates)) return(mask)
  gap <- ifelse(at_start, -along_stem, along_stem - len[stem])
  tip <- ifelse(at_start, contour[stem], length(sizes) + contour[stem])

  # Each endpoint gets its nearest supported target. Decline equally near
  # alternatives when they require different bridges. All candidates see
  # the same original skeleton, so one addition cannot seed another.
  bridges <- lapply(split(candidates, tip[candidates]), function(rows) {
    best <- NULL
    for (i in rows[order(gap[rows])]) {
      if (!is.null(best) && gap[i] > best_gap + eps) break
      z <- drawing$contours[[contour[cross[i]]]]
      vertices <- rbind(z[, 1:2, drop = FALSE], z[nrow(z), 3:4])
      pixels <- .line_bridge_pixels(tips[tip[i], ], point[i, ], vertices, skeleton, lightness, max_gap)
      if (!length(pixels)) next
      if (!is.null(best) && !identical(best, pixels)) return(NULL)
      best <- pixels
      best_gap <- gap[i]
    }
    best
  })
  mask[unlist(bridges, use.names = FALSE)] <- TRUE
  mask
}

.line_bridge_pixels <- function(start, hit, vertices, skeleton, lightness, max_gap) {
  delta <- hit - start
  gap <- sqrt(sum(delta^2))
  forward <- delta / gap
  normal <- 1.25 * c(-forward[2L], forward[1L])
  at <- function(distance) sweep(outer(distance, forward), 2L, start, "+")
  # Compare signed contrast on the visible stem and across the empty gap.
  # Stop a pixel before the target, where its own edge mixes into the sample.
  xy <- at(c(-5:-1, seq(0.5, gap - 1, length.out = max(5L, ceiling(gap)))))
  contrast <- .line_sample(lightness, sweep(xy, 2L, normal, "+")) -
    .line_sample(lightness, sweep(xy, 2L, normal, "-"))
  if (anyNA(contrast)) return(NULL)
  stem <- stats::median(contrast[1:5])
  if (abs(stem) < 0.035 || any(contrast[-(1:5)] / stem < 0.4)) return(NULL)

  # Snap to a nearby pixel on the crossing contour. Reject ties after
  # preferring the shorter bridge. The cap also applies to the snapped gap.
  d2 <- rowSums(sweep(vertices, 2L, hit)^2)
  nearby <- unique(vertices[d2 <= min(d2) + 1e-8, , drop = FALSE])
  if (min(d2) > 1.5^2) return(NULL)
  distance2 <- rowSums(sweep(nearby, 2L, start)^2)
  closest <- which(distance2 <= min(distance2) + 1e-8)
  if (length(closest) != 1L || distance2[closest] > max_gap^2 + 1e-8) return(NULL)
  h <- nrow(skeleton)
  w <- ncol(skeleton)
  target <- nearby[closest, ]
  pixels <- cpp_raster_indices(matrix(c(start, target), nrow = 1L), h, w, 1)
  on_path <- arrayInd(pixels, dim(skeleton))[, 2:1, drop = FALSE]
  # A pre-existing stroke in the middle of the proposed bridge makes the
  # target ambiguous. Allow its own endpoint neighborhoods.
  interior <- rowSums(sweep(on_path, 2L, start)^2) > 2 &
    rowSums(sweep(on_path, 2L, target)^2) > 2
  if (any(skeleton[pixels[interior]])) return(NULL)
  added <- !skeleton[pixels]
  if (!any(added) || !.line_bridge_contacts(skeleton, on_path[added, , drop = FALSE], start, target))
    return(NULL)
  pixels[added]
}

# Connected patches of original ink touched by the new pixels. Touching a
# third patch by an edge or a diagonal can create an extra arm or a loop,
# even when the bridge itself passes through empty pixels.
.line_bridge_contacts <- function(skeleton, added, start, target) {
  offsets <- as.matrix(expand.grid(dx = -1:1, dy = -1:1))
  neighbors <- unique(added[rep(seq_len(nrow(added)), each = nrow(offsets)), , drop = FALSE] +
                        offsets[rep(seq_len(nrow(offsets)), times = nrow(added)), ])
  h <- nrow(skeleton)
  w <- ncol(skeleton)
  x <- neighbors[, 1L]
  y <- neighbors[, 2L]
  rc <- cbind(y, x)[x >= 1 & x <= w & y >= 1 & y <= h, , drop = FALSE]
  ink <- rc[skeleton[rc], , drop = FALSE]
  if (nrow(ink) < 2L) return(FALSE)
  groups <- cpp_cluster_points(ink, sqrt(2))
  ends <- rbind(start, target)
  contacts <- match(ends[, 2L] + (ends[, 1L] - 1) * h, ink[, 1L] + (ink[, 2L] - 1) * h)
  length(unique(groups)) == 2L && !anyNA(contacts) && groups[contacts[1L]] != groups[contacts[2L]]
}

.line_sample <- function(image, xy) {
  if (any(xy[, 1L] < 1 | xy[, 1L] > ncol(image) |
          xy[, 2L] < 1 | xy[, 2L] > nrow(image))) return(NA_real_)
  x <- floor(xy[, 1L])
  y <- floor(xy[, 2L])
  xn <- pmin(x + 1, ncol(image))
  yn <- pmin(y + 1, nrow(image))
  fx <- xy[, 1L] - x
  fy <- xy[, 2L] - y
  (1 - fx) * (1 - fy) * image[cbind(y, x)] + fx * (1 - fy) * image[cbind(y, xn)] +
    (1 - fx) * fy * image[cbind(yn, x)] + fx * fy * image[cbind(yn, xn)]
}
