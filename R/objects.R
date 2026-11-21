#' Construct a vectorized line drawing
#'
#' Contours are ordered matrices with columns x1, y1, x2, y2. Coordinates
#' are continuous image coordinates: x increases right, y increases down;
#' raster pixel centers begin at (1,1). Geometry may extend outside the image.
#' Zero-length segments require explicit cleaning; a segment whose endpoint
#' coordinates differ by at most 64 units in the last place of their own
#' magnitude (rounding residue), or whose length underflows to zero, counts as
#' zero-length. When such segments are dropped, the next kept segment starts
#' exactly where the previous one ends (cyclically for closed contours), and
#' cleaning repeats until no zero-length segment remains. Subpixel values are
#' otherwise retained.
#' @param contours List of numeric matrices with four columns.
#' @param width,height Positive integer image dimensions.
#' @param ids Unique character contour identifiers.
#' @param source Optional provenance metadata.
#' @param zero_length Reject or drop zero-length segments.
#' @return An S3 `mlv_drawing` with geometry and provenance.
#' @export
#' @examples
#' d <- mlv_drawing(list(rbind(c(1,1,11,1), c(11,1,11,6))), 20, 20)
#' as.data.frame(d)
mlv_drawing <- function(contours = list(), width, height, ids = NULL,
                        source = NULL, zero_length = c("error", "drop")) {
  zero_length <- match.arg(zero_length)
  .dimension <- function(z, name) {
    if (length(z) != 1L || !is.numeric(z) || !is.finite(z) ||
        z < 1 || z != floor(z) || z > .Machine$integer.max)
      stop(name, " must be a positive integer", call. = FALSE)
    as.integer(z)
  }
  width <- .dimension(width, "width"); height <- .dimension(height, "height")
  if (!is.list(contours) || is.data.frame(contours))
    stop("contours must be a list of segment matrices", call. = FALSE)
  if (is.null(ids)) ids <- as.character(seq_along(contours))
  if (length(ids) != length(contours) || anyNA(ids) ||
      any(!nzchar(as.character(ids))) || anyDuplicated(as.character(ids)))
    stop("ids must be unique, nonmissing contour identifiers", call. = FALSE)
  out <- vector("list", length(contours)); keep <- logical(length(contours))
  for (i in seq_along(contours)) {
    z <- contours[[i]]
    if (!is.matrix(z) || !is.numeric(z) || ncol(z) != 4L || any(!is.finite(z)))
      stop("each contour must be a finite numeric matrix with four columns", call. = FALSE)
    storage.mode(z) <- "double"
    # Dropping or snapping can itself create a zero-length segment (a closed
    # loop that collapses to one segment, or a gap snap that moves a start
    # point onto its end), so clean until no zero-length segment remains.
    while (nrow(z)) {
      len <- sqrt(rowSums((z[,3:4,drop=FALSE] - z[,1:2,drop=FALSE])^2))
      if (any(!is.finite(len))) stop("coordinate differences overflow", call. = FALSE)
      # Endpoints that differ only by rounding residue (64 ulp of their own
      # magnitude, with no absolute floor so small-scale geometry stays
      # valid) are the same point: such segments are zero-length (MLV-R009).
      # A segment whose length underflows to zero is zero-length as well.
      zero <- rowSums(abs(z[,3:4,drop=FALSE] - z[,1:2,drop=FALSE]) >
                      64*.Machine$double.eps*pmax(abs(z[,1:2,drop=FALSE]),abs(z[,3:4,drop=FALSE]))) == 0 |
              len == 0
      if (any(zero)) {
        if (zero_length == "error") stop("zero-length segment; use zero_length='drop' to clean", call. = FALSE)
        n0 <- nrow(z); kept <- which(!zero)
        closed0 <- n0 > 1L && all(abs(z[1L,1:2] - z[n0,3:4]) <=
                                  64*.Machine$double.eps*pmax(1,abs(z[1L,1:2]),abs(z[n0,3:4])))
        z <- z[kept,,drop=FALSE]
        # Dropping a segment leaves its residue as a gap between neighbours;
        # the next kept segment starts exactly where the previous one ends,
        # cyclically for closed contours, so connectivity and closure stay
        # exact and do not depend on the start vertex (MLV-R016).
        for (m in seq_along(kept)[-1L]) if (kept[m]-kept[m-1L]>1L) z[m,1:2] <- z[m-1L,3:4]
        if (length(kept) && closed0 && (kept[1L]>1L || kept[length(kept)]<n0))
          z[1L,1:2] <- z[length(kept),3:4]
        next
      }
      if (nrow(z)>1L && any(abs(z[-1L,1:2,drop=FALSE] -
                                   z[-nrow(z),3:4,drop=FALSE]) >
                            64*.Machine$double.eps*pmax(1,abs(z[-1L,1:2,drop=FALSE]),
                                                        abs(z[-nrow(z),3:4,drop=FALSE]))))
        stop("segments within a contour must be connected and ordered", call. = FALSE)
      # Gaps within the residue tolerance are made exact, including the closing
      # gap, so closure tests downstream do not depend on the start vertex or on
      # whether a segment was dropped (MLV-R016).
      if (nrow(z)>1L) {
        n1 <- nrow(z); z0 <- z
        z[-1L,1:2] <- z[-n1,3:4]
        if (all(abs(z[1L,1:2] - z[n1,3:4]) <=
                64*.Machine$double.eps*pmax(1,abs(z[1L,1:2]),abs(z[n1,3:4]))))
          z[1L,1:2] <- z[n1,3:4]
        if (!identical(z, z0)) next
      }
      break
    }
    colnames(z) <- c("x1","y1","x2","y2")
    out[[i]] <- z; keep[i] <- nrow(z)>0L
  }
  structure(list(contours=out[keep], width=width, height=height,
                 ids=as.character(ids)[keep], source=source), class="mlv_drawing")
}

.check_drawing <- function(x) {
  if (!inherits(x, "mlv_drawing")) stop("x must be an mlv_drawing", call. = FALSE)
  x
}
.rebuild <- function(x, contours=x$contours, ids=x$ids, operation=NULL) {
  out<-mlv_drawing(contours, x$width, x$height, ids,
              source=list(parent=x$source, operation=operation), zero_length="drop")
  if(!is.null(x$feature_spec))out$feature_spec<-x$feature_spec
  if(!is.null(x$curvature_estimator))out$curvature_estimator<-x$curvature_estimator
  out
}
.empty_segments <- function() data.frame(contour_id=character(), contour_index=integer(),
  segment_index=integer(), x1=double(), y1=double(), x2=double(), y2=double(),
  length=double(), orientation=double(), curvature=double())
.segment_table <- function(x) {
  .check_drawing(x)
  if (!length(x$contours)) return(.empty_segments())
  sizes <- vapply(x$contours,nrow,integer(1))
  z <- do.call(rbind,x$contours)
  dx <- z[,3]-z[,1]; dy <- z[,4]-z[,2]
  lens <- sqrt(dx^2+dy^2)
  ori <- (atan2(-dy,dx)*180/pi) %% 360
  n <- nrow(z); ends <- cumsum(sizes)
  neighbor <- seq_len(n)+1L
  neighbor[ends] <- pmax(ends-1L,1L)
  neighbor[ends[sizes==1L]] <- ends[sizes==1L]
  delta <- abs(ori-ori[neighbor])
  curv <- pmin(delta,360-delta)/lens
  data.frame(contour_id=x$ids[rep.int(seq_along(sizes),sizes)],
             contour_index=rep.int(seq_along(sizes),sizes),
             segment_index=sequence(sizes),z,length=lens,orientation=ori,
             curvature=curv,row.names=NULL)
}
#' @export
print.mlv_drawing <- function(x, ...) {
  cat("<mlv_drawing>", length(x$contours), "contours,",
      sum(vapply(x$contours,nrow,integer(1))), "segments,",
      x$width,"x",x$height,"pixels\n")
  invisible(x)
}
#' Drawing segment table
#' @param x A drawing.
#' @param row.names,optional,... Standard data frame method arguments.
#' @param spec Optional feature specification; defaults to an enriched drawing's recipe.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification; see
#'   [compute_contour_properties()]; exclusive of `spec`.
#' @return A data frame with one row per input segment: endpoints, length,
#'   orientation and curvature. Under the turning and arc-length estimators
#'   curvature is the segment's mean curvature (the exact integral of the
#'   estimator's density over the segment divided by its length) and
#'   `signed_curvature` its signed mean (NA where an exact reversal leaves the
#'   sign undefined); [curvature_profile()] returns the piecewise density itself.
#'   Lengths and curvature use the specification's units.
#' @export
as.data.frame.mlv_drawing <- function(x, row.names=NULL, optional=FALSE, ...,
                                    spec=x$feature_spec,pixels_per_degree=NULL,
                                    estimator=c("turning","mlv2023")) {
  estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
  if(!is.null(spec)&&estimator_given)stop("use spec or estimator, not both",call.=FALSE)
  .measurement_table(x,spec,pixels_per_degree,estimator)
}
#' Select contours and invalidate derived properties
#' @param x A drawing.
#' @param i Integer, logical, or character contour selection.
#' @return A drawing containing the selected contours.
#' @export
subset_contours <- function(x,i) {
  .check_drawing(x)
  if(is.character(i)) i <- match(i,x$ids)
  ii <- seq_along(x$contours)[i]
  if(anyNA(ii) || anyDuplicated(ii)) stop("invalid or repeated contour selection",call.=FALSE)
  .rebuild(x,x$contours[ii],x$ids[ii],operation="subset")
}
