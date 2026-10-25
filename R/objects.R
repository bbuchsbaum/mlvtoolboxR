#' Construct a vectorized line drawing
#'
#' Contours are ordered matrices with columns x1, y1, x2, y2. Coordinates
#' are continuous image coordinates: x increases right, y increases down;
#' raster pixel centers begin at (1,1). Geometry may extend outside the image.
#' Zero-length segments require explicit cleaning; subpixel values are retained.
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
    if (nrow(z)) {
      len <- sqrt(rowSums((z[,3:4,drop=FALSE] - z[,1:2,drop=FALSE])^2))
      if (any(!is.finite(len))) stop("coordinate differences overflow", call. = FALSE)
      if (any(len == 0)) {
        if (zero_length == "error") stop("zero-length segment; use zero_length='drop' to clean", call. = FALSE)
        z <- z[len > 0,,drop=FALSE]
      }
      if (nrow(z)>1L && any(abs(z[-1L,1:2,drop=FALSE] -
                                   z[-nrow(z),3:4,drop=FALSE]) >
                            64*.Machine$double.eps*pmax(1,abs(z[-1L,1:2,drop=FALSE]),
                                                        abs(z[-nrow(z),3:4,drop=FALSE]))))
        stop("segments within a contour must be connected and ordered", call. = FALSE)
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
  mlv_drawing(contours, x$width, x$height, ids,
              source=list(parent=x$source, operation=operation), zero_length="drop")
}
.empty_segments <- function() data.frame(contour_id=character(), contour_index=integer(),
  segment_index=integer(), x1=double(), y1=double(), x2=double(), y2=double(),
  length=double(), orientation=double(), curvature=double())
.segment_table <- function(x) {
  .check_drawing(x)
  if (!length(x$contours)) return(.empty_segments())
  rows <- lapply(seq_along(x$contours), function(i) {
    z <- x$contours[[i]]; n <- nrow(z)
    dx <- z[,3]-z[,1]; dy <- z[,4]-z[,2]
    lens <- sqrt(dx^2+dy^2)
    ori <- (atan2(-dy,dx)*180/pi) %% 360
    curv <- numeric(n)
    if(n>1L) {
      delta <- abs(ori - ori[c(2:n,n-1L)])
      curv <- pmin(delta,360-delta)/lens
    }
    data.frame(contour_id=x$ids[i],contour_index=i,segment_index=seq_len(n),
               z,length=lens,orientation=ori,curvature=curv,row.names=NULL)
  })
  do.call(rbind,rows)
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
#' @return A segment data frame with geometric measurements.
#' @export
as.data.frame.mlv_drawing <- function(x, row.names=NULL, optional=FALSE, ...) {
  .segment_table(x)
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
