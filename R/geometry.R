#' Compute contour measurements and junctions
#' @param x A drawing.
#' @param properties Requested properties; basic segment measures are cheap and
#'   computed together. Include "junctions" to detect contour intersections.
#' @param ... Arguments passed to `detect_junctions()`.
#' @param spec Optional feature specification, exclusive of properties, dots
#'   and `estimator`.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification. Both give
#'   one row per input segment. `"turning"` (default) spreads each vertex turn
#'   over its dual cell, counting it once; a segment's curvature is its mean
#'   curvature, the exact integral of that density over the segment divided by
#'   its length, and `signed_curvature` is the signed mean. The piecewise density
#'   is available from [curvature_profile()]. `"mlv2023"` gives the published
#'   value: the turn to the next segment divided by the segment's own length,
#'   with the last segment repeating the previous turn. Without a
#'   specification the choice is recorded on the drawing
#'   (`curvature_estimator`) and becomes the default of later calls on it that
#'   give neither `spec` nor `estimator` (statistics, feature tables, segment
#'   tables, plots, densities and contour splits); explicit arguments still win.
#' @return A drawing with segment and contour tables, and requested junctions.
#' @export
compute_contour_properties <- function(x,
  properties=c("orientation","length","curvature","junctions"), ...,
  spec=x$feature_spec,pixels_per_degree=NULL,estimator=c("turning","mlv2023")) {
  .check_drawing(x);estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
  if(!is.null(spec)) {
    validate_feature_spec(spec)
    if(!missing(properties)||length(list(...))||estimator_given)stop("use spec or individual contour settings, not both",call.=FALSE)
    if(!is.null(pixels_per_degree))spec$settings$pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
    properties<-c(intersect(spec$settings$properties,c("orientation","length","curvature")),
                  if(any(grepl("^junction_",spec$settings$properties)))"junctions")
  }
  if(any(!properties %in% c("orientation","length","curvature","junctions")))
    stop("unknown contour property",call.=FALSE)
  x$segment_properties <- .measurement_table(x,spec,pixels_per_degree,estimator)
  st <- x$segment_properties
  x$contour_properties <- data.frame(contour_id=x$ids,
    length=if(length(x$ids))unname(vapply(split(st$length,
      factor(st$contour_index,levels=seq_along(x$ids))),sum,double(1))) else numeric())
  if("junctions" %in% properties) x$junctions <- if(is.null(spec))detect_junctions(x,...) else
    .measurement_junctions(x,spec,pixels_per_degree)
  if(!is.null(spec))x$feature_spec<-spec
  x$curvature_estimator<-if(is.null(spec))estimator else NULL
  x
}
.empty_junctions <- function() {
  z <- data.frame(junction_id=integer(),x=double(),y=double(),type=character(),
                  min_angle=double(),max_angle=double())
  z$angles <- I(list()); z$arms <- I(list()); z$members <- I(list()); z
}
#' Detect and classify intersections between distinct contours
#'
#' Returns T, Y, X, Arrow and Star junctions. Same-contour self crossings and
#' collinear overlaps have no unique intersection and are excluded. Proximate
#' intersections merge by connected components, including transitive chains.
#'
#' For the default `rule = "exit"`, arms are defined on the junction
#' neighbourhood: the union of open discs of radius `endpoint_radius` around
#' every intersection merged into the junction, so each merged crossing keeps
#' its arms. Each incident contour contributes one arm for each direction in
#' which it leaves that region, following the contour past vertices inside
#' it; a contour that ends inside the region contributes no arm in that
#' direction. A segment that never enters the region and lies wholly beyond
#' its intersections (a gap wider than the radius) contributes only its arm
#' leaving away from them. Closed contours wrap around. Overshoots, gaps and
#' contour pieces shorter than the radius are therefore part of the junction,
#' and the result does not depend on contour direction, numbering or start
#' vertex. An arm's direction is the orientation of the segment through which
#' it leaves, pointing away from the junction, so angles between straight
#' arms do not depend on the merged position. Junctions with fewer than three
#' arms are dropped. An endpoint extension exactly equal to its limit is
#' accepted. Angles use the drawing convention: counterclockwise degrees,
#' image y down.
#' @param x A drawing.
#' @param absolute_epsilon Maximum endpoint extension in pixels, measured
#'   according to `gap_metric`.
#' @param relative_epsilon Maximum endpoint extension as a segment fraction;
#'   the stricter of the two limits applies.
#' @param merge_radius Radius for merging nearby intersections.
#' @param endpoint_radius Radius of the junction neighbourhood that arms must leave.
#' @param method Exact sweep candidate filtering or exhaustive pair enumeration.
#' @param gap_metric How `absolute_epsilon` is measured along a segment.
#'   `"euclidean"` (default) measures length and is rotation invariant.
#'   `"chebyshev"` (as in MLV) counts the larger coordinate change, so an
#'   extension reaches one raster step: up to `sqrt(2)` pixels along a
#'   diagonal.
#' @param rule Junction definition. `"exit"` (default) is the definition
#'   described in Details. `"mlv2023"` reproduces the MLV MATLAB toolbox
#'   (commit f12418a) including its recorded defects (MLV-M006 to M008, M019,
#'   M021 to M026 in UPSTREAM_BUGS.md): endpoint snapping, discarded negative
#'   coordinates, order-dependent merging, start-first arm attachment and
#'   adjacent-segment arms, no wrapping of closed contours, and its
#'   interior-arm angle convention. It always uses the Chebyshev gap and the
#'   exhaustive pair order; `method` is ignored and `gap_metric` may only be
#'   `"chebyshev"`. Its members list the adjacent segments MLV adds. Use it
#'   only to reproduce published MLV values.
#' @return A data frame with one row per junction: `junction_id`, position
#'   `x`, `y`, `type`, `min_angle` and `max_angle` (degrees), and three list
#'   columns:
#'   \describe{
#'   \item{`arms`}{The orientation of each arm, in degrees from 0 to 360,
#'     counterclockwise from the positive x axis in the drawing convention
#'     (image y pointing down, so 90 points up the screen), sorted
#'     increasing. Two arms can share an orientation.}
#'   \item{`angles`}{The gaps between consecutive arms. With `n` arms,
#'     `angles[k] = arms[k + 1] - arms[k]` for `k < n`, and the last gap wraps
#'     around: `angles[n] = arms[1] + 360 - arms[n]`. Under the default rule
#'     the gaps sum to 360.}
#'   \item{`members`}{The incident contour and segment indices.}
#'   }
#'
#'   Under `rule = "mlv2023"`, `arms` are the directions MLV itself assigns,
#'   not directions in the drawing convention. A segment that starts or ends
#'   within `endpoint_radius` of the junction position contributes its
#'   orientation pointing away from that end, and its neighbouring segment in
#'   the contour contributes one more arm (MLV-M024). A segment that passes
#'   the junction contributes, for each of its two ends, the direction from
#'   that end towards the junction position computed with y pointing up
#'   (MLV-M007). On diagonal segments these arms are therefore
#'   mirrored, and they can coincide with other arms. `angles` are computed
#'   from these arms exactly as above, modulo 360.
#' @export
detect_junctions <- function(x,absolute_epsilon=1,relative_epsilon=0.3,
                            merge_radius=2,endpoint_radius=2,
                            method=c("sweep","exhaustive"),
                            gap_metric=c("euclidean","chebyshev"),
                            rule=c("exit","mlv2023")) {
  .check_drawing(x); method <- match.arg(method); rule <- match.arg(rule)
  explicit_metric <- !missing(gap_metric); gap_metric <- match.arg(gap_metric)
  for(v in list(absolute_epsilon,relative_epsilon,merge_radius,endpoint_radius))
    if(length(v)!=1L || !is.finite(v) || v<0) stop("junction tolerances must be finite and nonnegative",call.=FALSE)
  if(rule=="mlv2023") {
    if(explicit_metric && gap_metric!="chebyshev")
      stop("rule = \"mlv2023\" uses the Chebyshev gap; gap_metric must be \"chebyshev\"",call.=FALSE)
    return(.junctions_mlv2023(x,absolute_epsilon,relative_epsilon,merge_radius,endpoint_radius))
  }
  st <- .segment_table(x)
  if(nrow(st)<2L) return(.empty_junctions())
  raw <- cpp_intersections(as.matrix(st[,c("x1","y1","x2","y2")]),
                            st$contour_index,absolute_epsilon,relative_epsilon,method=="sweep",
                            gap_metric=="euclidean")
  if(!nrow(raw)) return(.empty_junctions())
  groups <- cpp_cluster_points(raw[,3:4,drop=FALSE],merge_radius)
  result <- list()
  for(g in unique(groups)) {
    rr <- raw[groups==g,,drop=FALSE]; p <- colMeans(rr[,3:4,drop=FALSE])
    indices <- sort(unique(as.integer(rr[,1:2])))
    members <- unique(st[indices,c("contour_index","segment_index"),drop=FALSE])
    rays <- list(); keys <- character(); pts <- unname(rr[,3:4,drop=FALSE])
    # The junction neighbourhood is the union of open discs of radius
    # endpoint_radius around every intersection merged into the junction, so
    # each merged crossing keeps its own arms (MLV-R014). Each incident
    # contour contributes one arm per direction in which it leaves that
    # region; a contour ending inside it contributes none in that direction
    # (MLV-R006). The rule is symmetric in contour direction.
    r2 <- max(endpoint_radius^2,.Machine$double.eps)
    inside <- function(v) any((pts[,1]-v[1])^2+(pts[,2]-v[2])^2 < r2)
    add_arm <- function(z, ci, si, direction) {
      # Closed contours (last end == first start) wrap around, so the result
      # does not depend on which vertex starts the polygon (MLV-R007).
      n <- nrow(z); j <- si; closed <- n>1L && all(z[1,1:2]==z[n,3:4])
      repeat {
        if(!inside(if(direction>0) z[j,3:4] else z[j,1:2])) break
        if(!closed && ((direction>0 && j==n) || (direction<0 && j==1L))) return(invisible(NULL))
        j <- if(closed) (j-1L+direction)%%n+1L else j+direction
        if(j==si) return(invisible(NULL)) # the whole closed contour lies inside
      }
      key <- paste(ci,j,direction,sep=":")
      if(!key %in% keys) {
        keys <<- c(keys,key)
        # The arm points along the segment through which the contour leaves
        # the neighbourhood. Using that segment, not the chord from the merged
        # junction position, keeps straight-arm angles exact (MLV-R005).
        d <- unname(z[j,3:4]-z[j,1:2])*direction
        rays[[length(rays)+1L]] <<- (atan2(-d[2],d[1])*180/pi)%%360
      }
    }
    for(i in seq_len(nrow(members))) {
      ci <- members$contour_index[i]; si <- members$segment_index[i]
      z <- x$contours[[ci]]; a <- unname(z[si,1:2]); u <- unname(z[si,3:4])-a
      # Parameters of the intersections this segment takes part in. A segment
      # that never enters the neighbourhood and lies wholly beyond those
      # intersections (a gap wider than the radius) leaves only away from
      # them; it has no arm pointing back across the junction (MLV-R008).
      g <- which(st$contour_index==ci & st$segment_index==si)
      own <- pts[rr[,1]==g | rr[,2]==g,,drop=FALSE]
      t <- ((own[,1]-a[1])*u[1]+(own[,2]-a[2])*u[2])/sum(u^2)
      tc <- pmin(pmax(((pts[,1]-a[1])*u[1]+(pts[,2]-a[2])*u[2])/sum(u^2),0),1)
      reaches <- any((a[1]+tc*u[1]-pts[,1])^2+(a[2]+tc*u[2]-pts[,2])^2 < r2)
      if(reaches || all(t<=0) || any(t>0 & t<1)) add_arm(z,ci,si,1L)
      if(reaches || all(t>=1) || any(t>0 & t<1)) add_arm(z,ci,si,-1L)
    }
    ori <- sort(unlist(rays)); n <- length(ori)
    if(n<3L) next
    angles <- (c(ori[-1L],ori[1L]+360)-ori)
    mx <- max(angles)
    type <- if(n==3L) {if(mx<160) "Y" else if(mx<=200) "T" else "Arrow"} else if(n==4L) "X" else "Star"
    z <- data.frame(junction_id=length(result)+1L,x=p[1],y=p[2],type=type,
                    min_angle=min(angles),max_angle=mx)
    members$contour_id <- x$ids[members$contour_index]
    z$angles <- I(list(angles)); z$arms <- I(list(ori)); z$members <- I(list(members))
    result[[length(result)+1L]] <- z
  }
  if(!length(result)) .empty_junctions() else do.call(rbind,result)
}
#' Mean geometric or medial-axis property
#' @param x A drawing, with medial-axis properties when requested.
#' @param property Property name.
#' @param spec Optional feature specification.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification; see
#'   [compute_contour_properties()]. Mean curvature is integrated absolute
#'   turning per unit length.
#' @return One numeric summary, or NA when undefined.
#' @export
average_property <- function(x,property=c("length","curvature","orientation","junctions",
                                         "mirror","parallelism","separation","taper"),
                             spec=x$feature_spec,pixels_per_degree=NULL,
                             estimator=c("turning","mlv2023")) {
  .check_drawing(x); property <- match.arg(property)
  estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
  if(!is.null(spec)&&estimator_given)stop("use spec or estimator, not both",call.=FALSE)
  # Without a specification, only curvature depends on the estimator; the other
  # summaries measure the input segments directly.
  st <- .measurement_table(x,spec,pixels_per_degree,if(property=="curvature")estimator else "mlv2023")
  total <- sum(st$length)
  if(property=="length") {
    if(!length(x$ids)) return(NA_real_)
    return(total/length(x$ids))
  }
  if(total==0) return(NA_real_)
  if(property=="curvature") return(sum(st$curvature*st$length)/total)
  if(property=="orientation") {
    # Explicit axial circular mean. An isotropic resultant has no direction.
    a <- st$orientation*pi/90
    sx <- sum(st$length*cos(a)); sy <- sum(st$length*sin(a))
    if(sqrt(sx^2+sy^2) <= 64*.Machine$double.eps*total) return(NA_real_)
    return((atan2(sy,sx)*90/pi)%%180)
  }
  if(property=="junctions") return(nrow(.measurement_junctions(x,spec,pixels_per_degree))/total*10000)
  .checked_mat_measurements(x,spec,pixels_per_degree)
  values <- x$mat_properties[[property]]
  if(is.null(values)) stop("requested MAT property has not been mapped",call.=FALSE)
  values <- unlist(values,use.names=FALSE)
  if(!any(is.finite(values))) NA_real_ else mean(values[is.finite(values)])
}
