#' Compute contour measurements and junctions
#' @param x A drawing.
#' @param properties Requested properties; basic segment measures are cheap and
#'   computed together. Include "junctions" to detect contour intersections.
#' @param ... Arguments passed to `detect_junctions()`.
#' @return A drawing with segment and contour tables, and requested junctions.
#' @export
compute_contour_properties <- function(x,
  properties=c("orientation","length","curvature","junctions"), ...) {
  .check_drawing(x)
  if(any(!properties %in% c("orientation","length","curvature","junctions")))
    stop("unknown contour property",call.=FALSE)
  x$segment_properties <- .segment_table(x)
  st <- x$segment_properties
  x$contour_properties <- data.frame(contour_id=x$ids,
    length=if(length(x$ids))unname(vapply(split(st$length,
      factor(st$contour_index,levels=seq_along(x$ids))),sum,double(1))) else numeric())
  if("junctions" %in% properties) x$junctions <- detect_junctions(x,...)
  x
}
.empty_junctions <- function() {
  z <- data.frame(junction_id=integer(),x=double(),y=double(),type=character(),
                  min_angle=double(),max_angle=double())
  z$angles <- I(list()); z$members <- I(list()); z
}
#' Detect and classify intersections between distinct contours
#'
#' Returns T, Y, X, Arrow and Star junctions. Same-contour self crossings and
#' collinear overlaps have no unique intersection and are excluded. Proximate
#' intersections merge by connected components, including transitive chains.
#' @param x A drawing.
#' @param absolute_epsilon Maximum endpoint extension in pixels.
#' @param relative_epsilon Maximum endpoint extension as a segment fraction.
#' @param merge_radius Radius for merging nearby intersections.
#' @param endpoint_radius Radius for treating adjacent segments as junction arms.
#' @param method Exact sweep candidate filtering or exhaustive pair enumeration.
#' @return A data frame with positions, types, angles and incident segment lists.
#' @export
detect_junctions <- function(x,absolute_epsilon=1,relative_epsilon=0.3,
                            merge_radius=2,endpoint_radius=2,
                            method=c("sweep","exhaustive")) {
  .check_drawing(x); method <- match.arg(method)
  for(v in list(absolute_epsilon,relative_epsilon,merge_radius,endpoint_radius))
    if(length(v)!=1L || !is.finite(v) || v<0) stop("junction tolerances must be finite and nonnegative",call.=FALSE)
  st <- .segment_table(x)
  if(nrow(st)<2L) return(.empty_junctions())
  raw <- cpp_intersections(as.matrix(st[,c("x1","y1","x2","y2")]),
                            st$contour_index,absolute_epsilon,relative_epsilon,method=="sweep")
  if(!nrow(raw)) return(.empty_junctions())
  groups <- cpp_cluster_points(raw[,3:4,drop=FALSE],merge_radius)
  result <- list()
  for(g in unique(groups)) {
    rr <- raw[groups==g,,drop=FALSE]; p <- colMeans(rr[,3:4,drop=FALSE])
    indices <- sort(unique(as.integer(rr[,1:2])))
    members <- unique(st[indices,c("contour_index","segment_index"),drop=FALSE])
    rays <- list(); keys <- character()
    add_ray <- function(ci, si, direction) {
      z <- x$contours[[ci]]; n <- nrow(z); j <- si
      repeat {
        target <- if(direction>0) z[j,3:4] else z[j,1:2]
        if(sum((target-p)^2)>max(endpoint_radius^2, .Machine$double.eps) ||
           (direction>0 && j==n) || (direction<0 && j==1L)) break
        j <- j+direction
      }
      if(sum((target-p)^2)<.Machine$double.eps) return(invisible(NULL))
      key <- paste(ci,j,direction,sep=":")
      if(!key %in% keys) {
        keys <<- c(keys,key)
        rays[[length(rays)+1L]] <<- (atan2(p[2]-target[2],target[1]-p[1])*180/pi)%%360
      }
    }
    for(i in seq_len(nrow(members))) {
      ci <- members$contour_index[i]; si <- members$segment_index[i]
      z <- x$contours[[ci]]; a <- z[si,1:2]; b <- z[si,3:4]
      da <- sum((a-p)^2); db <- sum((b-p)^2)
      if(da<endpoint_radius^2 && da<=db) {
        add_ray(ci,si,1L)
        if(si>1L) add_ray(ci,si-1L,-1L)
      } else if(db<endpoint_radius^2) {
        add_ray(ci,si,-1L)
        if(si<nrow(z)) add_ray(ci,si+1L,1L)
      } else {
        add_ray(ci,si,-1L); add_ray(ci,si,1L)
      }
    }
    ori <- sort(unlist(rays)); n <- length(ori)
    if(n<3L) next
    angles <- (c(ori[-1L],ori[1L]+360)-ori)
    mx <- max(angles)
    type <- if(n==3L) {if(mx<160) "Y" else if(mx<=200) "T" else "Arrow"} else if(n==4L) "X" else "Star"
    z <- data.frame(junction_id=length(result)+1L,x=p[1],y=p[2],type=type,
                    min_angle=min(angles),max_angle=mx)
    members$contour_id <- x$ids[members$contour_index]
    z$angles <- I(list(angles)); z$members <- I(list(members))
    result[[length(result)+1L]] <- z
  }
  if(!length(result)) .empty_junctions() else do.call(rbind,result)
}
#' Mean geometric or medial-axis property
#' @param x A drawing, with medial-axis properties when requested.
#' @param property Property name.
#' @return One numeric summary, or NA when undefined.
#' @export
average_property <- function(x,property=c("length","curvature","orientation","junctions",
                                         "mirror","parallelism","separation","taper")) {
  .check_drawing(x); property <- match.arg(property)
  st <- .segment_table(x); total <- sum(st$length)
  if(property=="length") {
    if(!length(x$ids)) return(NA_real_)
    return(mean(vapply(x$contours,function(z)sum(sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))),double(1))))
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
  if(property=="junctions") return(nrow(detect_junctions(x))/total*10000)
  values <- x$mat_properties[[property]]
  if(is.null(values)) stop("requested MAT property has not been mapped",call.=FALSE)
  values <- unlist(values,use.names=FALSE)
  if(!any(is.finite(values))) NA_real_ else mean(values[is.finite(values)])
}
