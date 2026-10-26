.canonical_curve <- function(z) {
 xy<-rbind(z[1,1:2],z[,3:4,drop=FALSE]);closed<-all(xy[1,]==xy[nrow(xy),])
 if(closed)xy<-xy[-nrow(xy),,drop=FALSE]
 repeat {
  n<-nrow(xy);vertices<-if(closed)seq_len(n) else if(n>2)2:(n-1L) else integer()
  drop<-integer()
  for(i in vertices) {
   a<-xy[i,]-xy[if(i==1)n else i-1L,]
   b<-xy[if(i==n)1L else i+1L,]-xy[i,]
   a<-a/sqrt(sum(a*a));b<-b/sqrt(sum(b*b))
   cross<-a[1]*b[2]-a[2]*b[1];dot<-sum(a*b)
   if(abs(cross)<=64*.Machine$double.eps) {
    if(dot<0)stop("180-degree reversal has undefined signed turning",call.=FALSE)
    if(dot>0)drop<-c(drop,i)
   }
  }
  if(!length(drop))break
  xy<-xy[-drop,,drop=FALSE]
 }
 if(closed)xy<-rbind(xy,xy[1,])
 delta<-xy[-1,,drop=FALSE]-xy[-nrow(xy),,drop=FALSE]
 lengths<-sqrt(rowSums(delta^2));s<-c(0,cumsum(lengths));n<-length(lengths)
 vertices<-if(closed)seq_len(n) else if(n>1L)2:n else integer()
 previous<-ifelse(vertices==1,n,vertices-1L)
 turn<-vapply(seq_along(vertices),function(j) {
  a<-delta[previous[j],]/lengths[previous[j]];b<-delta[vertices[j],]/lengths[vertices[j]]
  atan2(a[2]*b[1]-a[1]*b[2],sum(a*b))*180/pi
 },double(1))
 list(xy=xy,lengths=lengths,s=s,total=tail(s,1),closed=closed,
      vertices=vertices,previous=previous,turn=turn)
}
.curvature_units <- function(units,pixels_per_degree) {
 units<-match.arg(units,c("pixels","visual_degrees"))
 if(units=="pixels")return(1)
 if(is.null(pixels_per_degree))stop("visual_degrees requires pixels_per_degree",call.=FALSE)
 1/.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
}
#' Exact signed and absolute vertex turning
#' @param x A drawing.
#' @param units Spatial units for arc positions; angles remain geometric degrees.
#' @param pixels_per_degree Required calibration for visual-degree positions.
#' @return A table of canonical vertices, positions and turns. x/y coordinates
#'   remain in input pixels. Collinear forward vertices are removed; exact
#'   retracing is rejected because its signed turn is undefined.
#' @export
contour_turning <- function(x,units=c("pixels","visual_degrees"),pixels_per_degree=NULL) {
 .check_drawing(x);units<-match.arg(units);factor<-.curvature_units(units,pixels_per_degree)
 empty<-data.frame(contour_id=character(),contour_index=integer(),vertex_index=integer(),
  x=double(),y=double(),arclength=double(),signed_turn=double(),absolute_turn=double())
 out<-lapply(seq_along(x$contours),function(i) {
  c<-.canonical_curve(x$contours[[i]]);v<-c$vertices
  if(!length(v))return(empty)
  data.frame(contour_id=x$ids[i],contour_index=i,vertex_index=v,
   x=c$xy[v,1],y=c$xy[v,2],arclength=c$s[v]*factor,
   signed_turn=c$turn,absolute_turn=abs(c$turn))
 })
 if(length(out))do.call(rbind,out) else empty
}
.profile_cuts <- function(x,L) {
 x<-sort(unique(pmax(0,pmin(L,x))))
 keep<-c(TRUE,diff(x)>64*.Machine$double.eps*L)
 x<-x[keep];x[1]<-0;x[length(x)]<-L;x
}
.tent_values <- function(at,positions,turn,h,L,closed) {
 out<-numeric(length(at))
 for(i in seq_along(positions)) {
  distance<-abs(at-positions[i])
  if(closed) {distance<-pmin(distance,L-distance);area<-h} else {
   a<-min(h,positions[i]);b<-min(h,L-positions[i])
   # Stable even when h is much larger than the entire open contour.
   area<-a-a*(a/h)/2+b-b*(b/h)/2
  }
  out<-out+turn[i]*pmax(0,1-distance/h)/area
 }
 if(any(!is.finite(out)))stop("curvature exceeds the supported numerical range",call.=FALSE)
 out
}
.empty_curvature_profile <- function()data.frame(contour_id=character(),contour_index=integer(),
 segment_index=integer(),x1=double(),y1=double(),x2=double(),y2=double(),length=double(),
 orientation=double(),curvature=double(),signed_curvature=double(),curvature_start=double(),
 curvature_end=double(),signed_start=double(),signed_end=double(),s_start=double(),s_end=double())
#' Piecewise exact polygon turning-density profiles
#'
#' Turning uses constant densities on vertex dual cells. Arclength uses
#' source-normalized triangular kernels with physical half-width `scale`.
#' Open-end kernels retain their turning mass; closed kernels are periodic with
#' half-width at most half the contour length. This smooths turning density,
#' not curve coordinates. Opposite turns can cancel in the smoothed field.
#' @param x A drawing.
#' @param estimator Arclength tent density or vertex dual-cell density.
#' @param scale Positive kernel half-width in target spatial units.
#' @param units Spatial units for lengths and inverse-curvature dimensions.
#' @param pixels_per_degree Required for visual-degree units.
#' @return A table partitioned at geometric vertices, kernel knots and signed
#'   zero crossings. x/y coordinates remain in input pixels. Curvature start/end
#'   values are unsigned; signed_start/end retain direction. Midpoint values
#'   integrate means exactly, but histograms require affine interval integration.
#' @export
curvature_profile <- function(x,estimator=c("arclength","turning"),scale=1,
 units=c("pixels","visual_degrees"),pixels_per_degree=NULL) {
 .check_drawing(x);estimator<-match.arg(estimator);units<-match.arg(units)
 scale<-.spec_number(scale,"scale",strict=TRUE);factor<-.curvature_units(units,pixels_per_degree)
 out<-lapply(seq_along(x$contours),function(i) {
  curve<-.canonical_curve(x$contours[[i]])
  s<-curve$s*factor;L<-curve$total*factor
  if(!is.finite(L)||L<=0)stop("contour length exceeds the supported numerical range",call.=FALSE)
  positions<-s[curve$vertices];turn<-curve$turn
  cuts<-c(0,L,s)
  if(estimator=="arclength") {
   h<-if(curve$closed)min(scale,L/2) else scale
   knots<-c(positions-h,positions,positions+h)
   if(curve$closed)knots<-knots%%L
   cuts<-.profile_cuts(c(cuts,knots),L)
   values<-.tent_values(cuts,positions,turn,h,L,curve$closed)
   left<-head(values,-1);right<-tail(values,-1)
   cross<-which(sign(left)*sign(right)<0)
   if(length(cross)) {
    zeros<-cuts[cross]+diff(cuts)[cross]*(-left[cross])/(right[cross]-left[cross])
    cuts<-.profile_cuts(c(cuts,zeros),L)
    values<-.tent_values(cuts,positions,turn,h,L,curve$closed)
   }
   v0<-head(values,-1);v1<-tail(values,-1)
  } else {
   lo<-positions-curve$lengths[curve$previous]*factor/2
   hi<-positions+curve$lengths[curve$vertices]*factor/2
   knots<-c(lo,hi);if(curve$closed)knots<-knots%%L
   cuts<-.profile_cuts(c(cuts,knots),L)
   mid<-(head(cuts,-1)+tail(cuts,-1))/2
   v0<-numeric(length(mid))
   for(j in seq_along(turn)) {
    inside<-if(curve$closed)(mid-lo[j])%%L<hi[j]-lo[j] else mid>lo[j]&mid<hi[j]
    v0[inside]<-v0[inside]+turn[j]/(hi[j]-lo[j])
   }
   v1<-v0
  }
  a<-head(cuts,-1);b<-tail(cuts,-1);mid<-(a+b)/2
  index<-findInterval(mid,s,all.inside=TRUE)
  point<-function(t) cbind(stats::approx(s,curve$xy[,1],xout=t,rule=2)$y,
                           stats::approx(s,curve$xy[,2],xout=t,rule=2)$y)
  pa<-point(a);pb<-point(b);delta<-curve$xy[index+1,,drop=FALSE]-curve$xy[index,,drop=FALSE]
  data.frame(contour_id=x$ids[i],contour_index=i,segment_index=seq_along(a),
   x1=pa[,1],y1=pa[,2],x2=pb[,1],y2=pb[,2],length=b-a,
   orientation=(atan2(-delta[,2],delta[,1])*180/pi)%%360,
   curvature=(abs(v0)+abs(v1))/2,signed_curvature=(v0+v1)/2,
   curvature_start=abs(v0),curvature_end=abs(v1),signed_start=v0,signed_end=v1,
   s_start=a,s_end=b)
 })
 result<-if(length(out))do.call(rbind,out) else .empty_curvature_profile()
 attr(result,"definition")<-list(estimator=estimator,version=if(estimator=="arclength")
  "tent-source-normalized-1" else "vertex-dual-1",scale=scale,units=units)
 result
}
