#' Rotate a drawing in image coordinates
#' @param x A drawing.
#' @param angle Counterclockwise angle in degrees.
#' @param center Two-coordinate rotation center.
#' @return A drawing with derived properties invalidated.
#' @export
rotate_drawing <- function(x,angle,center=c(x$width,x$height)/2) {
 .check_drawing(x)
 if(length(angle)!=1L||!is.finite(angle)||length(center)!=2L||any(!is.finite(center)))
  stop("invalid rotation",call.=FALSE)
 a<-angle*pi/180
 transform<-function(z) {
  for(cols in list(1:2,3:4)) {
   p<-sweep(z[,cols,drop=FALSE],2L,center,"-")
   z[,cols]<-cbind(cos(a)*p[,1]+sin(a)*p[,2],-sin(a)*p[,1]+cos(a)*p[,2])
   z[,cols]<-sweep(z[,cols,drop=FALSE],2L,center,"+")
  };z
 }
 .rebuild(x,lapply(x$contours,transform),operation=list(rotate=angle,center=center))
}
#' Shift contours, optionally sampling translations within the image
#' @param x A drawing.
#' @param offsets A contour-by-two translation matrix, or NULL for random shifts.
#' @param within Require translated contours to fit inside the image.
#' @return A drawing; R's ordinary random seed controls random translations.
#' @export
shift_contours <- function(x,offsets=NULL,within=is.null(offsets)) {
 .check_drawing(x);n<-length(x$ids)
 if(is.null(offsets)) {
  offsets<-matrix(0,n,2)
  for(i in seq_len(n)) {
   z<-x$contours[[i]];lo<-c(1-min(z[,c(1,3)]),1-min(z[,c(2,4)]))
   hi<-c(x$width-max(z[,c(1,3)]),x$height-max(z[,c(2,4)]))
   if(any(hi<lo))stop("contour cannot fit inside the image",call.=FALSE)
   offsets[i,]<-stats::runif(2,lo,hi)
  }
 }
 if(!is.matrix(offsets)||!is.numeric(offsets)||!identical(dim(offsets),c(n,2L))||any(!is.finite(offsets)))
  stop("offsets must be a finite contour-by-two matrix",call.=FALSE)
 contours<-lapply(seq_len(n),function(i)sweep(x$contours[[i]],2L,rep(offsets[i,],2L),"+"))
 if(within&&any(vapply(contours,function(z)any(z[,c(1,3)]<1|z[,c(1,3)]>x$width) ||
      any(z[,c(2,4)]<1|z[,c(2,4)]>x$height),logical(1))))
  stop("translated contour is outside the image",call.=FALSE)
 .rebuild(x,contours,operation=list(shift=offsets))
}
.piece_drawing <- function(x,pieces,parents,operation) {
 if(!length(pieces))return(.rebuild(x,list(),character(),operation))
 ids<-make.unique(paste0(parents,".part"),sep=".")
 out<-.rebuild(x,pieces,ids,operation);out$parent_ids<-parents;out
}
#' Clip all segment-circle intersections
#'
#' Handles segments with both endpoints outside but crossing the aperture.
#' Adjacent retained segments remain in the same contour.
#' @param x A drawing.
#' @param radius Positive aperture radius.
#' @param center Aperture center.
#' @return A drawing with parent contour IDs and no derived-property cache.
#' @export
clip_circle <- function(x,radius=min(x$width,x$height)/2,center=c(x$width,x$height)/2) {
 .check_drawing(x)
 if(length(radius)!=1L||!is.finite(radius)||radius<=0||length(center)!=2L||any(!is.finite(center)))
  stop("invalid aperture",call.=FALSE)
 pieces<-list();parents<-character()
 for(i in seq_along(x$contours)) {
  z<-x$contours[[i]];current<-matrix(numeric(),0,4)
  flush<-function() {
   if(nrow(current)){pieces[[length(pieces)+1L]]<<-current;parents<<-c(parents,x$ids[i]);current<<-matrix(numeric(),0,4)}
  }
  for(s in seq_len(nrow(z))) {
   a<-z[s,1:2];v<-z[s,3:4]-a;d<-a-center
   length<-sqrt(sum(v^2));unit<-v/length
   along<-sum((center-a)*unit);foot<-a+along*unit
   perpendicular<-sum((foot-center)^2)
   if(perpendicular>=radius^2){flush();next}
   half<-sqrt(max(0,radius^2-perpendicular))
   lower<-along-half;upper<-along+half
   lo<-max(0,lower);hi<-min(length,upper)
   if(hi<=lo){flush();next}
   first<-if(lo==0)a else foot-half*unit
   last<-if(hi==length)a+v else foot+half*unit
   piece<-c(first,last)
   if(nrow(current)&&sqrt(sum((current[nrow(current),3:4]-piece[1:2])^2))>1e-9)flush()
   current<-rbind(current,piece)
   if(hi<length)flush()
  };flush()
 }
 .piece_drawing(x,pieces,parents,list(aperture=list(radius=radius,center=center)))
}
#' Simplify ordered contours by line fitting
#' @param x A drawing.
#' @param tolerance Maximum mean point-to-line distance.
#' @return A drawing with the same contour IDs.
#' @export
simplify_contours <- function(x,tolerance=1) {
 .check_drawing(x)
 contours<-lapply(x$contours,function(z) {
  xy<-rbind(z[,1:2,drop=FALSE],z[nrow(z),3:4])
  spans<-fit_line_segments(xy,tolerance)$spans
  cbind(xy[spans[,1],,drop=FALSE],xy[spans[,2],,drop=FALSE])
 })
 .rebuild(x,contours,operation=list(simplify=tolerance))
}
#' Split contours at detected junction positions
#' @param x A drawing.
#' @param junctions Optional junction data frame.
#' @return A drawing whose pieces retain parent IDs and total segment geometry.
#' @export
segment_at_junctions <- function(x,junctions=NULL) {
 .check_drawing(x)
 if(is.null(junctions))junctions<-detect_junctions(x)
 pieces<-list();parents<-character()
 for(i in seq_along(x$contours)) {
  z<-x$contours[[i]];current<-matrix(numeric(),0,4)
  flush<-function() {
   if(nrow(current)){pieces[[length(pieces)+1L]]<<-current;parents<<-c(parents,x$ids[i]);current<<-matrix(numeric(),0,4)}
  }
  for(s in seq_len(nrow(z))) {
   a<-z[s,1:2];v<-z[s,3:4]-a;times<-numeric()
   for(j in seq_len(nrow(junctions))) {
    m<-junctions$members[[j]]
    if(any(m$contour_index==i & m$segment_index==s)) {
     t<-sum((c(junctions$x[j],junctions$y[j])-a)*v)/sum(v^2)
     if(t>=-1e-8&&t<=1+1e-8)times<-c(times,pmax(0,pmin(1,t)))
    }
   }
   cuts<-sort(unique(c(0,times,1)))
   if(0 %in% times)flush()
   for(k in seq_len(length(cuts)-1L)) {
    lo<-cuts[k];hi<-cuts[k+1L]
    if(hi>lo)current<-rbind(current,c(a+lo*v,a+hi*v))
    if(hi %in% times)flush()
   }
  };flush()
 }
 .piece_drawing(x,pieces,parents,"segment_at_junctions")
}
.arc_slice <- function(z,from,to) {
 lens<-sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2));total<-sum(lens)
 start<-c(0,head(cumsum(lens),-1));end<-cumsum(lens)
 out<-matrix(numeric(),0,4)
 for(s in seq_len(nrow(z))) {
  lo<-max(start[s],from*total);hi<-min(end[s],to*total)
  if(hi>lo) {
   v<-z[s,3:4]-z[s,1:2]
   out<-rbind(out,c(z[s,1:2]+(lo-start[s])/lens[s]*v,z[s,1:2]+(hi-start[s])/lens[s]*v))
  }
 };out
}
#' Separate middle halves from end quarters
#' @param x A drawing.
#' @return List with middle and junction-region drawings. Each first/last
#' quarter is defined by arc length after segmentation at junctions.
#' @export
split_junction_regions <- function(x) {
 x<-segment_at_junctions(x)
 middle<-lapply(x$contours,.arc_slice,from=.25,to=.75)
 ends<-unlist(lapply(x$contours,function(z)list(.arc_slice(z,0,.25),.arc_slice(z,.75,1))),recursive=FALSE)
 list(middle=.piece_drawing(x,middle,x$ids,"middle_half"),
      junction=.piece_drawing(x,ends,rep(x$ids,each=2),"end_quarters"))
}
.contour_criterion <- function(x,property) {
 st<-.measurement_table(x,x$feature_spec)
 ends<-cumsum(tabulate(st$contour_index,nbins=length(x$ids)))
 starts<-c(1L,head(ends,-1L)+1L)
 junctions<-if(property=="junctions") {
  .measurement_junctions(x,x$feature_spec)
 } else NULL
 vapply(seq_along(x$ids),function(i) {
  s<-st[starts[i]:ends[i],,drop=FALSE];ll<-sum(s$length)
  switch(property,
   length=ll,curvature=sum(s$curvature*s$length)/ll,
   integrated_curvature=sum(s$curvature*s$length),
   orientation=sum((abs(cos(s$orientation*pi/180))-abs(sin(s$orientation*pi/180)))*s$length)/ll,
   junctions=sum(vapply(junctions$members,function(m)i %in% m$contour_index,logical(1))),
   random=stats::runif(1),
   {.checked_mat_measurements(x,x$feature_spec);z<-x$mat_properties[[property]][[i]]
    if(is.null(z))stop("missing contour criterion: ",property,call.=FALSE)
    # A contour with no defined mapped score has no criterion value (NA).
    if(any(is.finite(z)))mean(z[is.finite(z)]) else NA_real_})
 },double(1))
}
#' Split ranked contours under a whole-contour length budget
#'
#' Each group is a ranked prefix from one end that fits within fraction times
#' total contour length. A contour that would exceed the budget stops that
#' prefix. Ties use contour IDs. Curvature ranks by a length-weighted mean;
#' integrated_curvature explicitly exposes the upstream length-weighted sum.
#' Medial-axis criteria rank by the mean of a contour's finite mapped scores. A
#' contour with no finite mapped score (no defined branch score reaches any of
#' its pixels) cannot be ranked: it is left unassigned, its score is NA, and
#' it is listed in `unranked`. Ranks are taken among the rankable contours
#' only; the length budget remains a fraction of the total length of all
#' contours.
#' @param x A drawing.
#' @param property One or more criterion names.
#' @param fraction Length fraction in [0,1] for each group.
#' @param weights Criterion rank weights.
#' @param scores Optional finite per-contour scores, bypassing rank combination.
#' @return List with top, bottom, unassigned drawings, overlap IDs, scores
#'   (NA for unrankable contours) and `unranked`, the IDs of contours without
#'   a criterion value (always empty when `scores` is supplied).
#' @export
split_contours <- function(x,property="length",fraction=.5,weights=NULL,scores=NULL) {
 .check_drawing(x);n<-length(x$ids)
 if(length(fraction)!=1L||!is.finite(fraction)||fraction<0||fraction>1)
  stop("fraction must be in [0,1]",call.=FALSE)
 if(is.null(scores)) {
  if(!length(property))stop("supply at least one property",call.=FALSE)
  if(is.null(weights))weights<-rep(1,length(property))
  if(length(weights)!=length(property)||any(!is.finite(weights)))stop("invalid weights",call.=FALSE)
  values<-lapply(property,function(p).contour_criterion(x,p))
  rankable<-Reduce(`&`,lapply(values,is.finite))
  scores<-numeric(n)
  for(p in seq_along(property)) {
   ord<-order(values[[p]][rankable],x$ids[rankable])
   ranks<-numeric(sum(rankable));ranks[ord]<-seq_len(sum(rankable))
   scores[rankable]<-scores[rankable]+weights[p]*ranks
  }
  scores[!rankable]<-NA_real_
 } else {
  if(!is.numeric(scores)||length(scores)!=n||any(!is.finite(scores)))stop("invalid scores",call.=FALSE)
  rankable<-rep(TRUE,n)
 }
 lens<-.contour_criterion(x,"length");budget<-fraction*sum(lens)
 take<-function(order) {
  if(fraction==1)return(order)
  if(!length(order))return(integer())
  order[cumsum(lens[order])<=budget+64*.Machine$double.eps*max(1,budget)]
 }
 ok<-which(rankable)
 bottom<-take(ok[order(scores[ok],x$ids[ok])]);top<-take(ok[order(-scores[ok],x$ids[ok])])
 unused<-setdiff(seq_len(n),union(top,bottom))
 list(top=subset_contours(x,top),bottom=subset_contours(x,bottom),
      unassigned=subset_contours(x,unused),overlap=x$ids[intersect(top,bottom)],
      scores=stats::setNames(scores,x$ids),unranked=x$ids[!rankable])
}
#' Split contours using weighted feature histograms
#' @param x A drawing.
#' @param weights Named numeric vector matching feature-table columns.
#' @param properties,bins Feature settings.
#' @param fraction Whole-contour length budget.
#' @return The same structure as split_contours.
#' @export
split_by_histogram <- function(x,weights,properties=c("orientation","length","curvature"),
                               bins=8,fraction=.5) {
 f<-feature_table(x,properties,bins,normalize=FALSE)
 if(is.null(names(weights))||any(!names(weights)%in% names(f))||
    !is.numeric(weights)||any(!is.finite(weights)))stop("weights must name feature-table columns",call.=FALSE)
 score<-if(length(weights))as.numeric(as.matrix(f[,names(weights),drop=FALSE])%*%weights) else rep(0,nrow(f))
 split_contours(x,fraction=fraction,scores=score)
}
