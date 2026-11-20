# The specification-free curvature estimator of a drawing: an explicit choice
# wins; otherwise the one recorded by compute_contour_properties(estimator=),
# otherwise "turning".
.drawing_estimator <- function(x,estimator="turning",given=FALSE) {
 if(given||is.null(x$curvature_estimator))return(estimator)
 match.arg(x$curvature_estimator,c("turning","mlv2023"))
}
# Piecewise curvature profile (vertex turning or arc-length estimators) used by
# histograms, which integrate the density exactly within each bin. A
# specification overrides `estimator`.
.profile_table <- function(x,spec=NULL,pixels_per_degree=NULL,estimator="turning") {
 if(is.null(spec))return(if(estimator=="mlv2023").segment_table(x) else curvature_profile(x,estimator))
 validate_feature_spec(spec);settings<-spec$settings
 factor<-.spec_scale(spec,pixels_per_degree)
 if(settings$estimator!="mlv2023")return(curvature_profile(x,settings$estimator,
  settings$curvature_scale,settings$units,if(is.null(pixels_per_degree))settings$pixels_per_degree else pixels_per_degree))
 st<-.segment_table(x);st$length<-st$length*factor;st$curvature<-st$curvature/factor;st
}
# One row per input segment. For "mlv2023" curvature is the published
# per-segment value; for "turning" and "arclength" it is the segment's mean
# curvature, the exact integral of the estimator's density over the segment
# divided by its length, and signed_curvature is the signed mean (NA where an
# exact reversal makes the sign undefined). The piecewise density itself is
# available from curvature_profile().
.measurement_table <- function(x,spec=NULL,pixels_per_degree=NULL,estimator=.drawing_estimator(x)) {
 if(!is.null(spec)) {
  validate_feature_spec(spec);settings<-spec$settings
  if(settings$estimator=="mlv2023")return(.profile_table(x,spec,pixels_per_degree))
  return(.segment_means(x,settings$estimator,settings$curvature_scale,settings$units,
   if(is.null(pixels_per_degree))settings$pixels_per_degree else pixels_per_degree))
 }
 if(estimator=="mlv2023").segment_table(x) else .segment_means(x,estimator)
}
.segment_means <- function(x,estimator,scale=1,units="pixels",pixels_per_degree=NULL) {
 st<-.segment_table(x)
 factor<-.curvature_units(units,pixels_per_degree)
 scale<-.spec_number(scale,"scale",strict=TRUE)
 st$length<-st$length*factor
 st$signed_curvature<-numeric(nrow(st))
 for(i in seq_along(x$contours)) {
  rows<-which(st$contour_index==i)
  cp<-.contour_profile(x,i,estimator,scale,factor);p<-cp$profile;L<-cp$total;closed<-cp$curve$closed
  a<-p$s_start;w<-p$s_end-a
  k0<-p$curvature_start;k1<-p$curvature_end
  z0<-p$signed_start;z1<-p$signed_end;na<-is.na(z0)|is.na(z1);z0[na]<-0;z1[na]<-0
  cumulative<-function(v0,v1)c(0,cumsum((v0+v1)/2*w))
  # Integral of the piecewise-affine density from the canonical origin to t.
  primitive<-function(t,C,v0,v1) {
   F<-function(t) {
    j<-pmax(1L,findInterval(t,a));d<-pmin(pmax(t-a[j],0),w[j])
    f<-ifelse(w[j]>0,d/w[j],0)
    C[j]+d*(v0[j]+(v0[j]+(v1[j]-v0[j])*f))/2
   }
   if(closed){q<-floor(t/L);q*C[length(C)]+F(t-q*L)} else F(pmin(pmax(t,0),L))
  }
  # Input arc positions in canonical coordinates. A closed contour's canonical
  # origin is its first vertex that survives the removal of collinear vertices.
  s_in<-c(0,cumsum(st$length[rows]))
  origin<-s_in[cp$curve$kept[1]]
  u<-s_in[-length(s_in)]-origin;v<-s_in[-1]-origin
  integral<-function(v0,v1){C<-cumulative(v0,v1);primitive(v,C,v0,v1)-primitive(u,C,v0,v1)}
  len<-st$length[rows]
  st$curvature[rows]<-pmax(0,integral(k0,k1))/len
  signed<-integral(z0,z1)/len
  undefined<-integral(as.double(na),as.double(na))>0
  signed[undefined]<-NA_real_
  st$signed_curvature[rows]<-signed
 }
 st
}
.measurement_junctions <- function(x,spec=NULL,pixels_per_degree=NULL) {
 if(is.null(spec))return(if(is.null(x$junctions))detect_junctions(x) else x$junctions)
 args<-spec$settings$junction;factor<-.spec_scale(spec,pixels_per_degree)
 for(k in c("absolute_epsilon","merge_radius","endpoint_radius"))args[[k]]<-args[[k]]/factor
 do.call(detect_junctions,c(list(x=x),args))
}
.geometry_snapshot <- function(x)list(contours=x$contours,ids=x$ids,width=x$width,height=x$height)
.checked_mat_measurements <- function(x,spec,pixels_per_degree=NULL) {
 if(is.null(spec))return(invisible(NULL))
 validate_feature_spec(spec)
 if(is.null(x$mat_feature_spec)||is.null(x$mat_geometry))
  stop("MAT measurements lack a feature specification; recompute with analyse_mat(spec=...)",call.=FALSE)
 if(!identical(x$mat_geometry,.geometry_snapshot(x)))stop("stale MAT measurements: drawing geometry changed",call.=FALSE)
 if(!identical(x$mat_feature_spec,.mat_semantics(spec,pixels_per_degree)))
  stop("MAT measurement specification mismatch; recompute analyse_mat()",call.=FALSE)
 invisible(NULL)
}
.mat_semantics <- function(spec,pixels_per_degree=NULL) {
 validate_feature_spec(spec);s<-spec$settings
 list(branch_estimator=s$branch_estimator,branch_definition=spec$definitions[['branch']],
  units=s$units,physical_pixel_size=.spec_scale(spec,pixels_per_degree),window=s$branch_window,
  power=s$branch_power,reference_radius=s$reference_radius,fit_tolerance=s$fit_tolerance,
  mat=s$mat,mat_definitions=spec$definitions[c("aof","morphology","branch_rule")],mapping_radius=s$mapping_radius,raster_line_width=s$raster_line_width)
}
.linear_histogram <- function(left,right,lengths,groups,edges,n_groups) {
 k<-length(edges)-1L;out<-matrix(0,n_groups,k)
 if(!length(left)||!n_groups)return(out)
 lo<-pmin(left,right);hi<-pmax(left,right);constant<-hi==lo
 edge<-edges;edge[1]<--Inf;edge[length(edge)]<-Inf
 const_bin<-.bin_index(lo,edges)
 for(j in seq_len(k)) {
  mass<-numeric(length(lo));mass[constant & const_bin==j]<-lengths[constant & const_bin==j]
  moving<-which(!constant)
  mass[moving]<-lengths[moving]*pmax(0,pmin(hi[moving],edge[j+1L])-pmax(lo[moving],edge[j]))/(hi[moving]-lo[moving])
  out[,j]<-.grouped_bins(rep(0,length(mass)),mass,groups,c(-Inf,Inf),n_groups)[,1]
 }
 out
}
.measurement_metadata <- function(x,blocks,normalize,spec=NULL,pixels_per_degree=NULL,by="contour",estimator="turning") {
 coverage<-lapply(intersect(vapply(blocks,`[[`,character(1),"property"),names(x$mat_properties)),function(p) {
  pieces<-x$mat_properties[[p]]
  data.frame(contour_id=x$ids,property=p,mapped=vapply(pieces,function(v)sum(is.finite(v)),integer(1)),
   sampled=lengths(pieces),coverage=vapply(pieces,function(v)if(length(v))mean(is.finite(v)) else NA_real_,double(1)))
 })
 list(schema=.feature_schema_from_blocks(blocks,normalize,spec,estimator),aggregation=by,
  calibration=if(is.null(spec))NULL else list(units=spec$settings$units,
   pixels_per_degree=if(is.null(pixels_per_degree))spec$settings$pixels_per_degree else pixels_per_degree),
  input_provenance=x$source,coverage=if(length(coverage))do.call(rbind,coverage) else NULL)
}
