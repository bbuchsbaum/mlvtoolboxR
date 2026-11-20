.measurement_table <- function(x,spec=NULL,pixels_per_degree=NULL) {
 if(is.null(spec))return(.segment_table(x))
 validate_feature_spec(spec);settings<-spec$settings
 factor<-.spec_scale(spec,pixels_per_degree)
 if(settings$estimator!="mlv2023")return(curvature_profile(x,settings$estimator,
  settings$curvature_scale,settings$units,if(is.null(pixels_per_degree))settings$pixels_per_degree else pixels_per_degree))
 st<-.segment_table(x);st$length<-st$length*factor;st$curvature<-st$curvature/factor;st
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
 const_bin<-findInterval(lo,edges,all.inside=TRUE)
 for(j in seq_len(k)) {
  mass<-numeric(length(lo));mass[constant & const_bin==j]<-lengths[constant & const_bin==j]
  moving<-which(!constant)
  mass[moving]<-lengths[moving]*pmax(0,pmin(hi[moving],edge[j+1L])-pmax(lo[moving],edge[j]))/(hi[moving]-lo[moving])
  out[,j]<-.grouped_bins(rep(0,length(mass)),mass,groups,c(-Inf,Inf),n_groups)[,1]
 }
 out
}
.measurement_metadata <- function(x,blocks,normalize,spec=NULL,pixels_per_degree=NULL,by="contour") {
 coverage<-lapply(intersect(vapply(blocks,`[[`,character(1),"property"),names(x$mat_properties)),function(p) {
  pieces<-x$mat_properties[[p]]
  data.frame(contour_id=x$ids,property=p,mapped=vapply(pieces,function(v)sum(is.finite(v)),integer(1)),
   sampled=lengths(pieces),coverage=vapply(pieces,function(v)if(length(v))mean(is.finite(v)) else NA_real_,double(1)))
 })
 list(schema=.feature_schema_from_blocks(blocks,normalize,spec),aggregation=by,
  calibration=if(is.null(spec))NULL else list(units=spec$settings$units,
   pixels_per_degree=if(is.null(pixels_per_degree))spec$settings$pixels_per_degree else pixels_per_degree),
  input_provenance=x$source,coverage=if(length(coverage))do.call(rbind,coverage) else NULL)
}
