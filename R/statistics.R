.bins_count <- function(bins) {
 if(length(bins)!=1L || !is.finite(bins) || bins<1 || bins!=floor(bins) ||
    bins>.Machine$integer.max) stop("bins must be a positive integer",call.=FALSE)
 as.integer(bins)
}
.empty_statistics <- function() data.frame(contour_id=character(),property=character(),
 bin=integer(),label=character(),lower=double(),upper=double(),center=double(),
 value=double(),units=character())
.grouped_bins <- function(values,weights,groups,edges,n_groups,index=NULL) {
 k<-length(edges)-1L
 if(as.double(n_groups)*k>.Machine$integer.max)stop("histogram exceeds index capacity",call.=FALSE)
 if(!length(values)||!n_groups)return(matrix(0,n_groups,k))
 if(is.null(index))index<-.bin_index(values,edges)
 cpp_grouped_sums(weights,groups,index,n_groups,k)
}
# Bin-edge snapping. A value within .bin_edge_tolerance(edge) of a bin edge is
# replaced by the edge itself before binning, so that values that differ from
# an edge only by floating-point rounding bin as if exactly on it (and the
# upper-bin rule then applies). The tolerance is relative, 1e-9 * max(1, |edge|),
# in the coordinate in which the bins are equally spaced (degrees, log10(1 + x),
# or the raw score). It is far above accumulated rounding of the angle and
# length computations (~1e-13 relative) and far below any measurement
# difference that matters. src/statistics.cpp implements the same rule with the
# same arithmetic; the two must agree exactly.
.bin_edge_tolerance <- function(edges) 1e-9*pmax(1,abs(edges))
.snap_to_edges <- function(values,edges) {
 if(!length(values)||!length(edges))return(values)
 i<-findInterval(values,edges);k<-length(edges)
 for(candidate in list(pmax(i,1L),pmin(i+1L,k))) {
  e<-edges[candidate]
  hit<-!is.na(values)&is.finite(e)&abs(values-e)<=.bin_edge_tolerance(e)
  values[hit]<-e[hit]
 }
 values
}
.bin_index <- function(values,edges) findInterval(.snap_to_edges(values,edges),edges,all.inside=TRUE)
# Axial orientation bins are centred on 0 degrees: edges at (m - 1/2) * 180 / bins.
# Scale the axial angle directly. Shifting/subtracting half a bin can put 90
# and 270 degrees on different sides of the same edge for odd bin counts.
.orientation_index <- function(values,bins) {
 axial<-values%%180;scaled<-axial*bins/180+.5
 nearest<-floor(scaled+.5);edge<-(nearest-.5)*180/bins
 bin<-floor(scaled);snap<-!is.na(axial)&abs(axial-edge)<=.bin_edge_tolerance(edge)
 bin[snap]<-nearest[snap]
 (bin%%bins)+1L
}
.statistics_blocks <- function(x,properties,bins,by,normalize,length_range,curvature_range,
                               spec=NULL,pixels_per_degree=NULL,estimator="turning") {
 .check_drawing(x);bins<-.bins_count(bins)
 allowed<-c("orientation","length","curvature","horizontal_vertical",
            "junction_type","junction_angle","mirror","parallelism","separation","taper")
 if(any(!properties %in% allowed))stop("unknown statistic property",call.=FALSE)
 for(r in list(length_range,curvature_range))
  if(length(r)!=2L||any(!is.finite(r))||r[1]<0||r[2]<=r[1])
   stop("log-bin ranges must be finite and increasing",call.=FALSE)
 ng<-if(by=="drawing")1L else length(x$ids)
 if(!length(properties))return(list())
 scale<-if(is.null(spec))1 else .spec_scale(spec,pixels_per_degree)
 length_unit<-if(is.null(spec)||spec$settings$units=="pixels")"pixels" else "visual_degrees"
 geometry_names<-c("orientation","length","curvature","horizontal_vertical")
 geometry<-unique(properties[properties %in% geometry_names])
 if(!is.null(spec))estimator<-spec$settings$estimator
 published<-estimator=="mlv2023"
 # The fused native pass measures input segments directly. It serves every
 # geometry property of the published estimator and of the specification-free
 # turning path (whose curvature it integrates over vertex dual cells, code 5);
 # specified turning and arc-length curvature come from the piecewise profile.
 fused<-published||is.null(spec)
 fused_geometry<-geometry
 fast_counts<-NULL;st<-NULL;segment_group<-NULL
 if(fused) {
  edge_list<-lapply(fused_geometry,function(prop)switch(prop,
   orientation=seq(0,180,length.out=bins+1L),
   length=seq(log10(length_range[1]+1),log10(length_range[2]+1),length.out=bins+1L),
   curvature={e<-seq(log10(curvature_range[1]+1),log10(curvature_range[2]+1),length.out=bins+1L)
              if(published)e else 10^e-1},
   horizontal_vertical=seq(-1,1,length.out=bins+1L)))
  codes<-match(fused_geometry,geometry_names);if(!published)codes[fused_geometry=="curvature"]<-5L
  fast<-cpp_geometry_histograms(x$contours,codes,edge_list,by=="drawing",scale)
  lens<-fast$lengths;fast_counts<-stats::setNames(fast$counts,fused_geometry)
 } else {
  lens<-cpp_contour_lengths(x$contours)*scale
  if(any(geometry!="length")) {
   st<-.profile_table(x,spec,pixels_per_degree)
   segment_group<-if(by=="drawing")rep.int(1L,nrow(st)) else st$contour_index
  }
 }
 totals<-if(by=="drawing")sum(lens) else lens
 contour_group<-if(by=="drawing")rep.int(1L,length(lens)) else seq_along(lens)
 if(any(grepl("^junction_",properties))) {
  j<-.measurement_junctions(x,spec,pixels_per_degree)
  membership<-lapply(j$members,function(m)unique(m$contour_index))
  ji<-if(by=="drawing")seq_len(nrow(j)) else rep.int(seq_len(nrow(j)),lengths(membership))
  jg<-if(by=="drawing")rep.int(1L,nrow(j)) else unlist(membership,use.names=FALSE)
 }
 blocks<-lapply(properties,function(prop) {
  labels<-NULL;unit<-length_unit;index<-NULL
  if(prop=="orientation") {
   width<-180/bins;edges<-seq(0,180,length.out=bins+1L)
   values<-st$orientation;weights<-st$length;groups<-segment_group
   index<-.orientation_index(values,bins)
   centers<-seq(0,180-width,length.out=bins)
  } else if(prop=="horizontal_vertical") {
   edges<-seq(-1,1,length.out=bins+1L)
   values<-pmax(-1,pmin(1,abs(cos(st$orientation*pi/180))-abs(sin(st$orientation*pi/180))))
   weights<-st$length;groups<-segment_group;centers<-(head(edges,-1)+tail(edges,-1))/2
  } else if(prop %in% c("length","curvature")) {
   range<-if(prop=="length")length_range else curvature_range
   edges<-seq(log10(range[1]+1),log10(range[2]+1),length.out=bins+1L)
   values<-log10(if(prop=="length")lens+1 else st$curvature+1)
   weights<-if(prop=="length")lens else st$length
   groups<-if(prop=="length")contour_group else segment_group
   centers<-10^((head(edges,-1)+tail(edges,-1))/2)-1
  } else if(grepl("^junction_",prop)) {
   unit<-"count";groups<-jg;weights<-rep.int(1,length(ji))
   if(prop=="junction_type") {
    labels<-c("T","Y","X","Arrow","Star");edges<-seq(.5,5.5)
    values<-match(j$type[ji],labels);centers<-seq_len(5L)
   } else {
    edges<-seq(0,120,length.out=bins+1L);values<-j$min_angle[ji]
    centers<-(head(edges,-1)+tail(edges,-1))/2
   }
  } else {
   .checked_mat_measurements(x,spec,pixels_per_degree)
   if(is.null(x$mat_properties[[prop]]))stop("MAT property not mapped: ",prop,call.=FALSE)
   pieces<-x$mat_properties[[prop]][seq_along(x$ids)]
   values<-unlist(pieces,use.names=FALSE)
   groups<-if(by=="drawing")rep.int(1L,length(values)) else rep.int(seq_along(pieces),lengths(pieces))
   finite<-is.finite(values);values<-values[finite];groups<-groups[finite]
   weights<-rep.int(1,length(values));edges<-seq(0,1,length.out=bins+1L)
   centers<-(head(edges,-1)+tail(edges,-1))/2;unit<-"count"
  }
  counts<-if(fused && prop %in% fused_geometry)fast_counts[[prop]] else
   if(prop=="curvature"&&!published)
   .linear_histogram(st$curvature_start,st$curvature_end,st$length,groups,10^edges-1,ng) else
   .grouped_bins(values,weights,groups,edges,ng,index)
  if(normalize) {
   counts<-counts/totals*10000
   counts[totals==0,]<-0
  }
  shown_edges<-if(prop %in% c("length","curvature"))10^edges-1 else edges
  if(prop=="orientation")shown_edges<-edges-90/bins
  list(counts=counts,property=prop,bin=seq_len(ncol(counts)),
       label=if(is.null(labels))as.character(signif(centers,8)) else labels,
       lower=head(shown_edges,-1),upper=tail(shown_edges,-1),center=centers,
       measurement_units=if(prop %in% c("orientation","junction_angle"))"degrees" else
        if(prop=="length")length_unit else if(prop=="curvature")paste0("degrees/",length_unit) else "dimensionless",
       units=if(normalize)paste0(unit," per 10000 contour ",length_unit) else unit)
 })
 blocks
}
#' Histograms of drawing and contour properties
#'
#' Orientation, curvature and horizontal/vertical histograms are weighted by
#' segment length. Length histograms collect total contour length. Junctions
#' count events; per-contour junctions count each participating contour once.
#' MAT histograms count mapped pixels, including valid zero scores. Values
#' outside the logarithmic range enter the first/last bin. Orientation is axial;
#' exact bin boundaries enter the upper bin, including odd bin counts.
#'
#' A value that differs from a bin edge by at most 1e-9 x max(1, |edge|) is
#' binned as if it lay exactly on the edge, in the coordinate in which the bins
#' are equally spaced (degrees for orientation and junction angles, log10(1 + x)
#' for length and the published curvature, the raw value otherwise). Angles that
#' are mathematically equal, such as a computed 59.999999999999993 and an exact
#' 60 degrees, therefore share a bin on every platform. Rounding residue is
#' many orders of magnitude smaller than this tolerance, and real measurement
#' differences many orders larger.
#'
#' Without a specification, curvature uses the vertex turning estimator
#' (`estimator = "turning"`): each vertex turn is spread uniformly over its dual
#' cell, half of each adjacent segment, so every turn is counted once. The
#' published estimator, `estimator = "mlv2023"`, assigns each segment the turn
#' to its neighbour divided by its own length and so counts every corner twice;
#' see [feature_spec()].
#' @param x A drawing.
#' @param properties Property names.
#' @param bins Positive number of bins (junction types always have five).
#' @param by Summarize the drawing or each contour.
#' @param normalize Scale by contour length to units per 10,000 pixels.
#' @param length_range,curvature_range Increasing nonnegative log-bin ranges.
#' @param spec Optional feature specification, exclusive of individual feature
#'   settings. Enriched drawings use their saved specification by default.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification: vertex
#'   turning (default) or the published `"mlv2023"` estimator.
#' @return A long data frame with bin boundaries, centers, values and units.
#' @export
contour_statistics <- function(x,properties=c("orientation","length","curvature",
 "horizontal_vertical","junction_type","junction_angle"),bins=8,
 by=c("drawing","contour"),normalize=FALSE,length_range=c(2,max(3,x$width+x$height)),
 curvature_range=c(0,90),spec=x$feature_spec,pixels_per_degree=NULL,
 estimator=c("turning","mlv2023")) {
 by<-match.arg(by);estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(properties)||!missing(bins)||!missing(normalize)||!missing(length_range)||!missing(curvature_range)||
     estimator_given)
   stop("use spec or individual feature settings, not both",call.=FALSE)
  s<-spec$settings;properties<-s$properties;bins<-s$bins;normalize<-s$normalize
  length_range<-s$length_range;curvature_range<-s$curvature_range;estimator<-s$estimator
 }
 blocks<-.statistics_blocks(x,properties,bins,by,normalize,length_range,curvature_range,spec,pixels_per_degree,estimator)
 if(!length(blocks)) {
  out<-.empty_statistics()
  attr(out,"measurement_metadata")<-.measurement_metadata(x,blocks,normalize,spec,pixels_per_degree,by,estimator)
  return(out)
 }
 counts<-do.call(cbind,lapply(blocks,`[[`,"counts"))
 widths<-vapply(blocks,function(b)ncol(b$counts),integer(1))
 flatten<-function(field)unlist(lapply(blocks,`[[`,field),use.names=FALSE)
 ng<-nrow(counts);per_group<-ncol(counts)
 out<-data.frame(contour_id=rep(if(by=="drawing")NA_character_ else x$ids,each=per_group),
  property=rep(rep(properties,widths),ng),bin=rep(flatten("bin"),ng),
  label=rep(flatten("label"),ng),lower=rep(flatten("lower"),ng),
  upper=rep(flatten("upper"),ng),center=rep(flatten("center"),ng),
  value=as.vector(t(counts)),units=rep(rep(flatten("units"),widths),ng),row.names=NULL)
 attr(out,"measurement_metadata")<-.measurement_metadata(x,blocks,normalize,spec,pixels_per_degree,by,estimator)
 out
}
#' Wide per-contour feature table
#' @param x A drawing.
#' @param properties Properties to include.
#' @param bins Number of bins.
#' @param normalize Normalize by contour length.
#' @param length_range,curvature_range Fixed ranges for logarithmic bins.
#' @param spec Optional feature specification, exclusive of individual settings.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification; see
#'   [contour_statistics()].
#' @return One row per contour, with named numeric predictors and a
#'   `feature_schema` attribute. Inspect it with `feature_schema()`.
#' @export
feature_table <- function(x,properties=c("orientation","length","curvature","junction_type"),
                          bins=8,normalize=TRUE,
                          length_range=c(2,max(3,x$width+x$height)),curvature_range=c(0,90),
                          spec=x$feature_spec,pixels_per_degree=NULL,
                          estimator=c("turning","mlv2023")) {
 estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(properties)||!missing(bins)||!missing(normalize)||!missing(length_range)||!missing(curvature_range)||
     estimator_given)
   stop("use spec or individual feature settings, not both",call.=FALSE)
  s<-spec$settings;properties<-s$properties;bins<-s$bins;normalize<-s$normalize
  length_range<-s$length_range;curvature_range<-s$curvature_range;estimator<-s$estimator
 }
 if(!is.logical(normalize)||length(normalize)!=1L||is.na(normalize))
  stop("normalize must be TRUE or FALSE",call.=FALSE)
 blocks<-.statistics_blocks(x,unique(properties),bins,"contour",normalize,
                            length_range,curvature_range,spec,pixels_per_degree,estimator)
 out<-data.frame(contour_id=x$ids)
 if(length(blocks)) {
 counts<-do.call(cbind,lapply(blocks,function(b) {
  colnames(b$counts)<-paste0(b$property,"_",b$bin);b$counts
 }))
 out<-data.frame(out,counts,check.names=FALSE,row.names=NULL)
 }
 metadata<-.measurement_metadata(x,blocks,normalize,spec,pixels_per_degree,estimator=estimator)
 attr(out,"feature_schema")<-metadata$schema
 attr(out,"measurement_metadata")<-metadata
 out
}
#' Predict one score per contour
#' @param x A drawing.
#' @param model An R model with a predict method, or a function of a feature table.
#' @param properties,bins Feature-table settings.
#' @param ... Additional model prediction arguments.
#' @param expected_spec Training schema, from `feature_schema()`. Required unless
#'   the model was bound with `bind_feature_model()`.
#' @param check_spec Check training/actual measurement compatibility. Setting
#'   FALSE is an explicit, warning-emitting legacy bypass.
#' @param normalize,length_range,curvature_range Feature-table settings.
#' @param spec Optional feature specification, exclusive of individual settings.
#' @param pixels_per_degree Optional input calibration override.
#' @param estimator Curvature definition without a specification; see
#'   [contour_statistics()].
#' @return A finite numeric vector in contour order.
#' @export
predict_contours <- function(x,model,properties=c("orientation","length","curvature","junction_type"),
                             bins=8,...,expected_spec=NULL,check_spec=TRUE,
                             normalize=TRUE,length_range=c(2,max(3,x$width+x$height)),
                             curvature_range=c(0,90),spec=x$feature_spec,pixels_per_degree=NULL,
                             estimator=c("turning","mlv2023")) {
 estimator_given<-!missing(estimator);estimator<-.drawing_estimator(x,match.arg(estimator),estimator_given)
 if(!is.logical(check_spec)||length(check_spec)!=1L||is.na(check_spec))
  stop("check_spec must be TRUE or FALSE",call.=FALSE)
 if(inherits(model,"mlv_feature_model")) {
  if(!is.null(expected_spec)).compare_feature_schemas(model$schema,expected_spec)
  expected_spec<-model$schema;model<-model$model
 }
 if(check_spec&&is.null(expected_spec))
  stop("missing training feature schema; bind_feature_model() or supply expected_spec",call.=FALSE)
 if(!is.null(spec)) {
  if(!missing(properties)||!missing(bins)||!missing(normalize)||!missing(length_range)||!missing(curvature_range)||
     estimator_given)
   stop("use spec or individual feature settings, not both",call.=FALSE)
  f<-feature_table(x,spec=spec,pixels_per_degree=pixels_per_degree)
 } else f <- feature_table(x,properties,bins,normalize,length_range,curvature_range,spec=NULL,estimator=estimator)
 if(check_spec).compare_feature_schemas(expected_spec,feature_schema(f)) else
  warning("unchecked prediction: training feature compatibility was not validated",call.=FALSE)
 predictors <- f[,setdiff(names(f),"contour_id"),drop=FALSE]
 result <- if(is.function(model))model(predictors,...) else stats::predict(model,newdata=predictors,...)
 if(!is.numeric(result) || length(result)!=length(x$ids) || any(!is.finite(result)))
  stop("model must return one finite numeric score per contour",call.=FALSE)
 as.numeric(result)
}
