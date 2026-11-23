.spec_number <- function(x,name,lower=0,strict=FALSE,integer=FALSE,upper=Inf) {
 if(integer)upper<-min(upper,.Machine$integer.max)
 if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||
    (if(strict)x<=lower else x<lower)||x>upper||(integer&&x!=floor(x)))
  stop("invalid ",name,call.=FALSE)
 if(integer)as.integer(x) else as.double(x)
}
.spec_options <- function(x,defaults,name) {
 if(!is.list(x)||is.null(names(x))&&length(x)||anyDuplicated(names(x))||
    any(!names(x)%in%names(defaults)))stop("invalid ",name," fields",call.=FALSE)
 defaults[names(x)]<-x
 defaults
}
#' Declare a reusable measurement specification
#'
#' Explicit specifications use fixed ranges, independent of the canvas. Each
#' method's definition is versioned separately from this specification's
#' serialization format and recorded in `definitions`, so a specification states
#' which estimator it uses whether or not that estimator was the default when it
#' was built. Calibration is input-specific and is excluded from semantic
#' compatibility after conversion to the target units.
#'
#' Curvature estimators. `"turning"` (the default) spreads each vertex's turning
#' angle uniformly over the vertex's dual cell, the half-segments on either side
#' of it, so every turn is counted once and a 90 degree corner integrates to 90
#' degrees. It has no tuning parameter, and its histograms are stable under
#' resampling: a circle of radius 10 drawn with 16 to 128 vertices, evenly or
#' with alternating 1:2 spacing, puts all its length in one 8-bin curvature
#' bin, near 180/(10 pi) degrees per pixel, where `"mlv2023"` splits the
#' unevenly sampled circles between two bins. `"arclength"` smooths the turning density with a
#' triangular kernel of half-width `curvature_scale`, giving curvature at a
#' stated spatial scale. `"mlv2023"` is the published MLV estimator: each segment
#' receives the whole turn to the next segment divided by its own length, and
#' the last segment of a contour, which has no next segment, repeats the turn
#' from the previous one. The last vertex of every contour is therefore counted
#' twice (the 90 degree apex of a two-segment roof integrates to 180 degrees),
#' the closing vertex of a closed contour is not counted, and the value at a
#' vertex depends on the length of one adjacent segment only. Use it to
#' reproduce published values (the \href{https://bbuchsbaum.github.io/mlvtoolboxR/articles/mlv-compatibility.html}{MLV compatibility} article).
#' @param properties Ordered feature names.
#' @param bins Histogram bin count; junction types always use five categories.
#' @param normalize Normalize length/count weights per 10,000 target length units.
#' @param length_range,curvature_range Fixed increasing nonnegative bin ranges.
#'   NULL length range means c(2,400) pixels or c(0,20) visual degrees.
#' @param units Spatial unit: pixels or degrees of visual angle. Tangent/turning
#'   angles are always measured in geometric degrees.
#' @param pixels_per_degree Positive calibration required for visual degrees.
#' @param estimator Contour curvature definition: vertex turning (default),
#'   published (`"mlv2023"`), or arc-length-scale curvature; see Details.
#' @param curvature_scale Half-width of the arc-length curvature support, in
#'   target spatial units. Applies only to the arclength estimator.
#' @param branch_estimator Published sample-based or arc-length-based scores.
#' @param branch_window Half-window: samples for mlv2023, target spatial units
#'   for arclength.
#' @param branch_power Positive score exponent.
#' @param reference_radius Reference radius in target spatial units for the
#'   new separation score; published separation retains its one-pixel reference.
#' @param fit_tolerance Straightness tolerance in target spatial units for the
#'   new branch estimator; published mirror fitting retains its one-pixel rule.
#' @param junction Named overrides of absolute_epsilon, relative_epsilon,
#'   merge_radius, endpoint_radius, method, rule (`"exit"` or `"mlv2023"`) and
#'   gap_metric (`"euclidean"` or `"chebyshev"`; defaults to `"chebyshev"` for
#'   rule `"mlv2023"`); see [detect_junctions()]. Spatial tolerances use target units.
#' @param mat Named overrides of threshold_angle, samples, epsilon,
#'   min_component, min_branch_points, aof (`"centered"` or `"mlv2023"`),
#'   morphology (`"minimal"`, `"hmt8"` or `"pratt"`) and branch_rule (`"invariant"` or
#'   `"mlv2023"`); see [compute_mat()]. Epsilon and
#'   component/point counts describe the raster method and remain in
#'   pixels/samples. The morphology also thins extracted contours.
#' @param mapping_radius,raster_line_width Raster mapping distance and ink width
#'   in pixels; rasterization is explicitly a pixel-grid operation.
#' @param extraction Named extraction settings: method, threshold, thin,
#'   min_points, simplify, provider, provider_version and lines. Simplification
#'   is in pixels. `lines` is a named list of [extract_lines()] parameters; the
#'   specification records the line extractor's definition in
#'   `definitions["lines"]` and the validated parameters. The definition is
#'   `"fdog-xdog-1"` by default and `"fdog-xdog-repair-1"` for a positive
#'   `junction_gap`. The provider defaults to `"xdog_lines"`, and [trace_rgb()]
#'   passes the parameters to the backend. A supplied provider version must
#'   match the selected definition.
#' @return A validated `mlv_feature_spec`. It declares methods; constructing it
#'   does not calculate features or infer extraction provenance for existing data.
#' @export
feature_spec <- function(properties=c("orientation","length","curvature","junction_type"),
 bins=8,normalize=TRUE,length_range=NULL,curvature_range=c(0,90),
 units=c("pixels","visual_degrees"),pixels_per_degree=NULL,
 estimator=c("turning","mlv2023","arclength"),curvature_scale=1,
 branch_estimator=c("mlv2023","arclength"),branch_window=5,branch_power=10,
 reference_radius=1,fit_tolerance=1,junction=list(),mat=list(),mapping_radius=4,
 raster_line_width=1,extraction=list()) {
 units<-match.arg(units);estimator<-match.arg(estimator);branch_estimator<-match.arg(branch_estimator)
 allowed<-c("orientation","length","curvature","horizontal_vertical","junction_type","junction_angle",
            "mirror","parallelism","separation","taper")
 if(!is.character(properties)||anyNA(properties)||anyDuplicated(properties)||any(!properties%in%allowed))
  stop("invalid or repeated properties",call.=FALSE)
 bins<-.bins_count(bins)
 if(!is.logical(normalize)||length(normalize)!=1L||is.na(normalize))stop("invalid normalize",call.=FALSE)
 if(is.null(length_range))length_range<-if(units=="pixels")c(2,400) else c(0,20)
 for(r in list(length_range,curvature_range))
  if(!is.numeric(r)||length(r)!=2L||any(!is.finite(r))||r[1]<0||r[2]<=r[1])
   stop("invalid histogram range",call.=FALSE)
 if(units=="visual_degrees"&&is.null(pixels_per_degree))stop("visual_degrees requires pixels_per_degree",call.=FALSE)
 if(!is.null(pixels_per_degree))pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
 curvature_scale<-.spec_number(curvature_scale,"curvature_scale",strict=TRUE)
 branch_window<-.spec_number(branch_window,"branch_window",strict=TRUE,integer=branch_estimator=="mlv2023")
 branch_power<-.spec_number(branch_power,"branch_power",strict=TRUE)
 reference_radius<-.spec_number(reference_radius,"reference_radius",strict=TRUE)
 fit_tolerance<-.spec_number(fit_tolerance,"fit_tolerance",strict=TRUE)
 if(branch_estimator=="mlv2023" && (reference_radius!=1||fit_tolerance!=1))
  stop("mlv2023 branch scores require reference_radius=1 and fit_tolerance=1",call.=FALSE)
 if(units=="visual_degrees"&&branch_estimator=="mlv2023"&&any(properties%in%c("mirror","parallelism","separation","taper")))
  stop("physical branch scores require branch_estimator='arclength'",call.=FALSE)
 junction<-.spec_options(junction,list(absolute_epsilon=1,relative_epsilon=.3,merge_radius=2,
  endpoint_radius=2,method="sweep",rule="exit",gap_metric=NULL),"junction")
 for(k in setdiff(names(junction),c("method","rule","gap_metric")))junction[[k]]<-.spec_number(junction[[k]],paste0("junction.",k))
 junction$method<-match.arg(junction$method,c("sweep","exhaustive"))
 for(k in c("rule","gap_metric"))if(!is.null(junction[[k]])&&(!is.character(junction[[k]])||length(junction[[k]])!=1L||is.na(junction[[k]])))
  stop("junction.",k," must be a single string",call.=FALSE)
 junction$rule<-match.arg(junction$rule,c("exit","mlv2023"))
 if(is.null(junction$gap_metric))junction$gap_metric<-if(junction$rule=="mlv2023")"chebyshev" else "euclidean"
 junction$gap_metric<-match.arg(junction$gap_metric,c("euclidean","chebyshev"))
 if(junction$rule=="mlv2023"&&junction$gap_metric!="chebyshev")
  stop("junction.rule = 'mlv2023' requires junction.gap_metric = 'chebyshev'",call.=FALSE)
 mat<-.spec_options(mat,list(threshold_angle=28,samples=60L,epsilon=1,min_component=NULL,min_branch_points=10L,
  aof="centered",morphology="minimal",branch_rule="invariant"),"mat")
 for(k in c("aof","morphology","branch_rule")) {
  allowed<-switch(k,aof=c("centered","mlv2023"),morphology=c("minimal","hmt8","pratt"),branch_rule=c("invariant","mlv2023"))
  if(!is.character(mat[[k]])||length(mat[[k]])!=1L||is.na(mat[[k]])||!mat[[k]]%in%allowed)
   stop("mat.",k," must be one of ",paste(allowed,collapse=", "),call.=FALSE)
 }
 mat$threshold_angle<-.spec_number(mat$threshold_angle,"threshold_angle",upper=90)
 mat$samples<-.spec_number(mat$samples,"samples",lower=4,upper=100000,integer=TRUE)
 mat$epsilon<-.spec_number(mat$epsilon,"epsilon",strict=TRUE)
 if(!is.null(mat$min_component))mat$min_component<-.spec_number(mat$min_component,"min_component",lower=1,integer=TRUE,upper=.Machine$integer.max)
 mat$min_branch_points<-.spec_number(mat$min_branch_points,"min_branch_points",lower=2,integer=TRUE,upper=.Machine$integer.max)
 .check_centered_samples(mat$aof,mat$samples)
 mapping_radius<-.spec_number(mapping_radius,"mapping_radius",strict=TRUE)
 raster_line_width<-.spec_number(raster_line_width,"raster_line_width",strict=TRUE)
 extraction<-.spec_options(extraction,list(method="pixel_graph_no_corner_cut",threshold=.5,
  thin=TRUE,min_points=2L,simplify=NULL,provider=NULL,provider_version=NULL,lines=NULL),"extraction")
 extraction$threshold<-.spec_number(extraction$threshold,"threshold",lower=-Inf)
 extraction$min_points<-.spec_number(extraction$min_points,"min_points",lower=2,integer=TRUE,upper=.Machine$integer.max)
 if(!is.null(extraction$simplify))extraction$simplify<-.spec_number(extraction$simplify,"simplify")
 if(!is.logical(extraction$thin)||length(extraction$thin)!=1L||is.na(extraction$thin))stop("invalid thin",call.=FALSE)
 for(k in c("method","provider","provider_version")) {
  v<-extraction[[k]]
  if((k=="method"||!is.null(v))&&(!is.character(v)||length(v)!=1L||is.na(v)||!nzchar(v)))stop("invalid extraction.",k,call.=FALSE)
 }
 extraction<-.spec_lines(extraction)
 if(xor(is.null(extraction$provider),is.null(extraction$provider_version)))stop("provider and provider_version must be supplied together",call.=FALSE)
 args<-list(properties=properties,bins=bins,normalize=normalize,length_range=as.double(length_range),
  curvature_range=as.double(curvature_range),units=units,pixels_per_degree=pixels_per_degree,
  estimator=estimator,curvature_scale=curvature_scale,branch_estimator=branch_estimator,
  branch_window=branch_window,branch_power=branch_power,reference_radius=reference_radius,
  fit_tolerance=fit_tolerance,junction=junction,mat=mat,mapping_radius=mapping_radius,
  raster_line_width=raster_line_width,extraction=extraction)
 structure(list(version=1L,definitions=c(contour=.contour_definition(estimator),
  branch=if(branch_estimator=="mlv2023")"mlv2023-port-0.1.1" else "physical-integrals-1",histogram="1",
  aof=.aof_definition(mat$aof),morphology=.morphology_definition(mat$morphology),
  branch_rule=.branch_rule_definition(mat$branch_rule),
  junction=if(junction$rule=="mlv2023")"mlv2023-port-f12418a" else paste0("exit-neighbourhood-1-",junction$gap_metric),
  tracing=.tracing_definition(mat$branch_rule),
  if(!is.null(extraction$lines))c(lines=extraction$provider_version)),settings=args),class="mlv_feature_spec")
}
# Canonical built-in parameters select the extraction definition.
.spec_lines <- function(extraction) {
 v<-extraction$provider_version;lines<-extraction$lines
 named<-!is.null(v)&&startsWith(v,"fdog-xdog-")
 if(named&&!v%in%.line_definitions)
  stop("unsupported line extraction definition '",v,"'",call.=FALSE)
 if(!is.null(extraction$provider)&&extraction$provider%in%.line_providers&&!is.null(v)&&!named)
  stop("provider '",extraction$provider,"' requires provider_version matching its line parameters",call.=FALSE)
 # Without line settings the field is absent, so specifications built before
 # it existed remain canonical.
 if(is.null(lines)&&!named){extraction$lines<-NULL;return(extraction)}
 if(is.null(lines))lines<-list()
 if(!is.list(lines))stop("extraction.lines must be a named list of line extraction parameters",call.=FALSE)
 extraction$lines<-.line_parameters(lines)
 definition<-.line_definition(extraction$lines)
 if(!is.null(v)&&!identical(v,definition))
  stop("extraction.lines requires provider_version '",definition,"', got '",v,"'",call.=FALSE)
 if(is.null(extraction$provider)&&is.null(v)){
  extraction$provider<-"xdog_lines";extraction$provider_version<-definition
 }else if(is.null(v))stop("extraction.lines requires provider_version '",definition,"'",call.=FALSE)
 extraction
}
#' Validate a feature specification
#' @param spec A feature specification.
#' @return The specification invisibly; malformed or unsupported definitions fail.
#' @export
validate_feature_spec <- function(spec) {
 if(!inherits(spec,"mlv_feature_spec")||!identical(names(spec),c("version","definitions","settings"))||
    !identical(spec$version,1L)||!is.list(spec$settings))stop("invalid or unsupported feature specification",call.=FALSE)
 expected_names<-names(formals(feature_spec))
 if(!identical(names(spec$settings),expected_names))stop("missing, reordered or unknown specification fields",call.=FALSE)
 rebuilt<-do.call(feature_spec,spec$settings)
 if(!identical(spec,rebuilt))stop("noncanonical or unsupported feature specification definitions",call.=FALSE)
 invisible(spec)
}
.spec_semantics <- function(spec) {
 validate_feature_spec(spec)
 spec$settings$pixels_per_degree<-NULL
 spec
}
.spec_scale <- function(spec,pixels_per_degree=NULL) {
 validate_feature_spec(spec)
 if(spec$settings$units=="pixels")return(1)
 ppd<-if(is.null(pixels_per_degree))spec$settings$pixels_per_degree else pixels_per_degree
 1/.spec_number(ppd,"pixels_per_degree",strict=TRUE)
}

# Versioned contour curvature definitions.
.contour_definition <- function(estimator) switch(estimator,mlv2023="mlv2023-port-0.1.1",
 turning="vertex-dual-1",arclength="tent-source-normalized-1")
# Versioned MAT method identifiers. Specifications and MATs built with other
# identifiers (including the 0.2.x Pratt default) fail canonical validation.
.aof_definition <- function(aof) switch(aof,centered="centered-tie-mean-1",mlv2023="mlv2023-computeAOF-lowest-index")
.branch_rule_definition <- function(rule) switch(rule,
 invariant="direction-symmetric-tiemax-graph-defined-3",mlv2023="mlv2023-edgelink-lastwriter-1")
.morphology_definition <- function(morphology) switch(morphology,
 minimal="hmt8-then-sequential-simple-pixel-removal-1",
 hmt8="sequential_hmt8_matlab_skel_order-1",pratt="octave_image_pratt_2.20.0")
.tracing_definition <- function(rule) if(rule=="invariant")"pixel_graph_no_corner_cut" else "kovesi_edgelink_cleanedgelist"
# The centred AOF is equivariant under all eight lattice symmetries only when the
# sample angles 2*pi*k/samples are mapped onto themselves by a quarter turn and
# by the diagonal reflection, i.e. when samples is a multiple of four.
.check_centered_samples <- function(aof,samples) {
 if(identical(aof,"centered")&&is.numeric(samples)&&length(samples)==1L&&is.finite(samples)&&samples%%4!=0)
  stop("aof = 'centered' requires samples to be a multiple of 4 (lattice equivariance); got ",samples,call.=FALSE)
 invisible(TRUE)
}
