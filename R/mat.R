.moving_mean <- function(x,k) {
 n<-length(x);if(!n)return(numeric())
 left<-floor(k/2);right<-k-left-1L
 lo<-pmax(1L,seq_len(n)-left);hi<-pmin(n,seq_len(n)+right)
 cs<-c(0,cumsum(x));(cs[hi+1L]-cs[lo])/(hi-lo+1L)
}
.pad_diff <- function(x) if(length(x)>1L){d<-diff(x);c(d,tail(d,1))} else rep(0,length(x))
#' Average outward flux of a signed-distance/nearest-boundary field
#' @param distance Numeric signed-distance matrix.
#' @param index Nearest-boundary indices in R column-major order.
#' @param samples Number of angular samples.
#' @param epsilon Sampling-circle radius.
#' @return Numeric flux matrix, without mutating the inputs.
#' @export
average_outward_flux <- function(distance,index,samples=60L,epsilon=1) {
 if(!is.matrix(distance)||!is.numeric(distance)||anyNA(distance) ||
    any(dim(distance)<1L)||prod(as.double(dim(distance)))>.Machine$integer.max ||
    !is.matrix(index)||!identical(dim(index),dim(distance)))
  stop("distance and index must be matching nonempty matrices",call.=FALSE)
 if(any(!is.na(index)&(!is.finite(index)|index<1|index>length(distance)|index!=floor(index))))
  stop("invalid nearest-boundary indices",call.=FALSE)
 if(length(samples)!=1L||!is.finite(samples)||samples<4||samples!=floor(samples)||samples>100000L ||
    length(epsilon)!=1L||!is.finite(epsilon)||epsilon<=0)
  stop("invalid samples or epsilon",call.=FALSE)
 storage.mode(index)<-"integer"
 cpp_aof(distance,index,as.integer(samples),epsilon)
}
#' Construct an AOF medial-axis representation
#'
#' Uses the MLV signed-distance construction and a 60-sample average outward
#' flux. The defaults are the corrected methods: a pixel-centred, tie-invariant
#' flux (`aof="centered"`), a topology-preserving thinning that leaves no
#' removable (simple) non-end pixel (`morphology="minimal"`) and order-invariant branch handling
#' (`branch_rule="invariant"`). The published-compatible recipe is
#' `aof="mlv2023", morphology="hmt8", branch_rule="mlv2023"`, with
#' `threshold_angle = asin(0.3*pi/2)*180/pi` (28.11 degrees, flux threshold 0.3)
#' to match the lab's datasets, which `computeAllMATfromVecLD` produced with the
#' omitted-threshold value 0.3 (MLV-M004; the default 28 degrees gives 0.2989).
#' The first two reproduce the published MATLAB MAT fields (distance map, AOF
#' and skeleton) exactly on the Forest_DBW reference, and all three together
#' reproduce the published per-pixel scores. `morphology="pratt"` is the Octave image backend.
#' Empty ink sets have Inf distances and no skeleton; an all-ink image has no
#' valid medial axis under this contour-background definition. Isolated ink
#' pixels retain signed distance -1. The `"mlv2023"` flux uses the lowest R
#' linear index among equidistant nearest boundary pixels, as MATLAB bwdist does.
#' @param mask Binary matrix, TRUE for ink; or a drawing.
#' @param threshold_angle Object angle in degrees, in [0,90].
#' @param samples,epsilon AOF sampling parameters. Angles use alpha=(2*pi)/samples followed by alpha*k;
#'   the order matters when samples are mapped to discrete grid cells. With
#'   `aof="centered"`, samples must be a multiple of 4 (the sample directions
#'   must be closed under quarter turns and the diagonal reflection).
#' @param min_component Minimum component size; default scales with image size.
#' @param min_branch_points Spur length threshold, in pixels (default
#'   `branch_rule`) or Euclidean length (`"mlv2023"`). Under `"invariant"` only
#'   spurs (paths from an end pixel to a junction) shorter than this are pruned,
#'   shortest first with tied spurs kept together, and never a component's last
#'   path; see [trace_skeleton()].
#' @param trace_branches Trace branches now; FALSE defers this work until scoring.
#' @param spec Optional feature specification for raster and MAT settings.
#' @param aof Flux sampling. `"centered"` (default) centres the sampling circle
#'   on the pixel, reads the nearest pixel cell (both cells for a point on a cell
#'   boundary) and averages unit directions over all equidistant nearest boundary
#'   pixels; it is independent of tie order and, for samples divisible by 4, the
#'   flux field of a given ink raster is equivariant under the eight lattice
#'   rotations/reflections. Pixels whose sampling circle leaves the image are
#'   set to 0 symmetrically on all four borders. `"mlv2023"` reproduces MATLAB computeAOF.m exactly,
#'   including lowest-index tie labels: its circle is centred on the pixel corner
#'   (i+0.5,j+0.5) while labels are read from `floor()` cells, which misplaces
#'   medial ridges by half a pixel and breaks that equivariance.
#' @param morphology Thinning used for both ink and flux skeletons; see
#'   [skeletonize()]. `"minimal"` (default; no simple pixel with two or more
#'   neighbours remains, 2x2 blocks only where all four pixels are non-simple),
#'   `"hmt8"` (MATLAB-exact on the reference; may keep simple non-end pixels and
#'   removable 2x2 blocks) or `"pratt"`.
#' @param branch_rule Branch tracing, scoring direction and contour mapping.
#'   `"invariant"` (default): pixel-graph tracing without repeated points and
#'   order-free spur pruning that keeps every skeleton component of two or more pixels (see
#'   [trace_skeleton()]),
#'   branch scores averaged over both traversal directions, tangent samples
#'   carrying the mean score of their segment's two points, the maximum over
#'   exact nearest-sample ties and over branches, and valid zeros retained;
#'   an isolated closed loop (no end or junction pixel) is scored with circular
#'   windows for local estimators (it has no start) and, for the mlv2023 mirror
#'   fit (a global greedy segmentation), as the mean over every start pixel and
#'   both directions;
#'   traced branches and per-pixel scores are unchanged by branch reversal,
#'   branch or contour renumbering, the traced start pixel (including that of
#'   an isolated closed loop), and lattice rotations/reflections of the MAT
#'   fields. The whole pipeline from a drawing
#'   is not guaranteed lattice-equivariant: thinning is sequential
#'   (direction-ordered) and the rasterizer resolves exact half-pixel ties by
#'   slope sign.
#'   `"mlv2023"`: the published chain (Kovesi edgelink/cleanedgelist tracing with
#'   duplicated merge points, forward scores, MLV fitLineSegments for mirror,
#'   last-writer skeleton rating image, MATLAB tangent enumeration order with
#'   lowest-index ties, zeros treated as missing). Together with
#'   `aof="mlv2023"` and `morphology="hmt8"` (the published-compatible recipe)
#'   this reproduces published MLV per-pixel scores.
#' @return An `mlv_mat` with distance_map, nearest_boundary, AOF, skeleton,
#'   branches, mask and method parameters.
#' @export
compute_mat <- function(mask,threshold_angle=28,samples=60L,epsilon=1,
                        min_component=NULL,min_branch_points=10L,trace_branches=TRUE,spec=NULL,
                        aof=c("centered","mlv2023"),morphology=c("minimal","hmt8","pratt"),
                        branch_rule=c("invariant","mlv2023")) {
 drawing<-if(inherits(mask,"mlv_drawing"))mask else NULL
 method_given<-!missing(aof)||!missing(morphology)||!missing(branch_rule)
 aof<-match.arg(aof);morphology<-match.arg(morphology);branch_rule<-match.arg(branch_rule)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(threshold_angle)||!missing(samples)||!missing(epsilon)||!missing(min_component)||!missing(min_branch_points)||method_given)
   stop("use spec or individual MAT settings, not both",call.=FALSE)
  s<-spec$settings$mat;threshold_angle<-s$threshold_angle;samples<-s$samples;epsilon<-s$epsilon
  min_component<-s$min_component;min_branch_points<-s$min_branch_points
  aof<-s$aof;morphology<-s$morphology;branch_rule<-s$branch_rule
 }
 .check_centered_samples(aof,samples)
 definitions<-list(aof_definition=.aof_definition(aof),morphology_definition=.morphology_definition(morphology),
  branch_rule=branch_rule,branch_rule_definition=.branch_rule_definition(branch_rule),
  tracing=if(branch_rule=="invariant")"pixel_graph_no_corner_cut" else "kovesi_edgelink_cleanedgelist")
 attach_spec<-function(out) {
  if(!is.null(spec)) {
   out$feature_spec<-spec
   if(!is.null(drawing))out$source_geometry<-.geometry_snapshot(drawing)
  }
  out
 }
 if(!is.null(drawing))mask<-rasterize(mask,line_width=if(is.null(spec))1 else spec$settings$raster_line_width)
 mask<-.check_mask(mask); h<-nrow(mask);w<-ncol(mask)
 if(!is.logical(trace_branches)||length(trace_branches)!=1L||is.na(trace_branches))stop("trace_branches must be TRUE or FALSE",call.=FALSE)
 if(length(threshold_angle)!=1L||!is.finite(threshold_angle)||threshold_angle<0||threshold_angle>90)
  stop("threshold_angle must be between 0 and 90 degrees",call.=FALSE)
 if(is.null(min_component))min_component<-max(floor(.025*max(h,w)),1L)
 if(!any(mask)||all(mask)) {
  D<-matrix(if(all(mask)) -Inf else Inf,h,w)
  idx<-matrix(NA_integer_,h,w)
  flux<-average_outward_flux(D,idx,samples,epsilon)
  sk<-skeletonize(matrix(FALSE,h,w),min_component,method=morphology)
  if(length(min_branch_points)!=1L||!is.finite(min_branch_points)||min_branch_points<2||min_branch_points!=floor(min_branch_points))
   stop("min_branch_points must be at least 2",call.=FALSE)
  return(attach_spec(structure(list(distance_map=D,nearest_boundary=idx,AOF=flux,skeleton=sk,
   mask=mask,branches=list(),parameters=list(threshold_angle=threshold_angle,
   threshold=2/pi*sin(threshold_angle*pi/180),samples=samples,epsilon=epsilon,
   min_component=min_component,min_branch_points=min_branch_points,aof=aof,
   morphology=morphology,aof_definition=definitions$aof_definition,morphology_definition=definitions$morphology_definition,
   branch_rule=branch_rule,branch_rule_definition=definitions$branch_rule_definition,
   skeleton_ink_excluded=0L,tracing=definitions$tracing)),class="mlv_mat")))
 }
 thin_ink<-skeletonize(mask,method=morphology)
 white<-!thin_ink
 outer<-matrix(FALSE,h,w);ink_neighbor<-matrix(FALSE,h,w)
 if(h>=3L&&w>=3L)for(dc in -1:1)for(dr in -1:1)if(dr!=0L||dc!=0L) {
  rows<-2:(h-1L);cols<-2:(w-1L)
  outer[rows,cols]<-outer[rows,cols] | (thin_ink[rows,cols]&white[rows+dr,cols+dc])
  ink_neighbor[rows,cols]<-ink_neighbor[rows,cols] | thin_ink[rows+dr,cols+dc]
 }
 outer<-outer&ink_neighbor
 d1<-distance_transform(thin_ink)
 d2<-distance_transform(white|outer)
 D<-d1$distance-d2$distance
 index<-d1$index
 on_ink<-d1$distance==0
 index[on_ink]<-d2$index[on_ink]
 index[outer]<-which(outer)
 flux<-if(aof=="mlv2023")average_outward_flux(D,index,samples,epsilon) else {
  if(length(samples)!=1L||!is.finite(samples)||samples<4||samples!=floor(samples)||samples>100000L ||
     length(epsilon)!=1L||!is.finite(epsilon)||epsilon<=0)stop("invalid samples or epsilon",call.=FALSE)
  cpp_aof_centered(D,thin_ink,outer,as.integer(samples),epsilon)
 }
 threshold<-2/pi*sin(threshold_angle*pi/180)
 sk<-skeletonize(flux>threshold*samples,min_component,method=morphology)
 out<-structure(list(distance_map=D,nearest_boundary=index,AOF=flux,skeleton=sk,
   mask=mask,branches=list(),parameters=list(threshold_angle=threshold_angle,
   threshold=threshold,samples=samples,epsilon=epsilon,min_component=min_component,
   min_branch_points=min_branch_points,aof=aof,morphology=morphology,aof_definition=definitions$aof_definition,morphology_definition=definitions$morphology_definition,
   branch_rule=branch_rule,branch_rule_definition=definitions$branch_rule_definition,
   skeleton_ink_excluded=sum(sk&!(D>0)),tracing=definitions$tracing)),class="mlv_mat")
 if(trace_branches)out$branches<-trace_skeleton(out,min_branch_points)
 attach_spec(out)
}
#' Trace ordered medial-axis branches
#'
#' Skeleton pixels with a nonpositive signed distance lie on ink and are not
#' medial points; they are excluded before tracing on every rule, and their
#' count is returned as attribute `skeleton_ink_excluded`.
#' `rule="invariant"` decomposes the pixel graph (diagonal steps only where no
#' orthogonal connection exists) into maximal paths between end/junction
#' pixels, then prunes spurs and merges to a fixed point. A *spur* is a path
#' from an end pixel to a junction; it is short when it has fewer than
#' `min_points` pixels (junction included). Each round removes the shortest
#' short spurs (all tied ones together) at every junction that keeps at least
#' two paths afterwards; a junction where removing its tied shortest spurs would
#' leave fewer than two paths is protected and keeps all its paths. Paths are
#' then joined through every pixel where exactly two remain. Without ties this
#' is sequential shortest-first spur pruning; with ties (the four equal arms of
#' a square) every tied spur is kept, so the result does not depend on
#' traversal order, start pixel or lattice orientation. Only spurs are pruned:
#' a path between two junctions, or a component that is a single path, is kept
#' whatever its length, so no skeleton component loses its last path. A closed
#' path whose start pixel joins no other path is an isolated loop. Junction
#' pixels are shared by the paths that meet there but never repeated within a
#' path.
#' `rule="mlv2023"` is the published MLV tracing: Kovesi `edgelink` with
#' `cleanedgelist(minlength=min_points)` (Euclidean length), ported line by
#' line, including its raster-scan traversal direction and the duplicated point
#' at every merged junction.
#' @param mat An mlv_mat object.
#' @param min_points Minimum spur size: pixels for `"invariant"` (spurs only),
#'   Euclidean length for `"mlv2023"` (cleanedgelist spurs and isolated edges).
#' @param rule Branch rule; defaults to the rule recorded in `mat`.
#' @return A list of branch data frames with row, column, x, y, radius and AOF.
#' @export
trace_skeleton <- function(mat,min_points=10L,rule=NULL) {
 if(!inherits(mat,"mlv_mat"))stop("mat must be an mlv_mat",call.=FALSE)
 if(length(min_points)!=1L||!is.finite(min_points)||min_points<2||min_points!=floor(min_points))
  stop("min_points must be at least 2",call.=FALSE)
 rule<-.branch_rule(mat,rule)
 sk<-.check_mask(mat$skeleton)
 D<-mat$distance_map
 medial<-if(is.null(D))sk else sk&!is.na(D)&D>0   # ink is never medial (MLV-R012)
 excluded<-sum(sk&!medial)
 paths<-if(rule=="invariant") {
  lapply(cpp_trace_graph_pruned(medial,as.integer(min_points)),function(index)arrayInd(index,dim(sk)))
 } else .trace_edgelink(medial,min_points)
 out<-lapply(paths,function(rc) {
  index<-rc[,1]+nrow(sk)*(rc[,2]-1)
  data.frame(row=rc[,1],column=rc[,2],x=rc[,2],y=rc[,1],
             radius=if(is.null(D))NA_real_ else D[index],AOF=if(is.null(mat$AOF))NA_real_ else mat$AOF[index])
 })
 attr(out,"skeleton_ink_excluded")<-excluded
 out
}
.branch_rule <- function(mat,rule=NULL) {
 if(is.null(rule))rule<-mat$parameters$branch_rule
 if(is.null(rule))rule<-"invariant"
 match.arg(rule,c("invariant","mlv2023"))
}
#' Greedy line fitting for ordered points
#' @param xy Finite two-column point matrix.
#' @param tolerance Maximum mean interior perpendicular distance.
#' @return List with endpoint-index spans and per-point straightness scores.
#' @export
fit_line_segments <- function(xy,tolerance=1) {
 if(!is.matrix(xy)||!is.numeric(xy)||ncol(xy)!=2L||any(!is.finite(xy)))
  stop("xy must be a finite two-column matrix",call.=FALSE)
 if(length(tolerance)!=1L||!is.finite(tolerance)||tolerance<0)stop("invalid tolerance",call.=FALSE)
 cpp_fit_lines(xy,tolerance)
}
#' Score one skeletal branch
#'
#' The taper method uses explicit moving means of width three for both
#' derivative stages; this replaces the upstream unspecified adaptive smoothing.
#' Final scores use a width-three moving mean and power ten. Separation radii
#' below one pixel are undefined and return NA.
#' @param branch A data frame/list with x, y and radius vectors.
#' @param property Parallelism, separation, mirror or taper.
#' @param window Local derivative integration half-window.
#' @param power Final score exponent.
#' @param estimator Published sample-based scores or physical arc-length scores.
#' @param reference_radius Positive reference radius for new separation scores.
#' @param fit_tolerance Positive geometric residual scale for new mirror scores.
#' @param units Spatial units for the new estimator's window and references.
#' @param pixels_per_degree Required calibration for visual-degree units.
#' @details The arclength estimator integrates the piecewise-linear axis and
#' radius exactly over clipped physical windows. Parallelism is axis arc length
#' divided by the joint axis/radius arc length. Mirror is chord/arc straightness
#' divided by one plus mean squared perpendicular error / fit_tolerance squared.
#' Taper is consistency with the best affine radius, using a dimensionless RMS
#' residual divided by window length. Separation is mean radius divided by
#' (mean radius + reference_radius). The exponent applies afterward. These are
#' versioned geometric descriptors, not claims of perceptual equivalence. A
#' zero-length window is undefined (NA); inconsistent radii at duplicate points
#' fail. The published estimator ignores none of these options: nondefault new
#' unit/reference settings require explicit arclength selection.
#' @return One numeric score per branch point.
#' @export
score_branch <- function(branch,property=c("parallelism","separation","mirror","taper"),
                         window=5L,power=10,estimator=c("mlv2023","arclength"),
                         reference_radius=1,fit_tolerance=1,units=c("pixels","visual_degrees"),
                         pixels_per_degree=NULL) {
 property<-match.arg(property)
 estimator<-match.arg(estimator);units<-match.arg(units)
 if(estimator=="arclength")return(.score_branch_arclength(branch,property,window,power,
  reference_radius,fit_tolerance,units,pixels_per_degree))
 if(units!="pixels"||!is.null(pixels_per_degree)||!identical(as.double(reference_radius),1)||
    !identical(as.double(fit_tolerance),1))stop("new units/references require estimator='arclength'",call.=FALSE)
 X<-branch$x;Y<-branch$y;R<-branch$radius;n<-length(X)
 if(length(Y)!=n||length(R)!=n||any(!is.finite(c(X,Y,R)))||any(R<0))
  stop("branch coordinates/radii must be finite, nonnegative and equally sized",call.=FALSE)
 if(length(window)!=1L||!is.finite(window)||window<1||window!=floor(window) ||
    length(power)!=1L||!is.finite(power)||power<=0)
  stop("invalid window or power",call.=FALSE)
 if(!n)return(numeric())
 R<-.moving_mean(R,2L);dx<-.pad_diff(X);dy<-.pad_diff(Y);dr<-.pad_diff(R)
 out<-numeric(n)
 if(property=="separation") {
  out<-ifelse(R>=1,1-1/R,NA_real_)
 } else if(property=="mirror") {
  if(n>3L)out<-fit_line_segments(cbind(X,Y))$scores
 } else {
  derivative<-if(property=="taper").moving_mean(.pad_diff(.moving_mean(dr,3L)),3L) else dr
  axis<-sqrt(dx^2+dy^2);arc<-sqrt(dx^2+dy^2+derivative^2)
  if(n>=3L) {
   ii<-2:(n-1L);k<-pmin(ii-1L,n-ii,window);lo<-ii-k;hi<-ii+k
   a<-c(0,cumsum(axis));b<-c(0,cumsum(arc))
   den<-b[hi+1L]-b[lo];out[ii]<-ifelse(den>0,(a[hi+1L]-a[lo])/den,0)
  }
 }
 # Missingness is not allowed to contaminate every later cumulative window.
 if(anyNA(out)) {
  out<-vapply(seq_len(n),function(i){v<-out[max(1,i-1):min(n,i+1)];if(anyNA(v))NA_real_ else mean(v)},double(1))
 } else out<-.moving_mean(out,3L)
 out^power
}
#' Compute medial-axis symmetry scores
#' @param mat An mlv_mat object.
#' @param properties Requested branch score names.
#' @param window,power Branch score parameters.
#' @param spec Optional feature specification, exclusive of individual settings.
#' @param pixels_per_degree Optional input calibration override.
#' @return An mlv_mat with branch scores and score images; unmapped pixels are NA.
#' @export
compute_mat_properties <- function(mat,properties=c("parallelism","separation","mirror"),
                                   window=5L,power=10,spec=mat$feature_spec,pixels_per_degree=NULL) {
 if(!inherits(mat,"mlv_mat"))stop("mat must be an mlv_mat",call.=FALSE)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(properties)||!missing(window)||!missing(power))stop("use spec or individual branch settings, not both",call.=FALSE)
  s<-spec$settings
  properties<-intersect(s$properties,c("parallelism","separation","mirror","taper"))
  window<-s$branch_window;power<-s$branch_power
  actual<-mat$parameters[names(s$mat)];expected<-s$mat
  if(is.null(expected$min_component))expected$min_component<-max(floor(.025*max(dim(mat$mask))),1L)
  if(!isTRUE(all.equal(actual,expected,check.attributes=FALSE)) ||
     !identical(mat$parameters$aof_definition,spec$definitions[["aof"]]) ||
     !identical(mat$parameters$morphology_definition,spec$definitions[["morphology"]]) ||
     !identical(mat$parameters$branch_rule_definition,spec$definitions[["branch_rule"]]))
   stop("MAT raster specification mismatch; recompute compute_mat()",call.=FALSE)
  if(!is.null(pixels_per_degree))spec$settings$pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
 }
 if(any(!properties %in% c("parallelism","separation","mirror","taper")))
  stop("unknown MAT property",call.=FALSE)
 rule<-.branch_rule(mat)
 if(!length(mat$branches)&&any(mat$skeleton))mat$branches<-trace_skeleton(mat,mat$parameters$min_branch_points)
 estimator<-if(is.null(spec))"mlv2023" else s$branch_estimator
 loop_window<-window*(if(!is.null(spec)&&identical(s$units,"visual_degrees")&&!is.null(spec$settings$pixels_per_degree))spec$settings$pixels_per_degree else 1)
 one<-function(b,prop)if(is.null(spec))score_branch(b,prop,window,power) else score_branch(b,prop,window,power,
    estimator=s$branch_estimator,reference_radius=s$reference_radius,fit_tolerance=s$fit_tolerance,
    units=s$units,pixels_per_degree=if(s$branch_estimator=="mlv2023")NULL else spec$settings$pixels_per_degree)
 isolated<-.isolated_loops(mat$branches)
 mat$score_images<-list()
 for(prop in properties) {
  image<-matrix(NA_real_,nrow(mat$skeleton),ncol(mat$skeleton))
  for(i in seq_along(mat$branches)) {
   b<-mat$branches[[i]]
   scores<-if(rule=="mlv2023") {
    if(prop=="mirror"&&estimator=="mlv2023").mirror_mlv(b,power) else one(b,prop)
   } else {
    # Direction-symmetric scores: the mean of the forward scores and of the
    # scores of the reversed branch, mapped back to the original order. An
    # isolated closed loop has no start: local estimators use circular windows
    # (both directions); the mlv2023 mirror fit is averaged over all starts.
    n<-nrow(b)
    if(isolated[i]).score_loop(b,function(z)one(z,prop),prop,estimator,power,loop_window) else {
     f<-one(b,prop)
     if(n>1L)(f+rev(one(b[n:1,,drop=FALSE],prop)))/2 else f
    }
   }
   mat$branches[[i]][[prop]]<-scores
   idx<-cbind(b$row,b$column)
   if(rule=="mlv2023") image[idx]<-scores else {
    # Maximum over all occurrences of a pixel, within the branch (a junction
    # pixel closing a loop occurs twice) and across branches.
    key<-b$row+nrow(image)*(b$column-1)
    best<-vapply(split(scores,key),function(v)if(all(is.na(v)))NA_real_ else max(v,na.rm=TRUE),double(1))
    at<-as.integer(names(best));old<-image[at]
    image[at]<-ifelse(is.na(old),best,pmax(old,best,na.rm=TRUE))
   }
  }
  mat$score_images[[prop]]<-image
 }
 mat$parameters$window<-window;mat$parameters$power<-power
 mat$parameters$taper_smoothing<-"moving_mean_3"
 if(!is.null(spec)) {
  mat$feature_spec<-spec
  mat$score_spec<-.mat_semantics(spec)
  if(s$branch_estimator!="mlv2023")mat$parameters$taper_smoothing<-"exact_affine_radius_residual"
 }
 mat
}
.tangent_points <- function(branch) {
 n<-nrow(branch)
 if(n<2L)return(list(points=matrix(numeric(),0,2),source=integer()))
 X<-branch$row;Y<-branch$column;R<-branch$radius
 x1<-head(X,-1);y1<-head(Y,-1);x2<-tail(X,-1);y2<-tail(Y,-1)
 base1<-head(R,-1);base2<-tail(R,-1)
 d<-sqrt((x1-x2)^2+(y1-y2)^2)
 mx<-(x1+x2)/2;my<-(y1+y2)/2;vx<-x1-mx;vy<-y1-my
 points<-list();sources<-integer()
 for(delta1 in c(-.1,0,.1))for(delta2 in c(-.1,0,.1)) {
  r1<-base1+delta1;r2<-base2+delta2
  # Coincident consecutive points (d == 0) define no tangent direction; they
  # contribute no tangent samples instead of producing NaN subscripts (MLV-R011).
  ca<-ifelse(d>0,(r2-r1)/d,0)
  bad<-abs(ca)>1;ca[bad]<-ca[bad]/(abs(ca[bad])+.00001)
  sa<-sqrt(pmax(0,1-ca^2));mr<-(r1+r2)/2
  for(sign in c(1,-1)) {
   xx<-ca*vx+sign*sa*vy;yy<--sign*sa*vx+ca*vy
   norm<-sqrt(xx^2+yy^2);valid<-d>0&is.finite(norm)&norm>0&mr>=0
   points[[length(points)+1L]]<-cbind(mx[valid]+xx[valid]*mr[valid]/norm[valid],
                                    my[valid]+yy[valid]*mr[valid]/norm[valid])
   sources<-c(sources,which(valid))
  }
 }
 list(points=do.call(rbind,points),source=sources)
}
#' Map branch ratings back to drawing contours
#'
#' Uses exact nearest tangent samples per branch, then the maximum across
#' contributing branches. An explicit mapping mask distinguishes valid zero
#' ratings from missing samples. Geometry is reused across requested properties.
#' @param x A drawing.
#' @param mat An mlv_mat with computed branch properties.
#' @param radius Maximum tangent mapping distance in pixels.
#' @param spec Optional feature specification retained by the scored MAT.
#' @param pixels_per_degree Optional input calibration override.
#' @return A drawing with per-contour MAT scores, score images and mapping masks.
#' @export
map_mat_properties <- function(x,mat,radius=4,spec=mat$feature_spec,pixels_per_degree=NULL) {
 .check_drawing(x)
 if(!inherits(mat,"mlv_mat")||is.null(mat$score_images))stop("compute MAT properties first",call.=FALSE)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(radius))stop("use spec or radius, not both",call.=FALSE)
  radius<-spec$settings$mapping_radius
  if(!is.null(pixels_per_degree))spec$settings$pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
  if(is.null(mat$score_spec)||!identical(mat$score_spec,.mat_semantics(spec)))
   stop("MAT score specification mismatch; recompute compute_mat_properties()",call.=FALSE)
  if(!is.null(mat$source_geometry)&&!identical(mat$source_geometry,.geometry_snapshot(x)))
   stop("MAT source geometry mismatch; recompute analyse_mat()",call.=FALSE)
 }
 if(!identical(dim(mat$mask),c(x$height,x$width)))stop("drawing and MAT image dimensions differ",call.=FALSE)
 if(length(radius)!=1L||!is.finite(radius)||radius<=0)stop("radius must be positive",call.=FALSE)
 ink<-which(rasterize(x));rc<-arrayInd(ink,dim(mat$mask))
 props<-names(mat$score_images);images<-lapply(props,function(p)rep(NA_real_,length(ink)));names(images)<-props
 rule<-.branch_rule(mat);tol<-1e-9
 for(branch in mat$branches) {
  fp<-if(rule=="mlv2023").tangent_points_mlv(branch) else .tangent_points(branch)
  if(!nrow(fp$points)||!length(ink))next
  # A point outside this expanded bounding box cannot be within radius of
  # any tangent sample. Keep equality here; the final predicate remains strict.
  near<-which(rc[,1]>=min(fp$points[,1])-radius & rc[,1]<=max(fp$points[,1])+radius &
              rc[,2]>=min(fp$points[,2])-radius & rc[,2]<=max(fp$points[,2])+radius)
  if(!length(near))next
  if(rule=="mlv2023") {
   # mapMATtoContour: knnsearch (lowest index among exact ties) and scores read
   # from the last-writer skeleton rating image at the segment's first point.
   nn<-cpp_nearest(fp$points,rc[near,,drop=FALSE])
   within<-which(nn$distance<radius);if(!length(within))next
   ok<-near[within];src<-fp$source[nn$index[within]];at_sk<-cbind(branch$row[src],branch$column[src])
   for(prop in props) {
    value<-mat$score_images[[prop]][at_sk];valid<-is.finite(value);at<-ok[valid];value<-value[valid]
    old<-images[[prop]][at];images[[prop]][at]<-ifelse(is.na(old),value,pmax(old,value))
   }
  } else {
   # A tangent sample of segment (i,i+1) carries the mean of the branch's own
   # scores at i and i+1; exact distance ties take the maximum; samples beyond
   # radius (less a rounding tolerance) are ignored.
   seg<-vapply(props,function(prop){v<-branch[[prop]];if(is.null(v))v<-rep(NA_real_,nrow(branch))
    (head(v,-1)+tail(v,-1))/2},double(nrow(branch)-1L))
   seg<-matrix(seg,ncol=length(props))
   nn<-cpp_nearest_max(fp$points,rc[near,,drop=FALSE],seg[fp$source,,drop=FALSE],tol)
   within<-which(nn$distance<radius-tol);if(!length(within))next
   ok<-near[within]
   for(q in seq_along(props)) {
    value<-nn$value[within,q];valid<-is.finite(value);at<-ok[valid];value<-value[valid]
    old<-images[[props[q]]][at];images[[props[q]]][at]<-ifelse(is.na(old),value,pmax(old,value))
   }
  }
 }
 # MATpropertiesToContours treats every zero as missing (MLV-M011); the
 # published-compatible rule reproduces that, the invariant rule keeps valid zeros.
 if(rule=="mlv2023")for(prop in props)images[[prop]][images[[prop]]==0]<-NA_real_
 x$mat_properties<-list();x$mat_images<-list();x$mat_mapping_masks<-list()
 contour_indices<-lapply(seq_along(x$ids),function(i).raster_indices(x,i))
 for(prop in props) {
  image<-matrix(NA_real_,x$height,x$width);image[ink]<-images[[prop]]
  x$mat_images[[prop]]<-image
  x$mat_mapping_masks[[prop]]<-is.finite(image)
  x$mat_properties[[prop]]<-lapply(contour_indices,function(index)image[index])
 }
 x$mat_parameters<-mat$parameters
 if(!is.null(spec)) {
  x$feature_spec<-spec;x$mat_feature_spec<-.mat_semantics(spec)
  x$mat_geometry<-.geometry_snapshot(x)
 }
 x
}
#' Complete vector drawing to medial-axis feature workflow
#' @param x A drawing.
#' @param properties MAT property names.
#' @param ... Arguments to compute_mat.
#' @param spec Optional feature specification, exclusive of properties and dots.
#' @param pixels_per_degree Optional input calibration override.
#' @return List containing the enriched drawing and medial-axis representation.
#' @export
analyse_mat <- function(x,properties=c("parallelism","separation","mirror"),...,
                        spec=x$feature_spec,pixels_per_degree=NULL) {
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(properties)||length(list(...)))stop("use spec or individual MAT settings, not both",call.=FALSE)
  if(!is.null(pixels_per_degree))spec$settings$pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
  mat<-compute_mat(x,spec=spec)
  mat<-compute_mat_properties(mat,spec=spec)
  return(list(drawing=map_mat_properties(x,mat,spec=spec),mat=mat))
 }
 mat<-compute_mat(x,...)
 mat<-compute_mat_properties(mat,properties)
 list(drawing=map_mat_properties(x,mat),mat=mat)
}
#' @export
print.mlv_mat <- function(x,...) {
 cat("<mlv_mat>",nrow(x$skeleton),"x",ncol(x$skeleton),",",
     sum(x$skeleton),"skeleton pixels,",length(x$branches),"branches\n")
 invisible(x)
}

# An isolated closed loop is a closed path (last point repeating the first)
# whose start pixel belongs to no other branch: no end or junction anchors its
# points, so it has no start. A closed path starting on a pixel shared with
# another branch (a junction) is anchored there (a lollipop). For tracer output
# this is the pruned-graph degree; it does not depend on the raw skeleton, whose
# pruned spurs still touch the former junction pixel.
.isolated_loops <- function(branches) {
 if(!length(branches))return(logical())
 key<-function(r,c)paste(r,c)
 count<-table(unlist(lapply(branches,function(b)unique(key(b$row,b$column)))))
 vapply(branches,function(b){n<-nrow(b)
  n>=4L&&b$row[1]==b$row[n]&&b$column[1]==b$column[n]&&count[[key(b$row[1],b$column[1])]]==1L},logical(1))
}
# Whether branch b (an element of branches) is an isolated closed loop.
.is_isolated_loop <- function(b,branches) {
 n<-nrow(b)
 n>=4L&&b$row[1]==b$row[n]&&b$column[1]==b$column[n]&&
  sum(vapply(branches,function(z)any(z$row==b$row[1]&z$column==b$column[1]),logical(1)))==1L
}
# Scores of an isolated closed loop (points 1..m, point m repeating point 1). A
# loop has no ends, so local estimators use circular windows: each point is
# scored on the periodic sequence, padded with P = window + 16 points on each
# side so that its window never reaches a padding end, in both directions (mean).
# This is a definition, not the average over cut points (which would introduce
# artificial ends). Points are scored in blocks of at most 64 so that running
# sums stay short and rounding does not grow with loop length. The mlv2023
# mirror fit is a global greedy segmentation without a local window; it is the
# mean over every start point and both directions (compiled).
.score_loop <- function(b,score,prop,estimator,power,window) {
 u<-b[-nrow(b),,drop=FALSE];m<-nrow(u)
 if(prop=="mirror"&&estimator=="mlv2023") {
  s<-cpp_mirror_loop(cbind(as.double(u$x),as.double(u$y)),power,1)
  return(c(s,s[1]))
 }
 P<-as.integer(ceiling(window))+16L;B<-64L;s<-numeric(m)
 for(start in seq(1L,m,by=B)) {
  end<-min(m,start+B-1L);pos<-(start-P):(end+P);idx<-((pos-1L)%%m)+1L
  mid<-P+seq_len(end-start+1L)
  f<-score(u[idx,,drop=FALSE])[mid]
  r<-rev(score(u[rev(idx),,drop=FALSE]))[mid]
  s[start:end]<-(f+r)/2
 }
 c(s,s[1])
}
