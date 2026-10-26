.cross2 <- function(a,b)a[1]*b[2]-a[2]*b[1]
.graph_components <- function(n,from,to) {
 parent<-seq_len(n)
 root<-function(i){while(parent[i]!=i)i<-parent[i];i}
 for(i in seq_along(from)) {
  a<-root(from[i]);b<-root(to[i]);parent[max(a,b)]<-min(a,b)
 }
 if(n)match(vapply(seq_len(n),root,integer(1)),unique(vapply(seq_len(n),root,integer(1)))) else integer()
}
.graph_arrangement <- function(z,owner,segment,closed,tolerance) {
 n<-nrow(z);cuts<-rep(list(c(0,1)),n);hits<-matrix(FALSE,n,n)
 self_crossing<-rep(FALSE,length(closed));sizes<-tabulate(owner,nbins=length(closed))
 numeric_tol<-64*.Machine$double.eps*max(1,abs(z))
 lengths<-sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))
 candidate_count<-2*n
 if(n>1L)for(i in seq_len(n-1L))for(j in seq.int(i+1L,n)) {
  a<-z[i,1:2];u<-z[i,3:4]-a;b<-z[j,1:2];v<-z[j,3:4]-b
  if(any(pmax(pmin(a,a+u),pmin(b,b+v))>pmin(pmax(a,a+u),pmax(b,b+v))+numeric_tol))next
  denominator<-.cross2(u,v);w<-b-a
  if(abs(denominator)<=64*.Machine$double.eps*lengths[i]*lengths[j]) {
   if(abs(.cross2(w,u))/lengths[i]>numeric_tol)next
   interval<-sort(c(sum(w*u),sum((b+v-a)*u))/sum(u*u))
   low<-max(0,interval[1]);high<-min(1,interval[2])
   if(high<low-numeric_tol/lengths[i])next
   if((high-low)*lengths[i]>numeric_tol)
    stop("collinear overlapping spans are unsupported (source segments ",i," and ",j,")",call.=FALSE)
   t<-pmax(0,pmin(1,(low+high)/2));q<-a+t*u;s<-sum((q-b)*v)/sum(v*v)
  } else {
   t<-.cross2(w,v)/denominator;s<-.cross2(w,u)/denominator
   if(t< -numeric_tol/lengths[i]||t>1+numeric_tol/lengths[i]||
      s< -numeric_tol/lengths[j]||s>1+numeric_tol/lengths[j])next
  }
  cuts[[i]]<-c(cuts[[i]],pmax(0,pmin(1,t)));cuts[[j]]<-c(cuts[[j]],pmax(0,pmin(1,s)))
  hits[i,j]<-hits[j,i]<-TRUE;candidate_count<-candidate_count+2L
  if(candidate_count>20000)stop("graph exceeds 20000 candidate split points",call.=FALSE)
  if(owner[i]==owner[j]) {
   adjacent<-abs(segment[i]-segment[j])==1L ||
    (closed[owner[i]]&&abs(segment[i]-segment[j])==sizes[owner[i]]-1L)
   if(!adjacent)self_crossing[owner[i]]<-TRUE
  }
 }
 cuts<-lapply(cuts,function(t)sort(unique(t)))
 points<-if(n)do.call(rbind,lapply(seq_len(n),function(i)
  cbind(z[i,1]+cuts[[i]]*(z[i,3]-z[i,1]),z[i,2]+cuts[[i]]*(z[i,4]-z[i,2])))) else matrix(numeric(),0,2)
 cluster<-if(nrow(points))cpp_cluster_points(points,tolerance) else integer()
 cluster<-match(cluster,unique(cluster));nn<-length(unique(cluster))
 centers<-if(nn)do.call(rbind,lapply(seq_len(nn),function(i)colMeans(points[cluster==i,,drop=FALSE]))) else matrix(numeric(),0,2)
 offset<-c(0L,cumsum(lengths(cuts)));from<-to<-which_segment<-integer();t0<-t1<-numeric()
 for(i in seq_len(n)) {
  at<-offset[i]+seq_along(cuts[[i]])
  from<-c(from,head(cluster[at],-1));to<-c(to,tail(cluster[at],-1))
  which_segment<-c(which_segment,rep(i,length(at)-1L))
  t0<-c(t0,head(cuts[[i]],-1));t1<-c(t1,tail(cuts[[i]],-1))
 }
 collapsed<-sum(from==to);keep<-from!=to
 edge<-data.frame(from=from[keep],to=to[keep],source_segment=which_segment[keep],t0=t0[keep],t1=t1[keep])
 edge$length<-(edge$t1-edge$t0)*lengths[edge$source_segment]
 component<-.graph_components(nn,edge$from,edge$to)
 nodes<-data.frame(node_id=seq_len(nn),x=centers[,1],y=centers[,2],
                   degree=tabulate(c(edge$from,edge$to),nbins=nn),component=component)
 list(nodes=nodes,edges=edge,hits=hits,self_crossing=self_crossing,collapsed=collapsed,
      first_node=if(n)cluster[offset[-length(offset)]+1L] else integer(),
      last_node=if(n)cluster[offset[-1L]] else integer(),numeric_tolerance=numeric_tol)
}
.graph_continuation <- function(z,owner,closed,arr,gap_max,angle_tolerance,tolerance) {
 endpoints<-list()
 for(i in seq_along(closed))if(!closed[i]) {
  rows<-which(owner==i)
  for(end in c("start","end")) {
   j<-if(end=="start")head(rows,1) else tail(rows,1)
   point<-if(end=="start")z[j,1:2] else z[j,3:4]
   direction<-if(end=="start")z[j,1:2]-z[j,3:4] else z[j,3:4]-z[j,1:2]
   node<-if(end=="start")arr$first_node[j] else arr$last_node[j]
   if(arr$nodes$degree[node]==1L)endpoints[[length(endpoints)+1L]]<-list(contour=i,end=end,
     point=point,direction=direction/sqrt(sum(direction^2)),node=node)
  }
 }
 out<-list();n<-length(endpoints)
 if(n>1L)for(i in seq_len(n-1L))for(j in seq.int(i+1L,n)) {
  a<-endpoints[[i]];b<-endpoints[[j]];v<-b$point-a$point;gap<-sqrt(sum(v*v))
  if(a$node==b$node||gap<=tolerance||gap>gap_max)next
  direction<-v/gap
  angle<-max(acos(pmax(-1,pmin(1,c(sum(a$direction*direction),-sum(b$direction*direction))))))*180/pi
  if(angle<=angle_tolerance+1e-10)out[[length(out)+1L]]<-data.frame(
   contour_a=a$contour,end_a=a$end,contour_b=b$contour,end_b=b$end,gap=gap,angle=angle)
 }
 if(length(out))do.call(rbind,out) else data.frame(contour_a=integer(),end_a=character(),contour_b=integer(),end_b=character(),gap=double(),angle=double())
}
.graph_parallel <- function(z,owner,hits,angle_tolerance,max_distance,min_extent) {
 out<-list();n<-nrow(z)
 if(n>1L)for(i in seq_len(n-1L))for(j in seq.int(i+1L,n)) {
  if(owner[i]==owner[j]||hits[i,j])next
  a<-z[i,1:2];b<-z[j,1:2];u<-z[i,3:4]-a;v<-z[j,3:4]-b
  nu<-u/sqrt(sum(u*u));nv<-v/sqrt(sum(v*v));cosine<-sum(nu*nv)
  angle<-acos(pmin(1,abs(cosine)))*180/pi
  if(angle>angle_tolerance+1e-10)next
  if(cosine<0)nv<--nv
  axis<-nu+nv;axis<-axis/sqrt(sum(axis*axis))
  ia<-sort(c(0,sum(u*axis)));ib<-sort(c(sum((b-a)*axis),sum((b+v-a)*axis)))
  low<-max(ia[1],ib[1]);high<-min(ia[2],ib[2]);extent<-high-low
  if(extent<=0||extent<min_extent-1e-10)next
  separation<-vapply(c(low,high),function(t) {
   p<-a+t/sum(u*axis)*u;q<-b+(t-sum((b-a)*axis))/sum(v*axis)*v
   abs(.cross2(axis,q-p))
  },double(1))
  if(mean(separation)>max_distance+1e-10)next
  out[[length(out)+1L]]<-data.frame(segment_a=i,segment_b=j,contour_a=owner[i],contour_b=owner[j],
    shared_extent=extent,separation=mean(separation),angle=angle)
 }
 if(length(out))do.call(rbind,out) else data.frame(segment_a=integer(),segment_b=integer(),contour_a=integer(),contour_b=integer(),shared_extent=double(),separation=double(),angle=double())
}
.graph_inside <- function(point,segments) {
 a<-segments[,1:2,drop=FALSE];b<-segments[,3:4,drop=FALSE]
 crossing<-which((a[,2]>point[2])!=(b[,2]>point[2]))
 if(!length(crossing))return(FALSE)
 x<-a[crossing,1]+(point[2]-a[crossing,2])*(b[crossing,1]-a[crossing,1])/(b[crossing,2]-a[crossing,2])
 sum(point[1]<x)%%2==1L
}
.graph_clearance <- function(a,b) {
 distance<-function(points,segments) {
  v<-segments[,3:4,drop=FALSE]-segments[,1:2,drop=FALSE];denominator<-rowSums(v*v)
  min(vapply(seq_len(nrow(points)),function(i) {
   w<-matrix(points[i,],nrow(v),2,byrow=TRUE)-segments[,1:2,drop=FALSE]
   t<-pmax(0,pmin(1,rowSums(w*v)/denominator))
   min(sqrt(rowSums((w-v*t)^2)))
  },double(1)))
 }
 min(distance(rbind(a[,1:2,drop=FALSE],a[,3:4,drop=FALSE]),b),
     distance(rbind(b[,1:2,drop=FALSE],b[,3:4,drop=FALSE]),a))
}

#' Construct a contour arrangement and relationship graph
#'
#' Exact crossings, including self-crossings, split source spans. Nearby
#' endpoints/intersections merge transitively; collinear overlaps error. Bends
#' remain degree-two nodes. Enclosure applies only to disjoint simple closed
#' polygons. Continuation links free endpoints; parallel relations use shared
#' projected extent and mean normal separation. These geometric relations do not
#' establish border ownership, global symmetry or perceptual validity.
#' @param x A drawing.
#' @param spec Optional feature specification; uses the drawing's saved recipe
#'   or a default pixel recipe. Distances use its target spatial units.
#' @param tolerance Positive node-clustering tolerance and enclosure clearance.
#' @param gap_max Maximum continuation gap.
#' @param angle_tolerance Tangent/axial angular tolerance in geometric degrees,
#'   nonnegative and strictly less than 90.
#' @param parallel_distance Maximum mean normal separation of parallel spans.
#' @param min_shared_extent Minimum projected overlap for parallel spans.
#' @param pixels_per_degree Optional input calibration override.
#' @return An `mlv_contour_graph` with nodes, source-span edges, continuation,
#'   parallel and enclosure tables, method specification and source geometry.
#'   At most 1,000 source segments and 20,000 candidate split points are accepted.
#'   Node coordinates are in target spatial units; edge lengths retain source
#'   geometry even when clustering moves the representative node centroid.
#' @export
contour_graph <- function(x,spec=NULL,tolerance=1e-6,gap_max=3,angle_tolerance=5,
 parallel_distance=6,min_shared_extent=1,pixels_per_degree=NULL) {
 .check_drawing(x)
 if(is.null(spec))spec<-if(is.null(x$feature_spec))feature_spec() else x$feature_spec
 validate_feature_spec(spec)
 if(!is.null(pixels_per_degree))spec$settings$pixels_per_degree<-.spec_number(pixels_per_degree,"pixels_per_degree",strict=TRUE)
 tolerance<-.spec_number(tolerance,"tolerance",strict=TRUE)
 gap_max<-.spec_number(gap_max,"gap_max");angle_tolerance<-.spec_number(angle_tolerance,"angle_tolerance",upper=90)
 if(angle_tolerance==90)stop("angle_tolerance must be less than 90",call.=FALSE)
 parallel_distance<-.spec_number(parallel_distance,"parallel_distance")
 min_shared_extent<-.spec_number(min_shared_extent,"min_shared_extent")
 st<-.segment_table(x)
 if(nrow(st)>1000)stop("graph supports at most 1000 source segments",call.=FALSE)
 factor<-.spec_scale(spec);z<-as.matrix(st[,c("x1","y1","x2","y2")])*factor
 owner<-st$contour_index
 closed<-vapply(x$contours,function(z)all(z[1,1:2]==z[nrow(z),3:4]),logical(1))
 arr<-.graph_arrangement(z,owner,st$segment_index,closed,tolerance)
 edges<-arr$edges;edges$edge_id<-seq_len(nrow(edges))
 edges$contour_index<-owner[edges$source_segment];edges$contour_id<-x$ids[edges$contour_index]
 edges$segment_index<-st$segment_index[edges$source_segment]
 continuation<-.graph_continuation(z,owner,closed,arr,gap_max,angle_tolerance,tolerance)
 parallel<-.graph_parallel(z,owner,arr$hits,angle_tolerance,parallel_distance,min_shared_extent)
 add_ids<-function(df) {
  df$contour_id_a<-x$ids[df$contour_a];df$contour_id_b<-x$ids[df$contour_b];df
 }
 continuation<-add_ids(continuation);parallel<-add_ids(parallel)
 simple<-closed & !arr$self_crossing;enclosures<-list()
 for(i in which(simple))for(j in setdiff(which(simple),i)) {
  ai<-which(owner==i);bi<-which(owner==j)
  if(any(arr$hits[ai,bi,drop=FALSE])||!.graph_inside(z[bi[1],1:2],z[ai,,drop=FALSE]))next
  clearance<-.graph_clearance(z[ai,,drop=FALSE],z[bi,,drop=FALSE])
  if(clearance<=tolerance)next
  area<-function(k)abs(sum(z[k,1]*z[k,4]-z[k,3]*z[k,2]))/2
  enclosures[[length(enclosures)+1L]]<-data.frame(outer=i,inner=j,outer_id=x$ids[i],inner_id=x$ids[j],
      outer_area=area(ai),inner_area=area(bi),clearance=clearance)
 }
 enclosure<-if(length(enclosures))do.call(rbind,enclosures) else
  data.frame(outer=integer(),inner=integer(),outer_id=character(),inner_id=character(),outer_area=double(),inner_area=double(),clearance=double())
 parameters<-list(tolerance=tolerance,gap_max=gap_max,angle_tolerance=angle_tolerance,
   parallel_distance=parallel_distance,min_shared_extent=min_shared_extent)
 structure(list(nodes=arr$nodes,edges=edges,continuation=continuation,parallel=parallel,enclosure=enclosure,
  contours=data.frame(contour_id=x$ids,closed=closed,simple_closed=simple),
  collapsed_spans=arr$collapsed,spec=spec,parameters=parameters,source_geometry=.geometry_snapshot(x),
  source_provenance=x$source,method="arrangement-relations-1",
  contract=list(parameters=parameters,spec=.spec_semantics(spec),pixel_size=factor)),class="mlv_contour_graph")
}
.check_contour_graph <- function(graph,x=NULL,spec=NULL) {
 if(!inherits(graph,"mlv_contour_graph")||!identical(graph$method,"arrangement-relations-1"))
  stop("unsupported contour graph",call.=FALSE)
 validate_feature_spec(graph$spec)
 if(!identical(graph$contract,list(parameters=graph$parameters,spec=.spec_semantics(graph$spec),pixel_size=.spec_scale(graph$spec))))
  stop("graph specification changed; recompute contour_graph()",call.=FALSE)
 if(!is.null(x)&&!identical(graph$source_geometry,.geometry_snapshot(x)))
  stop("stale graph: source geometry differs",call.=FALSE)
 if(!is.null(spec)&&(!identical(.spec_semantics(spec),.spec_semantics(graph$spec))||
                    !identical(.spec_scale(spec),.spec_scale(graph$spec))))
  stop("graph feature specification mismatch",call.=FALSE)
 invisible(graph)
}
#' Summarize topology and geometric contour relations
#' @param graph A contour graph.
#' @param x Optional source drawing to check for stale geometry.
#' @param spec Optional requested specification for compatibility validation.
#' @return A one-row data frame with graph sizes, components, cycle rank,
#'   degree counts, closure, continuation, parallel extent/separation and enclosure.
#'   Mean parallel separation is NA when no parallel relations exist. Raw vertex
#'   and edge counts depend on retained degree-two vertices; cycle rank does not.
#' @export
graph_features <- function(graph,x=NULL,spec=NULL) {
 .check_contour_graph(graph,x,spec)
 nodes<-graph$nodes;edges<-graph$edges;parallel<-graph$parallel
 components<-length(unique(nodes$component));extent<-sum(parallel$shared_extent)
 result<-data.frame(vertices=nrow(nodes),edges=nrow(edges),components=components,
  cycles=nrow(edges)-nrow(nodes)+components,endpoints=sum(nodes$degree==1),branch_nodes=sum(nodes$degree>=3),
  closed_contours=sum(graph$contours$closed),simple_closed_contours=sum(graph$contours$simple_closed),
  continuation_pairs=nrow(graph$continuation),parallel_contour_pairs=nrow(unique(parallel[,c("contour_a","contour_b")])),
  parallel_shared_extent=extent,mean_parallel_separation=if(extent>0)sum(parallel$shared_extent*parallel$separation)/extent else NA_real_,
  enclosures=nrow(graph$enclosure),collapsed_spans=graph$collapsed_spans)
 attr(result,"graph_contract")<-c(list(method=graph$method),graph$contract[c("parameters","spec")])
 attr(result,"input_pixel_size")<-graph$contract$pixel_size;result
}
#' Export graph spans as a data frame
#' @param x A contour graph.
#' @param ... Unused.
#' @return The edge table, retaining source contour and segment IDs.
#' @export
as.data.frame.mlv_contour_graph <- function(x,...) {
 .check_contour_graph(x);out<-x$edges
 attr(out,"graph_contract")<-c(list(method=x$method),x$contract);out
}
