.subset_prepare <- function(drawings,spec,target,length_bounds,count_bounds,nuisance,
                            count_imbalance,mode,pixels_per_degree=NULL) {
 if(inherits(drawings,"mlv_drawing"))drawings<-list(image=drawings)
 if(!is.list(drawings)||!length(drawings))stop("supply a nonempty drawing list",call.=FALSE)
 invisible(lapply(drawings,.check_drawing));validate_feature_spec(spec)
 allowed<-c("length","orientation","curvature","horizontal_vertical")
 if(spec$settings$normalize||any(!spec$settings$properties%in%allowed))
  stop("designer requires normalize=FALSE and additive contour features only",call.=FALSE)
 target<-match.arg(target,c("integrated_curvature","length","orientation"))
 mode<-match.arg(mode,c("set","per_drawing"))
 ids<-names(drawings)
 if(is.null(ids))ids<-as.character(seq_along(drawings))
 if(anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids))stop("drawing IDs must be nonempty and unique",call.=FALSE)
 names(drawings)<-ids
 if(!is.numeric(length_bounds)||length(length_bounds)!=2||any(!is.finite(length_bounds))||
    length_bounds[1]<0||length_bounds[2]<length_bounds[1])stop("invalid length_bounds",call.=FALSE)
 if(!is.numeric(count_bounds)||length(count_bounds)!=2||anyNA(count_bounds)||
    !is.finite(count_bounds[1])||count_bounds[1]<1||count_bounds[2]<count_bounds[1]||
    any(count_bounds!=floor(count_bounds)))stop("invalid count_bounds",call.=FALSE)
 if(!is.numeric(count_imbalance)||length(count_imbalance)!=1||is.na(count_imbalance)||
    count_imbalance<0||count_imbalance!=floor(count_imbalance))stop("invalid count_imbalance",call.=FALSE)
 if(!is.numeric(nuisance)||any(!is.finite(nuisance))||any(nuisance<0)||
    (length(nuisance)&&(is.null(names(nuisance))||anyDuplicated(names(nuisance))||
     any(!names(nuisance)%in%allowed))))stop("invalid nuisance tolerances",call.=FALSE)
 if(any(!setdiff(names(nuisance),"length")%in%spec$settings$properties))
  stop("nuisance histogram properties must be included in spec",call.=FALSE)
 if(!is.null(pixels_per_degree)&&(!is.numeric(pixels_per_degree)||
    length(pixels_per_degree)!=length(drawings)||any(!is.finite(pixels_per_degree))||any(pixels_per_degree<=0)))
  stop("pixels_per_degree must supply one positive calibration per drawing",call.=FALSE)
 n<-sum(vapply(drawings,function(x)length(x$ids),integer(1)))
 if(n>64)stop("bounded exact designer supports at most 64 eligible contours",call.=FALSE)
 prepared<-lapply(seq_along(drawings),function(i)compute_contour_properties(drawings[[i]],spec=spec,
   pixels_per_degree=if(is.null(pixels_per_degree))NULL else pixels_per_degree[i]))
 names(prepared)<-ids
 tables<-lapply(prepared,feature_table)
 for(f in tables).compare_feature_schemas(feature_schema(tables[[1]]),feature_schema(f))
 rows<-do.call(rbind,lapply(seq_along(prepared),function(i) {
  x<-prepared[[i]]
  data.frame(image_index=rep(i,length(x$ids)),image_id=rep(ids[i],length(x$ids)),contour_id=x$ids,
   length=.contour_criterion(x,"length"),score=.contour_criterion(x,target)*
     if(target=="orientation").contour_criterion(x,"length") else 1)
 }))
 features<-do.call(rbind,lapply(tables,function(f)f[,-1,drop=FALSE]))
 columns<-lapply(setdiff(names(nuisance),"length"),function(p)which(startsWith(names(features),paste0(p,"_"))))
 names(columns)<-setdiff(names(nuisance),"length")
 hist<-lapply(columns,function(j)as.matrix(features[,j,drop=FALSE]))
 list(drawings=prepared,rows=rows,hist=hist,spec=spec,target=target,
  length_bounds=length_bounds,count_bounds=count_bounds,nuisance=nuisance,count_imbalance=count_imbalance,
  mode=mode,groups=if(mode=="set")rep.int(1L,n) else rows$image_index,
  n_groups=if(mode=="set")1L else length(drawings),tolerance=1e-8)
}

.subset_residual <- function(problem,assignment) {
 p<-problem;r<-p$rows
 evaluate<-function(rows,label) {
  a<-rows[assignment[rows]==1L];b<-rows[assignment[rows]==2L]
  la<-sum(r$length[a]);lb<-sum(r$length[b]);na<-length(a);nb<-length(b)
  residual<-vapply(names(p$nuisance),function(prop) {
   if(prop=="length")return(abs(la-lb))
   h<-p$hist[[prop]];sum(abs(colSums(h[a,,drop=FALSE])-colSums(h[b,,drop=FALSE])))
  },double(1))
  valid<-all(c(la,lb)>=p$length_bounds[1]-p$tolerance)&all(c(la,lb)<=p$length_bounds[2]+p$tolerance)&
    all(c(na,nb)>=p$count_bounds[1])&all(c(na,nb)<=p$count_bounds[2])&
    abs(na-nb)<=p$count_imbalance&all(residual<=p$nuisance+p$tolerance)
  list(group=label,feasible=valid,length_A=la,length_B=lb,count_A=na,count_B=nb,
       objective=sum(r$score[a])-sum(r$score[b]),residual=residual)
 }
 groups<-lapply(seq_len(p$n_groups),function(g)evaluate(which(p$groups==g),
  if(p$mode=="set")"set" else names(p$drawings)[g]))
 per_drawing<-lapply(seq_along(p$drawings),function(i)evaluate(which(r$image_index==i),names(p$drawings)[i]))
 list(feasible=all(vapply(groups,`[[`,logical(1),"feasible")),
      objective=sum(r$score[assignment==1L])-sum(r$score[assignment==2L]),
      groups=groups,per_drawing=per_drawing)
}

.subset_solve <- function(problem,time_limit=5,max_nodes=100000) {
 time_limit<-.spec_number(time_limit,"time_limit")
 max_nodes<-.spec_number(max_nodes,"max_nodes",strict=TRUE,integer=TRUE)
 p<-problem;r<-p$rows;n<-nrow(r);ng<-p$n_groups;tol<-p$tolerance
 # Stable traversal; first optimum found wins numerical ties.
 order<-order(-abs(r$score),r$image_id,r$contour_id)
 score<-r$score[order];length<-r$length[order];group<-p$groups[order]
 remaining<-c(rev(cumsum(rev(abs(score)))),0)
 count_left<-len_left<-matrix(0,n+1L,ng)
 if(n)for(i in n:1L) {
  count_left[i,]<-count_left[i+1L,];len_left[i,]<-len_left[i+1L,]
  count_left[i,group[i]]<-count_left[i,group[i]]+1
  len_left[i,group[i]]<-len_left[i,group[i]]+length[i]
 }
 selected<-integer(n);best<-NULL;best_score<--Inf;nodes<-0L;stopped<-FALSE
 start<-proc.time()[['elapsed']]
 walk<-function(i,la,lb,ca,cb,value) {
  if(stopped)return(invisible(NULL))
  # The clock is read every 256 nodes; proc.time() dominates small nodes.
  if(nodes>=max_nodes||(nodes%%256L==0L&&proc.time()[['elapsed']]-start>=time_limit)) {stopped<<-TRUE;return(invisible(NULL))}
  nodes<<-nodes+1L
  if(any(la>p$length_bounds[2]+tol)||any(lb>p$length_bounds[2]+tol)||
     any(ca>p$count_bounds[2])||any(cb>p$count_bounds[2]))return(invisible(NULL))
  if(any(la+len_left[i,]<p$length_bounds[1]-tol)||any(lb+len_left[i,]<p$length_bounds[1]-tol)||
     any(ca+count_left[i,]<p$count_bounds[1])||any(cb+count_left[i,]<p$count_bounds[1]))return(invisible(NULL))
  if(value+remaining[i]<=best_score+tol)return(invisible(NULL))
  if(i>n) {
   original<-integer(n);original[order]<-selected
   result<-.subset_residual(p,original)
   if(result$feasible) {best<<-original;best_score<<-result$objective}
   return(invisible(NULL))
  }
  g<-group[i]
  choices<-if(score[i]>=0)c(1L,2L,0L) else c(2L,1L,0L)
  for(choice in choices) {
   selected[i]<<-choice
   aa<-la;bb<-lb;ac<-ca;bc<-cb
   if(choice==1L){aa[g]<-aa[g]+length[i];ac[g]<-ac[g]+1L}
   if(choice==2L){bb[g]<-bb[g]+length[i];bc[g]<-bc[g]+1L}
   walk(i+1L,aa,bb,ac,bc,value+if(choice==1L)score[i] else if(choice==2L)-score[i] else 0)
   if(stopped)break
  }
  invisible(NULL)
 }
 walk(1L,numeric(ng),numeric(ng),integer(ng),integer(ng),0)
 elapsed<-proc.time()[['elapsed']]-start
 status<-if(stopped)if(is.null(best))"budget_exhausted" else "feasible" else
   if(is.null(best))"infeasible" else "optimal"
 upper<-if(!stopped)if(is.null(best))NA_real_ else best_score else remaining[1]
 list(assignment=best,status=status,objective=if(is.null(best))NA_real_ else best_score,
  upper_bound=upper,absolute_gap=if(is.null(best))NA_real_ else max(0,upper-best_score),
  nodes=nodes,search_seconds=elapsed,time_limit=time_limit,max_nodes=max_nodes,
  working_data_bytes=as.numeric(utils::object.size(p)),solver="mlv-bounded-exact-1",
  numerical_tolerance=tol,tie_policy="absolute target descending, drawing/contour ID, preferred condition then unused")
}

.subset_verify_outputs <- function(p,conditions) {
 calibration<-if(p$spec$settings$units=="visual_degrees")
  vapply(p$drawings,function(x)x$feature_spec$settings$pixels_per_degree,double(1)) else NULL
 prepare<-function(drawings).subset_prepare(drawings,p$spec,p$target,p$length_bounds,
  p$count_bounds,p$nuisance,p$count_imbalance,p$mode,calibration)
 a<-prepare(conditions$A);b<-prepare(conditions$B)
 rebuilt<-p;rebuilt$rows<-rbind(a$rows,b$rows)
 rebuilt$hist<-Map(rbind,a$hist,b$hist)
 rebuilt$groups<-if(p$mode=="set")rep.int(1L,nrow(rebuilt$rows)) else rebuilt$rows$image_index
 assignment<-c(rep.int(1L,nrow(a$rows)),rep.int(2L,nrow(b$rows)))
 .subset_residual(rebuilt,assignment)
}

#' Design two matched conditions from whole contours
#'
#' Assigns contours to condition A, B, or unassigned, maximizing an additive
#' target difference. A deterministic bounded exact search is intended for small
#' problems (at most 64 eligible contours). It adds no solver dependency. Set
#' matching is primary; per-drawing mode applies all bounds separately to every
#' drawing and can be infeasible even when a pooled experiment is feasible.
#'
#' Nuisance constraints compare raw masses: length difference or L1 distance
#' between length-weighted histogram bins. Independently normalized histograms,
#' transformations, junction/MAT constraints and spatial union coverage are not
#' supported. Reconstruction recomputes features and verifies all constraints.
#' @param drawings Drawing or named list of drawings. All contours are eligible;
#'   subset inputs first to exclude any. Drawing IDs must be unique.
#' @param spec A shared feature specification with normalize=FALSE and additive
#'   contour properties only. It must include requested nuisance histograms.
#' @param target Additive contour target: integrated absolute curvature at the
#'   selected scale, length, or length-integrated horizontal-minus-vertical value.
#' @param length_bounds Common lower/upper retained length for each condition,
#'   in specification spatial units. Finite, increasing, nonnegative values.
#' @param count_bounds Common positive lower/upper whole-contour counts.
#' @param nuisance Named nonnegative tolerances. `length` bounds absolute length
#'   difference. `orientation`, `curvature`, and `horizontal_vertical` bound raw
#'   histogram L1 distance. An empty numeric vector omits nuisance constraints.
#' @param count_imbalance Maximum absolute count difference, or Inf.
#' @param mode Match over the complete set or separately within every drawing.
#' @param time_limit Search wall-time budget in seconds; zero performs no search.
#' @param max_nodes Search node budget. Preprocessing/verification are separate.
#' @param pixels_per_degree Optional positive calibration for each input drawing.
#' @param diagnostics Optional names: junctions, parallelism, separation, mirror,
#'   taper. These are recomputed on selected drawings and are never optimized.
#' @return An `mlv_subset_design` containing condition drawing sets, unassigned
#'   drawings, original contour/source IDs, specification and input provenance,
#'   recomputed residuals, and solver/resource metadata. Status is `optimal`,
#'   `infeasible`, `feasible` (budget stopped with a solution), or
#'   `budget_exhausted` (no solution found). A stopped run's upper bound may be
#'   loose; no solution is not an infeasibility proof. Unsupported inputs error.
#'   `working_data_bytes` is an object-size estimate, not peak process memory.
#'   Use saveRDS/readRDS to preserve the complete result.
#' @export
design_subsets <- function(drawings,spec,target="integrated_curvature",length_bounds,
 count_bounds=c(1,Inf),nuisance=c(length=0),count_imbalance=Inf,
 mode=c("set","per_drawing"),time_limit=5,max_nodes=100000,
 pixels_per_degree=NULL,diagnostics=character()) {
 started<-proc.time()[['elapsed']];mode<-match.arg(mode)
 if(!is.character(diagnostics)||anyNA(diagnostics)||anyDuplicated(diagnostics)||
    any(!diagnostics%in%c("junctions","parallelism","separation","mirror","taper")))
  stop("invalid diagnostics",call.=FALSE)
 p<-.subset_prepare(drawings,spec,target,length_bounds,count_bounds,nuisance,count_imbalance,mode,pixels_per_degree)
 # MAT diagnostics run on each drawing's own specification with the requested
 # properties; build those specifications now so an invalid combination fails
 # before the search rather than after it.
 mat_diagnostics<-setdiff(diagnostics,"junctions")
 diagnostic_spec<-function(x) {
  args<-x$feature_spec$settings;args$properties<-mat_diagnostics
  do.call(feature_spec,args)
 }
 if(length(mat_diagnostics))for(x in p$drawings)diagnostic_spec(x)
 prepared_at<-proc.time()[['elapsed']]
 result<-.subset_solve(p,time_limit,max_nodes)
 assignment<-result$assignment;result$assignment<-NULL
 selection<-p$rows[,c("image_id","contour_id"),drop=FALSE]
 selection$condition<-if(is.null(assignment))rep(NA_character_,nrow(selection)) else c("unassigned","A","B")[assignment+1L]
 conditions<-unassigned<-verified<-relational<-NULL
 if(!is.null(assignment)) {
  extract<-function(value) {
   out<-lapply(seq_along(p$drawings),function(i) {
    rows<-which(p$rows$image_index==i)
    subset_contours(p$drawings[[i]],which(assignment[rows]==value))
   });names(out)<-names(p$drawings);out
  }
  conditions<-list(A=extract(1L),B=extract(2L));unassigned<-extract(0L)
  verified<-.subset_verify_outputs(p,conditions)
  if(!verified$feasible||abs(verified$objective-result$objective)>p$tolerance)
   stop("internal error: reconstructed subsets failed independent verification",call.=FALSE)
  if(length(diagnostics))relational<-lapply(conditions,function(set)lapply(set,function(x) {
   out<-list()
   if("junctions"%in%diagnostics)out$junctions<-.measurement_junctions(x,x$feature_spec)
   if(length(mat_diagnostics))out$mat<-analyse_mat(x,spec=diagnostic_spec(x))
   out
  }))
 }
 result$preprocessing_seconds<-prepared_at-started
 result$total_seconds<-proc.time()[['elapsed']]-started
 structure(list(status=result$status,conditions=conditions,unassigned=unassigned,
  selection=selection,spec=spec,input_specifications=lapply(p$drawings,`[[`,"feature_spec"),
  input_provenance=lapply(p$drawings,`[[`,"source"),
  input_geometry=lapply(p$drawings,.geometry_snapshot),
  request=list(target=target,length_bounds=length_bounds,count_bounds=count_bounds,
   nuisance=nuisance,count_imbalance=count_imbalance,mode=mode,diagnostics=diagnostics),
  diagnostics=verified,relational_diagnostics=relational,solver=result),class="mlv_subset_design")
}

#' @export
print.mlv_subset_design <- function(x,...) {
 cat("Whole-contour subset design:",x$status,"\n")
 if(!is.null(x$diagnostics))cat("Target difference:",format(x$diagnostics$objective),
   " | upper bound:",format(x$solver$upper_bound),"\n")
 cat("Search nodes:",x$solver$nodes," | elapsed:",format(x$solver$search_seconds),"seconds\n")
 invisible(x)
}

#' @export
as.data.frame.mlv_subset_design <- function(x,...)x$selection
