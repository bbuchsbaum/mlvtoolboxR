source(system.file("fixtures", "designer-v1.R", package="mlvtoolbox"), local=TRUE)
# Independent evaluation/comparator for the frozen v1 benchmark.
# Does not call the optimizer, its constraint code, or its feature-matrix builder.
reference_measure <- function(drawings,spec,target) {
 result<-lapply(seq_along(drawings),function(i) {
  x<-drawings[[i]];st<-as.data.frame(x,spec=spec)
  score<-vapply(seq_along(x$ids),function(j) {
   z<-st[st$contour_index==j,]
   switch(target,integrated_curvature=sum(z$length*z$curvature),
    length=sum(z$length),orientation=sum(z$length*(abs(cos(z$orientation*pi/180))-abs(sin(z$orientation*pi/180)))))
  },double(1))
  f<-feature_table(x,spec=spec)
  length<-vapply(seq_along(x$ids),function(j)sum(st$length[st$contour_index==j]),double(1))
  data.frame(image=i,id=x$ids,length=length,score=score,f[,-1,drop=FALSE])
 })
 do.call(rbind,result)
}
reference_evaluate <- function(drawings,args,assignment) {
 f<-reference_measure(drawings,args$spec,args$target)
 stopifnot(length(assignment)==nrow(f),all(assignment%in%0:2))
 groups<-if(args$mode=="set")list(seq_len(nrow(f))) else split(seq_len(nrow(f)),f$image)
 checks<-lapply(groups,function(rows) {
  a<-intersect(rows,which(assignment==1));b<-intersect(rows,which(assignment==2))
  la<-sum(f$length[a]);lb<-sum(f$length[b]);ca<-length(a);cb<-length(b)
  residual<-vapply(names(args$nuisance),function(p) {
   if(p=="length")return(abs(la-lb))
   cols<-startsWith(names(f),paste0(p,"_"))
   sum(abs(colSums(f[a,cols,drop=FALSE])-colSums(f[b,cols,drop=FALSE])))
  },double(1))
  good<-all(c(la,lb)>=args$length_bounds[1]-1e-8)&all(c(la,lb)<=args$length_bounds[2]+1e-8)&
   all(c(ca,cb)>=args$count_bounds[1])&all(c(ca,cb)<=args$count_bounds[2])&
   abs(ca-cb)<=args$count_imbalance&all(residual<=args$nuisance+1e-8)
  list(feasible=good,length_A=la,length_B=lb,count_A=ca,count_B=cb,residual=residual)
 })
 list(feasible=all(vapply(checks,`[[`,logical(1),"feasible")),
  objective=sum(f$score[assignment==1])-sum(f$score[assignment==2]),groups=checks)
}
reference_ranked <- function(drawings,args) {
 f<-reference_measure(drawings,args$spec,args$target);assignment<-integer(nrow(f))
 groups<-if(args$mode=="set")list(seq_len(nrow(f))) else split(seq_len(nrow(f)),f$image)
 for(rows in groups) {
  # Whole-contour prefixes stop at the first over-budget contour. A is chosen
  # first, then B from remaining IDs; this explicitly prevents tied overlap.
  for(condition in 1:2) {
   eligible<-rows[assignment[rows]==0]
   order<-eligible[order(if(condition==1)-f$score[eligible] else f$score[eligible],
                         f$image[eligible],f$id[eligible])]
   order<-head(order,args$count_bounds[2])
   selected<-order[cumsum(f$length[order])<=args$length_bounds[2]+1e-8]
   assignment[selected]<-condition
  }
 }
 list(assignment=assignment,evaluation=reference_evaluate(drawings,args,assignment))
}
reference_exhaustive <- function(drawings,args) {
 n<-sum(vapply(drawings,function(x)length(x$ids),integer(1)))
 stopifnot(n<=8)
 candidates<-as.matrix(expand.grid(rep(list(0:2),n)))
 # This simple oracle uses its own vectorized constraint computation, independently
 # of the eventual depth-first optimizer. Rows enumerate all ternary assignments.
 f<-reference_measure(drawings,args$spec,args$target);a<-candidates==1;b<-candidates==2
 feasible<-rep(TRUE,nrow(candidates))
 groups<-if(args$mode=="set")list(seq_len(n)) else split(seq_len(n),f$image)
 for(rows in groups) {
  la<-as.numeric(a[,rows,drop=FALSE]%*%f$length[rows]);lb<-as.numeric(b[,rows,drop=FALSE]%*%f$length[rows])
  ca<-rowSums(a[,rows,drop=FALSE]);cb<-rowSums(b[,rows,drop=FALSE])
  feasible<-feasible & la>=args$length_bounds[1]-1e-8 & lb>=args$length_bounds[1]-1e-8 &
   la<=args$length_bounds[2]+1e-8 & lb<=args$length_bounds[2]+1e-8 &
   ca>=args$count_bounds[1]&cb>=args$count_bounds[1]&ca<=args$count_bounds[2]&cb<=args$count_bounds[2]&abs(ca-cb)<=args$count_imbalance
  for(p in names(args$nuisance)) {
   value<-if(p=="length")abs(la-lb) else {
    cols<-startsWith(names(f),paste0(p,"_"));m<-as.matrix(f[rows,cols,drop=FALSE])
    rowSums(abs((a[,rows,drop=FALSE]-b[,rows,drop=FALSE])%*%m))
   }
   feasible<-feasible & value<=args$nuisance[[p]]+1e-8
  }
 }
 scores<-as.numeric((a-b)%*%f$score)
 list(feasible=sum(feasible),optimum=if(any(feasible))max(scores[feasible]) else NA_real_)
}
