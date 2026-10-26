# Isolated candidates: checked archive namespace is never modified.
# Usage: Rscript tools/performance/histogram_spike.R LIB OUTPUT WORK
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
dest<-args[2];work<-args[3]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
dir.create(work,recursive=TRUE,showWarnings=FALSE)
native<-new.env(parent=asNamespace("mlvtoolbox"))
Rcpp::sourceCpp("tools/performance/histogram_candidates.cpp",env=native,
 cacheDir=file.path(work,"compile"),showOutput=TRUE)
baseline<-asNamespace("mlvtoolbox")
candidate<-function(lengths=FALSE) {
 env<-new.env(parent=baseline)
 env$grouped_sums<-native$grouped_sums
 env$contour_lengths<-native$contour_lengths
 for(name in c(".statistics_blocks","feature_table","contour_statistics")) {
  f<-get(name,baseline);environment(f)<-env;assign(name,f,env)
 }
 f<-function(values,weights,groups,edges,n_groups,index=NULL) {
  k<-length(edges)-1L
  if(as.double(n_groups)*k>.Machine$integer.max)stop("histogram exceeds index capacity",call.=FALSE)
  if(!length(values)||!n_groups)return(matrix(0,n_groups,k))
  if(is.null(index))index<-findInterval(values,edges,all.inside=TRUE)
  grouped_sums(weights,groups,index,n_groups,k)
 }
 environment(f)<-env;env$.grouped_bins<-f
 if(lengths) {
  b<-body(env$.statistics_blocks)
  at<-which(vapply(as.list(b),function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&
      identical(x[[2]],as.name("lens")),logical(1)))
  stopifnot(length(at)==1L);b[[at]]<-quote(lens<-contour_lengths(x$contours))
  body(env$.statistics_blocks)<-b
 }
 env
}
grouped<-candidate();both<-candidate(TRUE)
assertions<-0L
same<-function(a,b){stopifnot(identical(a,b));assertions<<-assertions+1L}
# Existing independent checked-output witnesses, including MAT/junction properties.
w<-readRDS(system.file("fixtures","statistics-before-optimization.rds",package="mlvtoolbox"))
for(e in list(grouped,both)) {
 for(case in w$cases)same(do.call(e$contour_statistics,c(list(w$drawing),case$args)),case$expected)
 same(e$feature_table(w$drawing,c("length","orientation","curvature","parallelism"),13),w$features)
}
# Adversarial accumulation: empty/unoccupied bins, missing weights, cancellations,
# repeated groups, input-order preservation and wide magnitude ranges.
set.seed(2751)
for(i in 1:150) {
 n<-sample(0:1500,1);ng<-sample(1:50,1);nb<-sample(1:20,1)
 g<-sample(seq_len(ng),n,TRUE);b<-sample(seq_len(nb),n,TRUE)
 weights<-sample(c(0,NA_real_,NaN,1,-1,1e16,-1e16,1e-12),n,TRUE)*runif(n)
 same(native$grouped_sums(weights,g,b,ng,nb),
      baseline$.grouped_bins(rep(0,n),weights,g,seq_len(nb+1),ng,b))
}
# All five executed contour fixtures, including axial bin boundaries.
for(name in c("integer0","integer1","integer2","subpixel","fractional")) {
 z<-as.matrix(read.csv(system.file("fixtures","parity",paste0(name,".input.csv"),package="mlvtoolbox"),header=FALSE))
 d<-mlv_drawing(lapply(unique(z[,1]),function(i)z[z[,1]==i,2:5,drop=FALSE]),256,192)
 for(bins in c(1,5,8,13,31))for(norm in c(FALSE,TRUE))for(by in c("drawing","contour")) {
  props<-c("orientation","length","curvature","horizontal_vertical")
  same(both$contour_statistics(d,props,bins,by,norm),contour_statistics(d,props,bins,by,norm))
 }
}
# Optional retained forest file is SHA-verified by the invoking runner.
forest<-read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat",zero_length="drop")
stopifnot(length(forest$contours)==562)
stopifnot(sum(vapply(forest$contours,nrow,integer(1)))==2072)
write.table(do.call(rbind,lapply(seq_along(forest$contours),function(i)cbind(i,forest$contours[[i]]))),
 file.path(work,"forest.csv"),sep=",",row.names=FALSE,col.names=FALSE)
write.table(c(forest$width,forest$height),file.path(work,"dimensions.csv"),sep=",",row.names=FALSE,col.names=FALSE)
make_drawing<-function(n,segments)mlv_drawing(lapply(seq_len(n),function(i) {
 x<-seq(2,20,length.out=segments+1)+(i-1)%%32*30
 y<-10+3*sin(seq(0,2*pi,length.out=segments+1))+(i-1)%/%32*30
 cbind(head(x,-1),head(y,-1),tail(x,-1),tail(y,-1))
}),1000,ceiling(n/32)*30+20)
drawings<-list(forest=forest,many_contours=make_drawing(1000,20),one_long_contour=make_drawing(1,20000))
props<-c("orientation","length","curvature")
jobs<-list(R_0.1.1=feature_table,cpp_grouping=grouped$feature_table,cpp_grouping_and_lengths=both$feature_table)
results<-list()
for(name in names(drawings)) {
 d<-drawings[[name]];expected<-feature_table(d,props,bins=8)
 for(f in jobs)same(f(d,props,bins=8),expected)
 same(native$contour_lengths(d$contours),
      vapply(d$contours,function(z)sum(sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))),double(1)))
 results[[name]]<-expected
}
# Alternate order over three rounds; five warm batches per entry.
timer<-function(f) {
 f();number<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(number))f())[3])
  if(elapsed>=.12||number>=512)break
  number<-number*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(number))f())[3])/number)
 c(median_seconds=median(times),min_seconds=min(times),max_seconds=max(times),iterations=number)
}
rows<-list()
for(round in 1:3)for(name in names(drawings)) {
 order<-if(round%%2)names(jobs) else rev(names(jobs))
 for(method in order) {
  f<-jobs[[method]];d<-drawings[[name]]
  rows[[length(rows)+1L]]<-data.frame(round,workload=name,method,
    as.list(timer(function()f(d,props,bins=8))))
  cat(round,name,method,tail(rows,1)[[1]]$median_seconds,"seconds\n")
 }
}
out<-do.call(rbind,rows);write.csv(out,file.path(dest,"timings.csv"),row.names=FALSE)
# Long enough profiles to distinguish the remaining overhead.
for(method in c("R_0.1.1","cpp_grouping_and_lengths")) {
 f<-jobs[[method]]
 Rprof(file.path(work,paste0(method,".Rprof")),interval=.001)
 for(i in 1:150)f(forest,props,bins=8)
 Rprof(NULL)
 capture.output(summaryRprof(file.path(work,paste0(method,".Rprof")))$by.total,
                file=file.path(dest,paste0(method,"-profile.txt")))
}
actual<-as.matrix(results$forest[,-1])
write.table(actual,file.path(work,"R-histogram.csv"),sep=",",row.names=FALSE,col.names=FALSE)
jsonlite::write_json(list(status="passed",exact_comparisons=assertions,
 accumulation_bytes=native$accumulator_bytes(),R=R.version.string,
 package_version=as.character(packageVersion("mlvtoolbox")),
 checked_source="d83ed952fdc9f639b98e849b0d8697259894de46",
 seed=2751,forest_contours=562,forest_segments=2072,
 note="Experimental closures and sourceCpp functions only; installed namespace unchanged."),
 file.path(dest,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
print(out)
