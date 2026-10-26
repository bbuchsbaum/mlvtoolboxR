# Compare independently installed packages and the earlier grouping prototype.
# Usage: Rscript tools/performance/fused_statistics.R LIB OUT METHOD [prototype]
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=3)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
dest<-args[2];method<-args[3];dir.create(dest,recursive=TRUE,showWarnings=FALSE)
ns<-asNamespace("mlvtoolbox");feature<-feature_table
if(length(args)>=4 && args[4]=="prototype") {
 env<-new.env(parent=ns);native<-new.env(parent=ns)
 Rcpp::sourceCpp("tools/performance/histogram_candidates.cpp",env=native,
                cacheDir=file.path(dest,"compile"),showOutput=TRUE)
 env$grouped_sums<-native$grouped_sums;env$contour_lengths<-native$contour_lengths
 for(name in c(".statistics_blocks","feature_table","contour_statistics")) {
  f<-get(name,ns);environment(f)<-env;assign(name,f,env)
 }
 f<-function(values,weights,groups,edges,n_groups,index=NULL) {
  k<-length(edges)-1L
  if(as.double(n_groups)*k>.Machine$integer.max)stop("histogram exceeds index capacity")
  if(!length(values)||!n_groups)return(matrix(0,n_groups,k))
  if(is.null(index))index<-findInterval(values,edges,all.inside=TRUE)
  grouped_sums(weights,groups,index,n_groups,k)
 }
 environment(f)<-env;env$.grouped_bins<-f
 b<-body(env$.statistics_blocks)
 at<-which(vapply(as.list(b),function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&
                       identical(x[[2]],as.name("lens")),logical(1)))
 stopifnot(length(at)==1L);b[[at]]<-quote(lens<-contour_lengths(x$contours))
 body(env$.statistics_blocks)<-b
 feature<-env$feature_table
}
forest<-read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat",zero_length="drop")
make<-function(n,segments)mlv_drawing(lapply(seq_len(n),function(i) {
 x<-seq(2,20,length.out=segments+1)+(i-1)%%32*30
 y<-10+3*sin(seq(0,2*pi,length.out=segments+1))+(i-1)%/%32*30
 cbind(head(x,-1),head(y,-1),tail(x,-1),tail(y,-1))
}),1000,ceiling(n/32)*30+20)
many<-make(1000,20);long<-make(1,20000);curved<-make(8,20)
props<-c("orientation","length","curvature")
physical<-feature_spec(props,units="visual_degrees",pixels_per_degree=20)
arc<-feature_spec(props,estimator="arclength",curvature_scale=1,bins=8)
jobs<-list(
 forest=function()feature(forest,props,bins=8),
 many_contours=function()feature(many,props,bins=8),
 one_long_contour=function()feature(long,props,bins=8),
 physical_forest=function()feature(forest,spec=physical),
 arclength=function()feature(curved,spec=arc))
outputs<-lapply(jobs,function(f)f())
timer<-function(f) {
 f();number<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(number))f())[3])
  if(elapsed>=.12||number>=512)break
  number<-number*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(number))f())[3])/number)
 list(median_seconds=median(times),min_seconds=min(times),max_seconds=max(times),iterations=number)
}
rows<-list();allocations<-list()
for(name in names(jobs)) {
 f<-jobs[[name]]
 rows[[name]]<-data.frame(workload=name,method,as.list(timer(f)),
                         package_version=as.character(packageVersion("mlvtoolbox")))
 # Allocation accounting is separate from timing and includes result creation.
 path<-file.path(dest,paste0(name,".Rprofmem"))
 Rprofmem(path,threshold=0);for(i in 1:5)f();Rprofmem(NULL)
 lines<-readLines(path);bytes<-suppressWarnings(as.numeric(sub(" .*","",lines)))
 allocations[[name]]<-data.frame(workload=name,method,iterations=5,
   R_allocated_bytes_per_call=sum(bytes,na.rm=TRUE)/5,R_allocations_per_call=sum(is.finite(bytes))/5)
 cat(name,rows[[name]]$median_seconds,"seconds\n")
}
write.csv(do.call(rbind,rows),file.path(dest,"timings.csv"),row.names=FALSE)
write.csv(do.call(rbind,allocations),file.path(dest,"allocations.csv"),row.names=FALSE)
saveRDS(outputs,file.path(dest,"outputs.rds"))
Rprof(file.path(dest,"forest.Rprof"),interval=.001)
for(i in 1:300)jobs$forest()
Rprof(NULL)
capture.output(summaryRprof(file.path(dest,"forest.Rprof"))$by.total,file=file.path(dest,"forest-profile.txt"))
# Same explicitly cleaned data for the fresh Python comparator.
write.table(do.call(rbind,lapply(seq_along(forest$contours),function(i)cbind(i,forest$contours[[i]]))),
 file.path(dest,"forest.csv"),sep=",",row.names=FALSE,col.names=FALSE)
write.table(c(forest$width,forest$height),file.path(dest,"dimensions.csv"),sep=",",row.names=FALSE,col.names=FALSE)
write.table(as.matrix(outputs$forest[,-1]),file.path(dest,"R-histogram.csv"),sep=",",row.names=FALSE,col.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
