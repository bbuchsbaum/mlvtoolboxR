# Run after installing an optimized package into the supplied library:
# Rscript tools/benchmark.R /path/to/library docs/evidence/benchmark-r.csv
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2L)stop("supply package library and CSV destination")
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
timer<-function(f) {
 cold<-unname(system.time(f())[3]);n<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(n))f())[3]);if(elapsed>=.12||n>=1024L)break;n<-n*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(n))f())[3])/n)
 c(cold_seconds=cold,median_seconds=median(times),min_seconds=min(times),
   max_seconds=max(times),iterations=n)
}
rows<-list()
for(n in c(32L,64L,128L,512L)) {
 feature<-matrix(FALSE,n,n);feature[round(n/2),round(n/2)]<-TRUE
 dt<-distance_transform(feature)
 # Warm kernel timing, no input generation/model loading inside the timed call.
 reference<-file.path(dirname(args[2]),paste0("aof-octave-",n,".csv"))
 if(file.exists(reference)) {
  expected<-unname(as.matrix(utils::read.csv(reference,header=FALSE)))
  stopifnot(max(abs(average_outward_flux(dt$distance,dt$index)-expected))<1e-10)
 }
 a<-timer(function()average_outward_flux(dt$distance,dt$index))
 rows[[length(rows)+1L]]<-data.frame(stage="aof",size=n,t(a))
 ink<-matrix(FALSE,n,n);ink[3:(n-2),round(n*.25)]<-TRUE;ink[3:(n-2),round(n*.75)]<-TRUE
 reference<-file.path(dirname(args[2]),paste0("mat-aof-octave-",n,".csv"))
 if(file.exists(reference)) {
  actual<-compute_mat(ink,trace_branches=FALSE)
  f<-function(prefix)unname(as.matrix(utils::read.csv(file.path(dirname(args[2]),paste0(prefix,n,".csv")),header=FALSE)))
  stopifnot(max(abs(actual$AOF-f("mat-aof-octave-")))<1e-10,
            max(abs(actual$distance_map-f("mat-distance-octave-")))<2e-6,
            identical(actual$skeleton,f("mat-skeleton-octave-")>0))
 }
 b<-timer(function()compute_mat(ink,trace_branches=FALSE))
 rows[[length(rows)+1L]]<-data.frame(stage="mat_core",size=n,t(b))
}
for(n in c(100L,1000L,10000L)) {
 contours<-lapply(seq_len(n),function(i)matrix(c(i*3,1,i*3+1,2),1))
 drawing<-mlv_drawing(contours,n*3+2L,5L)
 a<-timer(function()detect_junctions(drawing))
 rows[[length(rows)+1L]]<-data.frame(stage="junction_sparse",size=n,t(a))
}
result<-do.call(rbind,rows)
result$R<-R.version.string
result$machine<-Sys.info()[["machine"]]
result$source_commit<-system2("git",c("rev-parse","HEAD"),stdout=TRUE)
utils::write.csv(result,args[2],row.names=FALSE)
print(result[,1:7],row.names=FALSE)
