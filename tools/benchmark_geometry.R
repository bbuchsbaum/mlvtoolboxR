args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
timer<-function(f) {
 cold<-unname(system.time(f())[3]);number<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(number))f())[3]);if(elapsed>=.12||number>=1024L)break;number<-number*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(number))f())[3])/number)
 c(cold_seconds=cold,median_seconds=median(times),min_seconds=min(times),max_seconds=max(times),iterations=number)
}
rows<-list()
for(case in list(c(100,1),c(1000,1),c(1,10000))) {
 count<-case[1];segments<-case[2]
 contours<-lapply(seq_len(count),function(i){x<-seq_len(segments);cbind(x,i,x+1,i)})
 drawing<-mlv_drawing(contours,max(segments+2,count+2),count+2)
 f<-function()compute_contour_properties(drawing,c("orientation","length","curvature"))
 result<-f();stopifnot(all(result$segment_properties$length==1),
   all(result$segment_properties$orientation==0),all(result$segment_properties$curvature==0))
 rows[[length(rows)+1L]]<-data.frame(stage="contour_properties",size=count,segments_per_contour=segments,t(timer(f)))
 if(count==1000) {
  Rprof(args[3],interval=.001)
  for(k in 1:10)f()
  Rprof(NULL)
  capture.output(summaryRprof(args[3])$by.total,file=paste0(args[3],".summary.txt"))
 }
}
for(n in c(32L,64L,128L,512L)) {
 width<-n+7L;features<-matrix(FALSE,n,width);r<-round(n/2);c<-round(width/2);features[r,c]<-TRUE
 result<-distance_transform(features);expected<-sqrt((row(features)-r)^2+(col(features)-c)^2)
 stopifnot(max(abs(result$distance-expected))<1e-12,
           all(result$index==r+(c-1)*n))
 rows[[length(rows)+1L]]<-data.frame(stage="distance_transform",size=n,segments_per_contour=0,
    t(timer(function()distance_transform(features))))
}
out<-do.call(rbind,rows);out$source_commit<-system2("git",c("rev-parse","HEAD"),stdout=TRUE)
out$R<-R.version.string
utils::write.csv(out,args[2],row.names=FALSE)
print(out[,1:8],row.names=FALSE)
