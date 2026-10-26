args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
base<-".agent-work/parity-perf/natural"
read<-function(name)unname(as.matrix(read.csv(file.path(base,name),header=FALSE)))
z<-read("input.csv");size<-as.numeric(read("dimensions.csv"))
d<-mlv_drawing(lapply(unique(z[,1]),function(i)z[z[,1]==i,2:5,drop=FALSE]),size[1],size[2])
g<-compute_contour_properties(d,c("length","orientation","curvature"))
actual_g<-as.matrix(g$segment_properties[,c("length","orientation","curvature")])
actual_h<-as.matrix(feature_table(g,c("orientation","length","curvature"))[,-1])
rows<-list()
for(backend in c("octave","python"))for(stage in c("geometry","histogram")) {
 a<-if(stage=="geometry")actual_g else actual_h;b<-read(paste0(backend,"-",stage,".csv"))
 tolerance<-if(stage=="geometry"&&backend=="python")3e-8 else 1e-8
 rows[[length(rows)+1L]]<-data.frame(backend,stage,n=length(a),max_abs_error=max(abs(a-b)),
  tolerance,status=if(max(abs(a-b))<=tolerance)"match" else "different")
}
write.csv(do.call(rbind,rows),file.path(base,paste0(args[2],"-comparison.csv")),row.names=FALSE)
timer<-function(f) {
 first<-unname(system.time(f())[3]);number<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(number))f())[3]);if(elapsed>=.15||number>=128)break;number<-number*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(number))f())[3])/number)
 c(first_seconds=first,median_seconds=median(times),min_seconds=min(times),max_seconds=max(times),iterations=number)
}
geometry<-function()compute_contour_properties(d,c("length","orientation","curvature"))
histogram<-function()feature_table(g,c("orientation","length","curvature"))
times<-rbind(data.frame(stage="geometry",as.list(timer(geometry))),
             data.frame(stage="histogram",as.list(timer(histogram))))
write.csv(times,file.path(base,paste0(args[2],"-timings.csv")),row.names=FALSE)
print(do.call(rbind,rows));print(times)
