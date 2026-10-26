quantize_segments <- function(x,scale=1000) {
 stopifnot(length(scale)==1,is.finite(scale),scale>0)
 original<-as.matrix(as.data.frame(x,spec=NULL)[,c("x1","y1","x2","y2")])
 if(any(abs(original*scale)>.Machine$integer.max))stop("quantized coordinates exceed int32 bounds")
 graph<-contour_graph(x,tolerance=1e-8,gap_max=0,parallel_distance=0)
 edge<-graph$edges
 pieces<-t(vapply(seq_len(nrow(edge)),function(i) {
  z<-original[edge$source_segment[i],];u<-z[3:4]-z[1:2]
  c(z[1:2]+edge$t0[i]*u,z[1:2]+edge$t1[i]*u)
 },numeric(4)))
 rounded<-round(pieces*scale)
 if(any(rounded[,1]==rounded[,3]&rounded[,2]==rounded[,4]))stop("quantization collapses a span; topology changed")
 quantized<-mlv_drawing(lapply(seq_len(nrow(rounded)),function(i)rounded[i,,drop=FALSE]/scale),x$width,x$height)
 after<-contour_graph(quantized,tolerance=1e-8,gap_max=0,parallel_distance=0)
 signature<-function(g,source_map=NULL) {
  e<-g$edges;source<-if(is.null(source_map))e$source_segment else source_map[e$source_segment]
  sort(vapply(g$nodes$node_id,function(i)paste(sort(source[e$from==i|e$to==i]),collapse=":"),character(1)))
 }
 if(!identical(signature(graph),signature(after,edge$source_segment)))stop("quantization changes source-incidence topology")
 # Intersections introduced by rounding would require a fresh split. Refuse
 # them rather than passing illegal intersecting inputs to Boost.
 if(nrow(after$edges)!=nrow(rounded))stop("quantization changes split topology")
 list(segments=rounded,source_segment=edge$source_segment,original=original,
  quantized_original=round(original*scale)/scale,scale=scale,
  endpoint_error=max(sqrt(rowSums((rounded[,1:2,drop=FALSE]/scale-pieces[,1:2,drop=FALSE])^2)),
                     sqrt(rowSums((rounded[,3:4,drop=FALSE]/scale-pieces[,3:4,drop=FALSE])^2))),
  maximum_error_bound=sqrt(2)/(2*scale),graph_before=graph_features(graph),graph_after=graph_features(after))
}

clip_probe_edge <- function(z,width,height) {
 a<-z[1:2];v<-z[3:4]-a;lo<-0;hi<-1
 for(k in 1:2) {
  bound<-if(k==1)width else height
  if(v[k]==0) {if(a[k]<1||a[k]>bound)return(NULL)} else {
   interval<-sort(c((1-a[k])/v[k],(bound-a[k])/v[k]))
   lo<-max(lo,interval[1]);hi<-min(hi,interval[2])
  }
 }
 if(hi<=lo)return(NULL)
 c(a+lo*v,a+hi*v)
}
probe_segment_distance <- function(p,z) {
 a<-z[1:2];v<-z[3:4]-a;t<-max(0,min(1,sum((p-a)*v)/sum(v*v)))
 sqrt(sum((p-a-t*v)^2))
}
