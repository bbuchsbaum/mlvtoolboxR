library(mlvtoolbox)
source("tools/voronoi/fixtures.R")
source("tools/voronoi/preprocess.R")
destination<-"docs/evidence/voronoi-spike"
executable<-normalizePath(".agent-work/voronoi-spike/prototype")
cases<-voronoi_fixtures();results<-list();pairs<-list()
for(name in names(cases)) {
 case<-cases[[name]];x<-case$drawing
 prepared<-tryCatch(quantize_segments(x),error=function(e)e)
 if(inherits(prepared,"error")) {
  stopifnot(!is.null(case$error),grepl(case$error,conditionMessage(prepared)))
  writeLines(conditionMessage(prepared),file.path(destination,paste0(name,"-rejected.txt")))
  results[[name]]<-data.frame(case=name,status="rejected",reason=conditionMessage(prepared))
  next
 }
 stopifnot(is.null(case$error),prepared$endpoint_error<=prepared$maximum_error_bound+1e-12)
 input<-file.path(destination,paste0(name,"-input.txt"));rawfile<-file.path(destination,paste0(name,"-raw.csv"))
 logfile<-file.path(destination,paste0(name,"-native.log"))
 writeLines(c(as.character(nrow(prepared$segments)),apply(prepared$segments,1,function(z)paste(format(z,scientific=FALSE,trim=TRUE),collapse=" "))),input)
 write.csv(data.frame(piece=seq_len(nrow(prepared$segments)),source_segment=prepared$source_segment),
           file.path(destination,paste0(name,"-owners.csv")),row.names=FALSE)
 started<-proc.time()[['elapsed']]
 status<-system2("/usr/bin/time",c("-l",shQuote(executable),"1000"),stdin=input,stdout=rawfile,stderr=logfile)
 call_seconds<-proc.time()[['elapsed']]-started;stopifnot(status==0)
 raw<-read.csv(rawfile);raw$original_a<-prepared$source_segment[raw$site_a]
 raw$original_b<-prepared$source_segment[raw$site_b]
 eligible<-which(raw$primary==1 & raw$linear==1 & raw$finite==1 & raw$original_a!=raw$original_b)
 selected<-lapply(eligible,function(i) {
  clipped<-clip_probe_edge(as.numeric(raw[i,c("x1","y1","x2","y2")])/prepared$scale,x$width,x$height)
  if(is.null(clipped))return(NULL)
  data.frame(original_a=raw$original_a[i],original_b=raw$original_b[i],x1=clipped[1],y1=clipped[2],x2=clipped[3],y2=clipped[4])
 })
 selected<-Filter(Negate(is.null),selected)
 selected<-if(length(selected))do.call(rbind,selected) else data.frame(original_a=integer(),original_b=integer(),x1=double(),y1=double(),x2=double(),y2=double())
 write.csv(selected,file.path(destination,paste0(name,"-selected.csv")),row.names=FALSE)
 error<-violation<-numeric(nrow(selected))
 for(i in seq_len(nrow(selected))) {
  z<-as.numeric(selected[i,c("x1","y1","x2","y2")]);p<-(z[1:2]+z[3:4])/2
  distances<-apply(prepared$quantized_original,1,function(s)probe_segment_distance(p,s))
  a<-distances[selected$original_a[i]];b<-distances[selected$original_b[i]]
  error[i]<-abs(a-b);violation[i]<-max(a,b)-min(distances)
 }
 stopifnot(all(error<1e-8),all(violation<1e-8))
 axis_error<-axis_length<-NA_real_
 if(!is.null(case$axis_y)) {
  axis_error<-max(abs(c(selected$y1,selected$y2)-case$axis_y))
  axis_length<-sum(sqrt((selected$x2-selected$x1)^2+(selected$y2-selected$y1)^2))
  stopifnot(axis_error<=prepared$maximum_error_bound+1e-10,abs(axis_length-case$axis_length)<1e-8)
 }
 primary<-raw[raw$primary==1&raw$original_a!=raw$original_b,,drop=FALSE]
 pairs[[name]]<-sort(unique(vapply(seq_len(nrow(primary)),function(i)paste(sort(c(primary$original_a[i],primary$original_b[i])),collapse=":"),character(1))))
 timing<-readLines(logfile)
 native_seconds<-as.numeric(sub(".*construct_seconds=([^ ]+).*","\\1",timing[grepl("construct_seconds=",timing)][1]))
 rss_line<-timing[grepl("maximum resident set size",timing)]
 rss<-if(length(rss_line))as.numeric(trimws(sub("maximum resident set size.*","",rss_line[1]))) else NA_real_
 raster_times<-numeric(25)
 for(k in seq_along(raster_times)) {
  t<-proc.time()[['elapsed']];mat<-compute_mat(x,min_component=1,min_branch_points=2)
  raster_times[k]<-proc.time()[['elapsed']]-t
 }
 corridor<-NA_integer_;coincident_masks<-NA
 if(!is.null(case$axis_y)) {
  lines<-prepared$original;low<-min(lines[,c(2,4)]);high<-max(lines[,c(2,4)])
  xy<-arrayInd(which(mat$skeleton),dim(mat$skeleton))
  corridor<-sum(xy[,1]>low & xy[,1]<high & xy[,2]>=7 & xy[,2]<=23)
  coincident_masks<-identical(rasterize(subset_contours(x,1)),rasterize(subset_contours(x,2)))
 }
 write.csv(data.frame(row=which(mat$skeleton,arr.ind=TRUE)[,1],column=which(mat$skeleton,arr.ind=TRUE)[,2]),file.path(destination,paste0(name,"-raster-skeleton.csv")),row.names=FALSE)
 results[[name]]<-data.frame(case=name,status="measured",reason="",input_pieces=nrow(prepared$segments),
  endpoint_error=prepared$endpoint_error,primary_edges=nrow(primary),selected_linear_edges=nrow(selected),
  unsupported_curved=sum(primary$linear==0),unsupported_infinite=sum(primary$finite==0),
  max_equidistance_error=if(length(error))max(error) else NA_real_,max_nearest_violation=if(length(violation))max(violation) else NA_real_,
  axis_error=axis_error,axis_length=axis_length,native_construct_seconds=native_seconds,
  native_call_seconds=call_seconds,native_peak_rss_bytes=rss,raster_mean_seconds=mean(raster_times),
  raster_ink_pixels=sum(mat$mask),coincident_line_masks=coincident_masks,
  raster_skeleton_pixels=sum(mat$skeleton),raster_corridor_pixels=corridor)
}
for(name in c("shifted_subquantum","shifted_quantum","shifted_tenth","narrow"))stopifnot(identical(pairs[[name]],pairs$parallel))
stopifnot(results$narrow$raster_corridor_pixels==0,results$narrow$coincident_line_masks,
          results$curved$unsupported_curved>0)
columns<-unique(unlist(lapply(results,names)))
results<-lapply(results,function(x){for(n in setdiff(columns,names(x)))x[[n]]<-NA;x[,columns]})
summary<-do.call(rbind,results);rownames(summary)<-NULL
write.csv(summary,file.path(destination,"results.csv"),row.names=FALSE)
jsonlite::write_json(pairs,file.path(destination,"source-pairs.json"),pretty=TRUE,auto_unbox=TRUE)
print(summary)
