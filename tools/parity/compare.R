# Differential receipt; intentionally reports mismatches instead of hiding them.
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
base<-if(length(args)>=3)args[3] else "inst/fixtures/parity"
read<-function(name)unname(as.matrix(read.csv(file.path(base,name),header=FALSE)))
cases<-read.csv(file.path(base,"cases.csv"));rows<-list()
record<-function(case,backend,stage,actual,reference,atol=1e-10,rtol=1e-12,circular=FALSE) {
 shape_ok<-identical(dim(actual),dim(reference))||is.null(dim(actual))&&length(actual)==length(reference)
 a<-as.numeric(actual);b<-as.numeric(reference)
 if(!shape_ok||length(a)!=length(b)) {
  rows[[length(rows)+1L]]<<-data.frame(case,backend,stage,n=length(a),max_abs_error=Inf,
   tolerance=atol,disagreements=NA_integer_,status="shape_mismatch")
  return(invisible(NULL))
 }
 error<-abs(a-b);if(circular)error<-pmin(error,abs(360-error))
 tol<-atol+rtol*max(c(abs(b[is.finite(b)]),0))
 bad<-is.na(error)|error>tol
 rows[[length(rows)+1L]]<<-data.frame(case,backend,stage,n=length(a),
   max_abs_error=max(c(error,0),na.rm=TRUE),tolerance=tol,disagreements=sum(bad),
   status=if(any(bad))"different" else "match")
}
for(i in seq_len(nrow(cases))) {
 name<-cases$name[i]
 if(cases$family[i]=="geometry") {
  z<-read(paste0(name,".input.csv"))
  d<-mlv_drawing(lapply(unique(z[,1]),function(id)z[z[,1]==id,2:5,drop=FALSE]),256,192)
  st<-as.data.frame(d);actual<-as.matrix(st[,c("length","orientation","curvature")])
  for(backend in c("octave","python")) {
   ref<-read(paste0(name,".",backend,".geometry.csv"))
   for(k in 1:3)record(name,backend,colnames(actual)[k],actual[,k],ref[,k],
      atol=if(k==3&&backend=="python")3e-8 else 1e-10,circular=k==2)
   for(bins in c(5,8,13))for(prop in c("orientation","length","curvature"))for(norm in c(FALSE,TRUE)) {
    field<-paste0(if(norm)"norm_" else "",prop,bins)
    prefix<-if(backend=="python"&&prop=="length")"adapted_" else ""
    ref<-read(paste0(name,".",backend,".",prefix,field,".csv"))
    h<-contour_statistics(d,prop,bins,by="contour",normalize=norm)
    record(name,backend,paste0(prefix,field),matrix(h$value,ncol=bins,byrow=TRUE),ref,atol=1e-9)
   }
   record(name,backend,"rotate37_ccw",do.call(rbind,rotate_drawing(d,37)$contours),
          read(paste0(name,".",backend,".rotate37.csv")))
  }
 } else if(cases$family[i]=="mat") {
  m<-compute_mat(read(paste0(name,".input.csv"))>0,trace_branches=FALSE)
  D<-read(paste0(name,".octave.distance.csv"));I<-read(paste0(name,".octave.index.csv"))
  A<-read(paste0(name,".octave.aof.csv"));S<-read(paste0(name,".octave.skeleton.csv"))
  record(name,"octave","mat_distance",m$distance_map,D,2e-6)
  record(name,"octave","mat_index",m$nearest_boundary,I,0,0)
  record(name,"octave","mat_aof",m$AOF,A)
  record(name,"octave","mat_skeleton",m$skeleton,S,0,0)
  record(name,"octave","aof_common_fields",average_outward_flux(D,I),A)
  record(name,"octave","morphology_common_aof",skeletonize(A>(2/pi*sin(28*pi/180))*60),S,0,0)
  h<-nrow(D);r<-row(D);c<-col(D);a<-m$nearest_boundary-1;b<-I-1
  da<-(a%%h+1-r)^2+(a%/%h+1-c)^2;db<-(b%%h+1-r)^2+(b%/%h+1-c)^2
  record(name,"octave","nearest_label_squared_distance",da,db)
 } else if(cases$family[i]=="edt_aof") {
  mask<-read(paste0(name,".input.csv"))>0;dt<-distance_transform(mask)
  D<-read(paste0(name,".octave.distance.csv"));I<-read(paste0(name,".octave.index.csv"))
  for(backend in c("octave","python")) {
   record(name,backend,"edt_distance",dt$distance,read(paste0(name,".",backend,".distance.csv")),if(backend=="octave")2e-6 else 1e-12)
   for(eps in c(1,2))record(name,backend,paste0("aof_common_fields_e",eps),
     average_outward_flux(D,I,epsilon=eps),read(paste0(name,".",backend,".aof",eps,".csv")))
  }
 } else if(cases$family[i]=="branch") {
  z<-read(paste0(name,".branch.csv"));branch<-data.frame(x=z[,1],y=z[,2],radius=z[,3])
  ref<-read(paste0(name,".octave.scores.csv"))
  for(k in 1:3) {
   prop<-c("parallelism","separation","mirror")[k]
   record(name,"octave",paste0("branch_",prop),score_branch(branch,prop),ref[,k])
  }
 }
}
out<-do.call(rbind,rows)
explain<-function(i) {
 r<-out[i,]
 if(r$status=="match")return(if(grepl("^adapted_",r$stage))"matched_with_python_shape_adapter" else if(grepl("^branch_",r$stage))"matched_with_fixed_movmean_adapter" else "matched")
 if(r$status!="different")return("UNEXPLAINED")
 if(r$backend=="octave") {
  if(r$stage %in% c("mat_index","mat_aof","mat_skeleton")) {
   witness<-out[out$case==r$case & out$backend=="octave" & out$stage %in%
                  c("nearest_label_squared_distance","aof_common_fields","morphology_common_aof"),]
   if(nrow(witness)==3&&all(witness$status=="match"))return("equidistant_label_policy")
  }
  if(r$stage=="branch_mirror"&&r$case %in% c("curved","terminal_bend"))return("MLV-M016_terminal_fit")
  if(grepl("orientation13$",r$stage)&&grepl("^integer",r$case))return("MLV-M017_axial_boundary")
 } else if(r$backend=="python") {
  if(r$stage=="rotate37_ccw")return("MLV-P013_trigonometric_rounding")
  if(grepl("^aof_common_fields",r$stage))return("MLV-P012_sampling")
  if(r$case %in% c("subpixel","fractional")&&grepl("orientation|curvature",r$stage))return("MLV-P001_coordinate_truncation")
  if(grepl("orientation13$",r$stage)&&grepl("^integer",r$case))return("MLV-P015_axial_boundary")
 }
 "UNEXPLAINED"
}
out$explanation<-vapply(seq_len(nrow(out)),explain,character(1))
dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)

write.csv(out,args[2],row.names=FALSE)
print(with(out,table(backend,status)))
print(out[out$status!="match"&!grepl("orientation|curvature",out$stage),],row.names=FALSE)

if(any(out$explanation=="UNEXPLAINED"))stop("unexplained parity discrepancy; see receipt")
