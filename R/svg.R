.svg_numbers <- function(text) {
 if(is.na(text)||!nzchar(text))return(numeric())
 token<-regmatches(text,gregexpr("[-+]?(?:[0-9]*\\.[0-9]+|[0-9]+\\.?[0-9]*)(?:[eE][-+]?[0-9]+)?",text,perl=TRUE))[[1]]
 as.numeric(token)
}
.svg_transform <- function(text) {
 out<-diag(3)
 if(is.na(text)||!nzchar(text))return(out)
 terms<-regmatches(text,gregexpr("[A-Za-z]+\\s*\\([^)]*\\)",text,perl=TRUE))[[1]]
 residue<-gsub("[A-Za-z]+\\s*\\([^)]*\\)|[[:space:],]","",text,perl=TRUE)
 if(nzchar(residue))stop("invalid SVG transform",call.=FALSE)
 for(term in terms) {
  name<-sub("\\s*\\(.*","",term);v<-.svg_numbers(term);m<-diag(3)
  if(name=="matrix" && length(v)==6L)m<-matrix(c(v[1],v[2],0,v[3],v[4],0,v[5],v[6],1),3)
  else if(name=="translate" && length(v)%in%c(1L,2L))m[1:2,3]<-c(v[1],if(length(v)>1)v[2] else 0)
  else if(name=="scale" && length(v)%in%c(1L,2L))diag(m)<-c(v[1],if(length(v)>1)v[2] else v[1],1)
  else if(name=="rotate" && length(v)%in%c(1L,3L)) {
   a<-v[1]*pi/180;m[1:2,1:2]<-matrix(c(cos(a),sin(a),-sin(a),cos(a)),2)
   if(length(v)==3L)m[1:2,3]<-v[2:3]-m[1:2,1:2]%*%v[2:3]
  } else if(name=="skewX" && length(v)==1L)m[1,2]<-tan(v*pi/180)
  else if(name=="skewY" && length(v)==1L)m[2,1]<-tan(v*pi/180)
  else stop("unsupported or malformed SVG transform: ",name,call.=FALSE)
  out<-out%*%m
 };out
}
.svg_arc <- function(start,end,rx,ry,phi,large,sweep,step) {
 rx<-abs(rx);ry<-abs(ry)
 if(all(start==end))return(matrix(start,1))
 if(rx==0||ry==0)return(rbind(start,end))
 phi<-phi*pi/180;co<-cos(phi);si<-sin(phi);delta<-(start-end)/2
 xp<-co*delta[1]+si*delta[2];yp<--si*delta[1]+co*delta[2]
 scale<-(xp/rx)^2+(yp/ry)^2
 if(scale>1){rx<-rx*sqrt(scale);ry<-ry*sqrt(scale)}
 den<-rx^2*yp^2+ry^2*xp^2
 factor<-(if(large==sweep)-1 else 1)*sqrt(max(0,(rx^2*ry^2-den)/den))
 cp<-c(factor*rx*yp/ry,-factor*ry*xp/rx)
 center<-c(co*cp[1]-si*cp[2],si*cp[1]+co*cp[2])+(start+end)/2
 u<-c((xp-cp[1])/rx,(yp-cp[2])/ry);v<-c((-xp-cp[1])/rx,(-yp-cp[2])/ry)
 theta<-atan2(u[2],u[1]);dt<-atan2(u[1]*v[2]-u[2]*v[1],sum(u*v))
 if(!sweep&&dt>0)dt<-dt-2*pi
 if(sweep&&dt<0)dt<-dt+2*pi
 n<-max(1,ceiling(abs(dt)*max(rx,ry)/step))
 if(n>1000000L)stop("SVG curve exceeds the point budget; increase curve_step",call.=FALSE)
 t<-seq(theta,theta+dt,length.out=n+1L)
 out<-cbind(center[1]+co*rx*cos(t)-si*ry*sin(t),center[2]+si*rx*cos(t)+co*ry*sin(t))
 out[1,]<-start;out[nrow(out),]<-end;out
}
.svg_path <- function(text,step) {
 if(length(text)!=1L||is.na(text))stop("SVG path is missing d",call.=FALSE)
 pattern<-"[A-Za-z]|[-+]?(?:[0-9]*\\.[0-9]+|[0-9]+\\.?[0-9]*)(?:[eE][-+]?[0-9]+)?"
 tokens<-regmatches(text,gregexpr(pattern,text,perl=TRUE))[[1]]
 pos<-1L;p<-c(0,0);start<-p;current<-matrix(numeric(),0,2);pieces<-list();cmd<-NULL;previous<-NULL;control<-NULL;started<-FALSE
 flush<-function(){if(nrow(current)>1L)pieces[[length(pieces)+1L]]<<-current;current<<-matrix(numeric(),0,2)}
 append<-function(points) {
  if(!nrow(current))current<<-matrix(p,1)
  if(nrow(points)>1L)current<<-rbind(current,points[-1,,drop=FALSE])
  p<<-points[nrow(points),]
 }
 while(pos<=length(tokens)) {
  if(grepl("^[A-Za-z]$",tokens[pos])){cmd<-tokens[pos];pos<-pos+1L}
  if(is.null(cmd))stop("SVG path must begin with a command",call.=FALSE)
  op<-toupper(cmd);relative<-cmd!=op
  if(!started&&op!="M")stop("SVG path must start with M or m",call.=FALSE)
  if(op=="M")started<-TRUE
  if(op=="Z") {
   if(any(p!=start))append(rbind(p,start))
   flush();p<-start;previous<-"Z";control<-NULL;cmd<-NULL;next
  }
  count<-switch(op,M=2L,L=2L,H=1L,V=1L,Q=4L,T=2L,C=6L,S=4L,A=7L,
                stop("unsupported SVG path command: ",cmd,call.=FALSE))
  if(pos+count-1L>length(tokens)||any(grepl("^[A-Za-z]$",tokens[pos:(pos+count-1L)])))
   stop("incomplete SVG path command: ",cmd,call.=FALSE)
  v<-as.numeric(tokens[pos:(pos+count-1L)]);pos<-pos+count
  endpoint<-function(v)if(relative)v+p else v
  if(op=="M") {
   flush();p<-endpoint(v);start<-p;current<-matrix(p,1);cmd<-if(relative)"l" else "L";control<-NULL
  } else if(op %in% c("L","H","V")) {
   end<-if(op=="L")endpoint(v) else if(op=="H")c(if(relative)p[1]+v else v,p[2]) else c(p[1],if(relative)p[2]+v else v)
   append(rbind(p,end));control<-NULL
  } else if(op %in% c("Q","T","C","S")) {
   if(op=="Q"){ctrl<-rbind(endpoint(v[1:2]));end<-endpoint(v[3:4])}
   else if(op=="T"){ctrl<-matrix(if(previous %in% c("Q","T"))2*p-control else p,1);end<-endpoint(v)}
   else if(op=="C"){ctrl<-rbind(endpoint(v[1:2]),endpoint(v[3:4]));end<-endpoint(v[5:6])}
   else {ctrl<-rbind(if(previous %in% c("C","S"))2*p-control else p,endpoint(v[1:2]));end<-endpoint(v[3:4])}
   polygon<-rbind(p,ctrl,end);length<-sum(sqrt(rowSums(diff(polygon)^2)));n<-max(1,ceiling(length/step))
   if(n>1000000L)stop("SVG curve exceeds the point budget; increase curve_step",call.=FALSE)
   t<-seq(0,1,length.out=n+1L)
   if(nrow(ctrl)==1L)points<-outer((1-t)^2,p)+outer(2*(1-t)*t,ctrl[1,])+outer(t^2,end)
   else points<-outer((1-t)^3,p)+outer(3*(1-t)^2*t,ctrl[1,])+outer(3*(1-t)*t^2,ctrl[2,])+outer(t^3,end)
   append(points);control<-ctrl[nrow(ctrl),]
  } else if(op=="A") {
   if(!all(v[4:5] %in% c(0,1)))stop("invalid SVG arc flags",call.=FALSE)
   append(.svg_arc(p,endpoint(v[6:7]),v[1],v[2],v[3],v[4]==1,v[5]==1,step));control<-NULL
  }
  previous<-op
 }
 flush();pieces
}
#' Import vector contours from SVG
#'
#' Supports lines, polylines, polygons, rectangles, circles/ellipses, Bezier
#' paths, elliptical arcs and affine transforms. The viewBox origin maps to
#' (1,1) in image coordinates. Curves are sampled according to control-polygon
#' length; curve_step is a sampling parameter, not a certified Hausdorff bound.
#' Text, embedded images, clipping, masks and paint effects are not rendered.
#' @param file SVG file.
#' @param width,height Optional output dimensions.
#' @param curve_step Positive curve sampling step.
#' @param strict Error on unsupported visible content; otherwise warn.
#' @return A drawing with import provenance.
#' @export
read_mlv_svg <- function(file,width=NULL,height=NULL,curve_step=1,strict=TRUE) {
 if(!requireNamespace("xml2",quietly=TRUE))stop("install optional package 'xml2'",call.=FALSE)
 if(length(curve_step)!=1L||!is.finite(curve_step)||curve_step<=0)stop("invalid curve_step",call.=FALSE)
 doc<-xml2::read_xml(file,options="NONET");root<-xml2::xml_root(doc)
 attrnum<-function(node,name,default=0) {
  s<-xml2::xml_attr(node,name)
  if(is.na(s))return(default)
  if(!grepl("^[-+0-9.eE[:space:]]+(px)?$",s))stop("unsupported SVG unit in ",name,call.=FALSE)
  .svg_numbers(s)[1]
 }
 vb<-.svg_numbers(xml2::xml_attr(root,"viewBox"))
 if(length(vb)&&length(vb)!=4L)stop("viewBox needs four values",call.=FALSE)
 if(is.null(width))width<-ceiling(if(!is.na(xml2::xml_attr(root,"width"))&&!grepl("%",xml2::xml_attr(root,"width"),fixed=TRUE))attrnum(root,"width",NA_real_) else if(length(vb))vb[3] else NA_real_)
 if(is.null(height))height<-ceiling(if(!is.na(xml2::xml_attr(root,"height"))&&!grepl("%",xml2::xml_attr(root,"height"),fixed=TRUE))attrnum(root,"height",NA_real_) else if(length(vb))vb[4] else NA_real_)
 empty<-mlv_drawing(list(),width,height)
 transform<-diag(3);transform[1:2,3]<-1
 if(length(vb)) {
  if(any(vb[3:4]<=0))stop("viewBox dimensions must be positive",call.=FALSE)
  scale<-c(width,height)/vb[3:4]
  alignment<-c(0,0)
  aspect<-xml2::xml_attr(root,"preserveAspectRatio")
  if(is.na(aspect))aspect<-"xMidYMid meet"
  if(aspect!="none") {
   if(!grepl("^x(Min|Mid|Max)Y(Min|Mid|Max)( (meet|slice))?$",aspect))
    stop("unsupported preserveAspectRatio",call.=FALSE)
   uniform<-if(grepl("slice$",aspect))max(scale) else min(scale)
   extra<-c(width,height)-vb[3:4]*uniform
   ax<-if(grepl("^xMin",aspect))0 else if(grepl("^xMax",aspect))1 else .5
   ay<-if(grepl("YMin",aspect))0 else if(grepl("YMax",aspect))1 else .5
   alignment<-extra*c(ax,ay);scale<-rep(uniform,2)
  }
  transform[1,1]<-scale[1];transform[2,2]<-scale[2]
  transform[1:2,3]<-1-vb[1:2]*scale+alignment
 }
 pieces<-list()
 unsupported<-function(message){if(strict)stop(message,call.=FALSE) else warning(message,call.=FALSE)}
 visit<-function(node,parent) {
  name<-xml2::xml_name(node)
  if(name %in% c("defs","metadata","title","desc","style"))return(invisible(NULL))
  style<-xml2::xml_attr(node,"style")
  if(identical(xml2::xml_attr(node,"display"),"none") ||
     (!is.na(style)&&grepl("display\\s*:\\s*none",style,perl=TRUE)))return(invisible(NULL))
  if(name=="svg"&&xml2::xml_path(node)!=xml2::xml_path(root))
   unsupported("nested SVG viewports are unsupported")
  for(a in c("clip-path","mask","filter"))if(!is.na(xml2::xml_attr(node,a)))unsupported(paste("unsupported SVG effect:",a))
  tr<-parent%*%.svg_transform(xml2::xml_attr(node,"transform"));paths<-list()
  if(name %in% c("svg","g")) {
   for(child in xml2::xml_children(node))visit(child,tr)
   return(invisible(NULL))
  } else if(name=="line")paths<-list(rbind(c(attrnum(node,"x1"),attrnum(node,"y1")),c(attrnum(node,"x2"),attrnum(node,"y2"))))
  else if(name %in% c("polyline","polygon")) {
   p<-.svg_numbers(xml2::xml_attr(node,"points"))
   if(length(p)%%2L)stop("unpaired SVG point",call.=FALSE)
   p<-matrix(p,ncol=2,byrow=TRUE)
   if(name=="polygon"&&nrow(p)&&any(p[1,]!=p[nrow(p),]))p<-rbind(p,p[1,])
   paths<-list(p)
  } else if(name=="rect") {
   if(attrnum(node,"rx")!=0||attrnum(node,"ry")!=0)unsupported("rounded SVG rectangles are unsupported")
   a<-attrnum(node,"x");b<-attrnum(node,"y");w<-attrnum(node,"width");h<-attrnum(node,"height")
   if(w<0||h<0)stop("negative rectangle dimensions",call.=FALSE)
   paths<-list(rbind(c(a,b),c(a+w,b),c(a+w,b+h),c(a,b+h),c(a,b)))
  } else if(name %in% c("circle","ellipse")) {
   rx<-if(name=="circle")attrnum(node,"r") else attrnum(node,"rx")
   ry<-if(name=="circle")rx else attrnum(node,"ry")
   if(rx<0||ry<0)stop("negative ellipse radius",call.=FALSE)
   n<-max(12,ceiling(2*pi*max(rx,ry)/curve_step))
   if(n>1000000L)stop("SVG curve exceeds point budget",call.=FALSE)
   a<-seq(0,2*pi,length.out=n+1L)
   paths<-list(cbind(attrnum(node,"cx")+rx*cos(a),attrnum(node,"cy")+ry*sin(a)))
  } else if(name=="path")paths<-.svg_path(xml2::xml_attr(node,"d"),curve_step)
  else {unsupported(paste("unsupported SVG element:",name));return(invisible(NULL))}
  for(p in paths)if(nrow(p)>1L) {
   xy<-cbind(p,1)%*%t(tr);xy<-xy[,1:2,drop=FALSE]
   z<-cbind(xy[-nrow(xy),,drop=FALSE],xy[-1,,drop=FALSE])
   pieces[[length(pieces)+1L]]<<-z
  }
 }
 visit(root,transform)
 mlv_drawing(pieces,empty$width,empty$height,zero_length="drop",
   source=list(image=file,method="svg",curve_step=curve_step,viewBox=vb,origin=c(1,1)))
}
