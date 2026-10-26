#' Render a drawing with property or model-score colors
#' @param x A drawing.
#' @param property Optional segment or mapped medial-axis property name.
#' @param values Optional numeric score for each contour, for example model predictions.
#' @param width,height Output dimensions.
#' @param line_width Stroke width.
#' @param palette Color vector used from low to high values.
#' @return RGB array with white background. Missing scores are gray.
#' @export
render_drawing <- function(x,property=NULL,values=NULL,width=x$width,height=x$height,
                           line_width=1,palette=grDevices::hcl.colors(256,"Viridis")) {
 .check_drawing(x)
 mask<-rasterize(x,width,height,line_width)
 out<-array(1,c(nrow(mask),ncol(mask),3))
 if(is.null(property)&&is.null(values)){for(k in 1:3)out[,,k][mask]<-0;return(out)}
 if(!is.null(property)&&property=="junctions") {
  for(k in 1:3)out[,,k][mask]<-0
  j<-.measurement_junctions(x,x$feature_spec)
  for(i in seq_len(nrow(j))) {
   cx<-j$x[i]*width/x$width;cy<-j$y[i]*height/x$height
   if(cx+2<1||cx-2>width||cy+2<1||cy-2>height)next
   rr<-seq.int(max(1,ceiling(cy-2)),min(height,floor(cy+2)))
   cc<-seq.int(max(1,ceiling(cx-2)),min(width,floor(cx+2)))
   rr<-rr[rr>=1&rr<=height];cc<-cc[cc>=1&cc<=width]
   if(length(rr)&&length(cc))for(r in rr)for(c in cc)if((r-cy)^2+(c-cx)^2<=4)out[r,c,]<-c(1,0,0)
  }
  return(out)
 }
 if(!is.null(property)&&!is.null(x$mat_images[[property]])) {
  .checked_mat_measurements(x,x$feature_spec)
  image<-x$mat_images[[property]]
  rows<-pmax(1L,pmin(x$height,round(seq(1,x$height,length.out=height))))
  cols<-pmax(1L,pmin(x$width,round(seq(1,x$width,length.out=width))))
  image<-image[rows,cols,drop=FALSE];valid<-is.finite(image)
  index<-pmax(1L,pmin(length(palette),1L+as.integer(image*(length(palette)-1L))))
  colors<-array(grDevices::col2rgb(palette)/255,c(3,length(palette)))
  for(k in 1:3) {
   layer<-out[,,k];layer[mask&!valid]<-.7
   layer[valid]<-colors[k,index[valid]];out[,,k]<-layer
  }
  return(out)
 }
 st<-.measurement_table(x,x$feature_spec)
 if(!is.null(values)) {
  if(!is.numeric(values)||length(values)!=length(x$ids))stop("values must provide one numeric score per contour",call.=FALSE)
  v<-values[st$contour_index]
 } else if(property %in% c("length","orientation","curvature"))v<-st[[property]]
 else if(!is.null(x$mat_properties[[property]])) {
  .checked_mat_measurements(x,x$feature_spec)
  v<-vapply(x$mat_properties[[property]],function(z)if(any(is.finite(z)))mean(z[is.finite(z)]) else NA_real_,double(1))[st$contour_index]
 } else stop("unknown or unmapped render property",call.=FALSE)
 finite<-is.finite(v)
 colors<-rep("gray70",length(v))
 if(any(finite)) {
  bounds<-range(v[finite]);scaled<-if(diff(bounds)==0)rep(.5,sum(finite)) else (v[finite]-bounds[1])/diff(bounds)
  colors[finite]<-palette[1+pmin(length(palette)-1L,as.integer(scaled*(length(palette)-1L)))]
 }
 rgb<-grDevices::col2rgb(colors)/255
 scale<-rep(c(width/x$width,height/x$height),2)
 for(i in seq_len(nrow(st))) {
  seg<-as.matrix(st[i,c("x1","y1","x2","y2")]);seg<-sweep(seg,2,scale,"*")
  index<-cpp_raster_indices(seg,height,width,line_width)
  for(k in 1:3){layer<-out[,,k];layer[index]<-rgb[k,i];out[,,k]<-layer}
 }
 out
}
