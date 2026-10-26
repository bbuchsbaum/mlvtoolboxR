.check_mask <- function(mask) {
 if(!is.matrix(mask) || !(is.logical(mask)||is.numeric(mask)) || any(dim(mask)<1L))
  stop("mask must be a nonempty binary matrix; TRUE means ink",call.=FALSE)
 if(prod(as.double(dim(mask)))>.Machine$integer.max)
  stop("image exceeds the supported integer index capacity",call.=FALSE)
 if(is.logical(mask)) {
  if(anyNA(mask))stop("mask must be a nonempty binary matrix; TRUE means ink",call.=FALSE)
  return(unname(mask))
 }
 if(any(!is.finite(mask)) || any(!mask %in% c(0,1)))
  stop("mask must be a nonempty binary matrix; TRUE means ink",call.=FALSE)
 matrix(as.logical(mask),nrow(mask),ncol(mask))
}
#' Exact Euclidean distance to the nearest feature pixel
#' @param features Binary matrix; TRUE marks a feature.
#' @return List with distance and nearest one-based column-major index matrices.
#'   Empty feature sets produce Inf distances and NA indices. Ties select
#'   the smallest linear index.
#' @export
distance_transform <- function(features) cpp_distance_transform(.check_mask(features))

#' Thin a binary image using Pratt morphology
#' @param mask Binary ink matrix.
#' @param min_component Minimum eight-connected component size.
#' @return Logical skeleton, with backend metadata described in the package docs.
#' @export
skeletonize <- function(mask,min_component=1L) {
 mask <- .check_mask(mask)
 if(length(min_component)!=1L || !is.finite(min_component) ||
    min_component<1 || min_component!=floor(min_component) ||
    min_component>.Machine$integer.max) stop("invalid min_component",call.=FALSE)
 result <- cpp_skeletonize(mask)
 if(min_component>1L)result <- cpp_remove_small(result,as.integer(min_component))
 result
}

.raster_indices <- function(x,i,width=x$width,height=x$height,line_width=1) {
 z <- x$contours[[i]]
 scale <- rep(c(width/x$width,height/x$height),2L)
 cpp_raster_indices(sweep(z,2L,scale,"*"),height,width,line_width)
}
#' Render a drawing into an ink mask
#'
#' Uses deterministic clipped Bresenham lines for width one. Wider strokes
#' include pixels whose centers lie within half the stroke width of a segment.
#' This rasterizer is explicit and is not claimed equivalent to insertShape.
#' @param x A drawing.
#' @param width,height Output dimensions.
#' @param line_width Positive stroke width in output pixels.
#' @return Logical matrix in row/column order, TRUE for ink.
#' @export
rasterize <- function(x,width=x$width,height=x$height,line_width=1) {
 .check_drawing(x)
 tmp <- mlv_drawing(list(),width,height); width<-tmp$width; height<-tmp$height
 if(as.double(width)*height>.Machine$integer.max)stop("image exceeds index capacity",call.=FALSE)
 if(length(line_width)!=1L || !is.finite(line_width)||line_width<=0)
  stop("line_width must be positive",call.=FALSE)
 mask <- matrix(FALSE,height,width)
 for(i in seq_along(x$contours))mask[.raster_indices(x,i,width,height,line_width)]<-TRUE
 mask
}
#' Read an image as a numeric array
#' @param file PNG or JPEG file.
#' @return Numeric matrix or RGB(A) array in the range zero to one.
#' @export
read_mlv_image <- function(file) {
 ext <- tolower(tools::file_ext(file))
 package <- if(ext=="png")"png" else if(ext %in% c("jpg","jpeg"))"jpeg" else stop("supported formats are PNG and JPEG",call.=FALSE)
 if(!requireNamespace(package,quietly=TRUE))stop("install optional package '",package,"'",call.=FALSE)
 if(package=="png")png::readPNG(file) else jpeg::readJPEG(file)
}
#' Export a drawing or image to PNG
#' @param x Drawing, binary ink mask, or numeric image.
#' @param file Destination PNG path.
#' @param ... Arguments to rasterize for drawings.
#' @return The destination path invisibly.
#' @export
write_mlv_png <- function(x,file,...) {
 if(!requireNamespace("png",quietly=TRUE))stop("install optional package 'png'",call.=FALSE)
 img <- if(inherits(x,"mlv_drawing"))1*rasterize(x,...) else x
 if(is.logical(img) || inherits(x,"mlv_drawing"))img <- 1-img
 if(!(is.matrix(img)||length(dim(img))==3L) || !is.numeric(img) ||
    any(!is.finite(img))||any(img<0|img>1))stop("image values must be in [0,1]",call.=FALSE)
 png::writePNG(img,file)
 invisible(file)
}
#' Trace contours from an ink mask or grayscale line drawing
#'
#' The graph uses eight-neighbor links but suppresses redundant diagonal
#' shortcuts across an existing orthogonal link. Junctions split paths and
#' closed loops retain their closure. This is an explicit graph method,
#' rather than an unqualified emulation of Kovesi edge linking.
#' @param image Binary ink mask, grayscale matrix, RGB array, or image filename.
#' @param threshold Values below this threshold are ink in numeric images.
#' @param min_points Minimum traced path size.
#' @param thin Thin the mask before tracing.
#' @param simplify Optional nonnegative line-fit distance tolerance.
#' @param spec Optional feature/extraction specification, exclusive of individual
#'   extraction settings. The resulting drawing retains the specification.
#' @return A drawing with extraction provenance.
#' @export
trace_contours <- function(image,threshold=0.5,min_points=2L,thin=TRUE,simplify=NULL,spec=NULL) {
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(threshold)||!missing(min_points)||!missing(thin)||!missing(simplify))stop("use spec or individual extraction settings, not both",call.=FALSE)
  s<-spec$settings$extraction
  if(s$method!="pixel_graph_no_corner_cut")stop("unsupported extraction method",call.=FALSE)
  threshold<-s$threshold;min_points<-s$min_points;thin<-s$thin;simplify<-s$simplify
 }
 source <- NULL
 if(is.character(image)){source<-image;image<-read_mlv_image(image)}
 if(length(min_points)!=1L || !is.finite(min_points)||min_points<2 ||
    min_points!=floor(min_points))stop("min_points must be at least 2",call.=FALSE)
 if(length(threshold)!=1L||!is.finite(threshold))stop("invalid threshold",call.=FALSE)
 if(length(dim(image))==3L) {
  if(dim(image)[3]<3L)stop("RGB image needs three channels",call.=FALSE)
  image <- .2989*image[,,1]+.5870*image[,,2]+.1140*image[,,3]
 }
 mask <- if(is.logical(image)) .check_mask(image) else {
  if(!is.matrix(image)||any(!is.finite(image)))stop("invalid image",call.=FALSE)
  .check_mask(image<threshold)
 }
 if(thin)mask<-skeletonize(mask)
 paths<-cpp_trace_graph(mask,as.integer(min_points))
 contours<-lapply(paths,function(indices) {
  rc<-arrayInd(indices,dim(mask)); xy<-cbind(rc[,2],rc[,1]); n<-nrow(xy)
  cbind(xy[-n,,drop=FALSE],xy[-1,,drop=FALSE])
 })
 x<-mlv_drawing(contours,ncol(mask),nrow(mask),
    source=list(image=source,method="pixel_graph_no_corner_cut",thin=thin,threshold=threshold,
                morphology="octave_image_pratt_2.20.0",min_points=min_points,simplify=simplify))
 if(!is.null(simplify))x<-simplify_contours(x,simplify)
 if(!is.null(spec))x$feature_spec<-spec
 x
}
#' Trace RGB images with an explicit extraction provider
#' @param image Numeric RGB array or image filename.
#' @param backend Function returning a logical ink mask from the image.
#' @param ... Arguments passed to the provider.
#' @param min_points Minimum traced path size.
#' @param spec Optional specification, including explicit provider name/version.
#' @return A drawing. Model weights and preprocessing belong to the provider.
#' @export
trace_rgb <- function(image,backend,...,min_points=2L,spec=NULL) {
 if(!is.function(backend))stop("backend must be an explicitly supplied extraction function",call.=FALSE)
 if(!is.null(spec)) {
  validate_feature_spec(spec)
  if(!missing(min_points))stop("use spec or min_points, not both",call.=FALSE)
  if(is.null(spec$settings$extraction$provider))stop("RGB specification requires provider and provider_version",call.=FALSE)
 }
 if(is.character(image))image<-read_mlv_image(image)
 if(length(dim(image))!=3L||dim(image)[3]<3L||any(!is.finite(image)))
  stop("image must be a finite RGB array",call.=FALSE)
 mask<-.check_mask(backend(image,...))
 if(!identical(dim(mask),dim(image)[1:2]))stop("backend mask dimensions must match the image",call.=FALSE)
 x<-if(is.null(spec))trace_contours(mask,min_points=min_points) else trace_contours(mask,spec=spec)
 x$source$extraction_provider<-if(is.null(spec))deparse(substitute(backend)) else spec$settings$extraction$provider
 x$source$provider_version<-if(is.null(spec))NULL else spec$settings$extraction$provider_version
 x
}
#' Plot vector contours or properties
#' @param x A drawing.
#' @param property Optional segment property or mapped MAT property.
#' @param col Default contour color.
#' @param lwd Line width.
#' @param ... Arguments to the initial plot.
#' @return The drawing invisibly.
#' @export
plot.mlv_drawing <- function(x,property=NULL,col="black",lwd=1,...) {
 .check_drawing(x)
 arguments<-utils::modifyList(list(x=NA_real_,xlim=c(0,x$width),ylim=c(x$height,0),
                                   asp=1,xlab="",ylab="",type="n"),list(...))
 do.call(graphics::plot,arguments)
 st<-.measurement_table(x,x$feature_spec)
 colors<-rep(col,nrow(st))
 if(!is.null(property)) {
  if(property %in% c("length","orientation","curvature")) {
   value<-st[[property]]
   if(length(value)) {
    rr<-range(value,finite=TRUE)
    scaled<-if(diff(rr)==0)rep(0.5,length(value)) else (value-rr[1])/diff(rr)
    colors<-grDevices::hcl.colors(256,"Viridis")[1L+pmin(255L,as.integer(255*scaled))]
   }
  } else if(property=="junctions") {
   j<-.measurement_junctions(x,x$feature_spec)
   graphics::segments(st$x1,st$y1,st$x2,st$y2,col=colors,lwd=lwd)
   graphics::points(j$x,j$y,pch=21,bg="red")
   return(invisible(x))
  } else if(!is.null(x$mat_properties[[property]])) {
   .checked_mat_measurements(x,x$feature_spec)
   image<-x$mat_images[[property]];where<-which(is.finite(image));rc<-arrayInd(where,dim(image))
   graphics::segments(st$x1,st$y1,st$x2,st$y2,col="gray70",lwd=lwd)
   palette<-grDevices::hcl.colors(256,"Viridis")
   mapped_colors<-palette[pmax(1L,pmin(256L,1L+as.integer(image[where]*255)))]
   graphics::points(rc[,2],rc[,1],pch=15,cex=.4,col=mapped_colors)
   return(invisible(x))
  } else stop("unknown or unmapped plot property",call.=FALSE)
 }
 graphics::segments(st$x1,st$y1,st$x2,st$y2,col=colors,lwd=lwd)
 invisible(x)
}
#' Feature density on a pixel grid
#' @param x A drawing.
#' @param property Segment property, junctions, or a mapped MAT property.
#' @param sigma Gaussian smoothing standard deviation in pixels; zero disables.
#' @return Numeric image matrix. Kernel boundaries are zero padded. With a
#'   specification, curvature distributes integrated absolute turning over
#'   raster pixels; the historical path distributes segment curvature scores.
#' @export
feature_density <- function(x,property="length",sigma=0) {
 .check_drawing(x)
 if(length(sigma)!=1L||!is.finite(sigma)||sigma<0)stop("sigma must be nonnegative",call.=FALSE)
 out<-matrix(0,x$height,x$width);st<-.measurement_table(x,x$feature_spec)
 if(property=="junctions") {
  j<-.measurement_junctions(x,x$feature_spec)
  if(nrow(j))for(k in seq_len(nrow(j))) {
   rr<-round(j$y[k]);cc<-round(j$x[k])
   if(rr>=1&&rr<=x$height&&cc>=1&&cc<=x$width)out[rr,cc]<-out[rr,cc]+1
  }
 } else if(property %in% c("length","orientation","curvature")) {
  for(i in seq_len(nrow(st))) {
   index<-cpp_raster_indices(as.matrix(st[i,c("x1","y1","x2","y2")]),x$height,x$width,1)
   value<-st[[property]][i]
   if(property=="curvature"&&!is.null(x$feature_spec))value<-value*st$length[i]
   if(length(index))out[index]<-out[index]+value/length(index)
  }
 } else if(!is.null(x$mat_properties[[property]])) {
  .checked_mat_measurements(x,x$feature_spec)
  for(i in seq_along(x$ids)) {
   index<-.raster_indices(x,i);v<-x$mat_properties[[property]][[i]]
   ok<-is.finite(v);out[index[ok]]<-pmax(out[index[ok]],v[ok])
  }
 } else stop("unknown or unmapped property",call.=FALSE)
 if(sigma>0) {
  radius<-ceiling(3*sigma);g<-exp(-seq(-radius,radius)^2/(2*sigma^2));g<-g/sum(g)
  smooth<-function(v) {
   padded<-c(rep(0,radius),v,rep(0,radius))
   as.numeric(stats::filter(padded,g,sides=2))[(radius+1L):(radius+length(v))]
  }
  out<-t(apply(out,1L,smooth));out<-apply(out,2L,smooth)
  dim(out)<-c(x$height,x$width)
 }
 out
}
