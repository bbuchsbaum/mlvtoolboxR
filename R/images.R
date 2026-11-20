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

#' Thin a binary image
#'
#' `"minimal"` (default) is `"hmt8"` thinning followed by sequential removal
#' (column-major scan, repeated to a fixed point) of every 8-simple pixel that
#' has at least two neighbours. No simple pixel with two or more neighbours
#' remains; a 2x2 block can remain only where every one of its pixels is
#' non-simple (removing any of them would change the 8-connected topology).
#' The result has the same 8-connected topology as the input, keeps branch end
#' points and is idempotent.
#' `"hmt8"` is sequential eight-direction hit-or-miss thinning with edge and
#' corner structuring elements, applied with background to the S, SW, E, SE, W,
#' NW, NE and N in turn. It preserves topology but is not always one pixel wide
#' (a 4x4 square leaves a 2x2 block; 4-connected staircase corners remain). That
#' order reproduces MATLAB `bwmorph(...,"skel",Inf)` exactly on the Forest_DBW
#' reference (flux skeleton and ink thinning) and is the published-compatible
#' choice; MATLAB's own algorithm is undocumented, so wider equivalence is not
#' claimed. `"pratt"` reproduces the Octave image 2.20.0
#' `bwmorph(...,"skel",Inf)` lookup tables; its skeletons are not one pixel
#' wide (strokes of even width keep two-pixel ribbons and 2x2 blocks).
#' All three are sequential or direction-ordered, so none is equivariant under
#' lattice rotations/reflections for general inputs.
#' @param mask Binary ink matrix.
#' @param min_component Minimum eight-connected component size.
#' @param method Thinning algorithm: `"minimal"` (default), `"hmt8"` or `"pratt"`.
#' @return Logical skeleton, with backend metadata described in the package docs.
#' @export
skeletonize <- function(mask,min_component=1L,method=c("minimal","hmt8","pratt")) {
 mask <- .check_mask(mask)
 method <- match.arg(method)
 if(length(min_component)!=1L || !is.finite(min_component) ||
    min_component<1 || min_component!=floor(min_component) ||
    min_component>.Machine$integer.max) stop("invalid min_component",call.=FALSE)
 result <- switch(method,pratt=cpp_skeletonize(mask),hmt8=cpp_skeletonize_hmt(mask),minimal=cpp_skeletonize_minimal(mask))
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
#' @param thin Thin the mask before tracing (the default `"minimal"` thinning, or the specification's MAT morphology).
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
 morphology<-if(is.null(spec))"minimal" else spec$settings$mat$morphology
 if(thin)mask<-skeletonize(mask,method=morphology)
 paths<-cpp_trace_graph(mask,as.integer(min_points))
 contours<-lapply(paths,function(indices) {
  rc<-arrayInd(indices,dim(mask)); xy<-cbind(rc[,2],rc[,1]); n<-nrow(xy)
  cbind(xy[-n,,drop=FALSE],xy[-1,,drop=FALSE])
 })
 x<-mlv_drawing(contours,ncol(mask),nrow(mask),
    source=list(image=source,method="pixel_graph_no_corner_cut",thin=thin,threshold=threshold,
                morphology=.morphology_definition(morphology),min_points=min_points,simplify=simplify))
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
#'
#' Draws the drawing on its canvas in image coordinates (x rightward, y
#' downward, equal aspect ratio). The canvas spans the raster extent, 0.5 to
#' `width + 0.5` and 0.5 to `height + 0.5`, because pixel centres lie at
#' integers starting at 1; the plot coordinates are those of
#' [plot.mlv_mat()], so MAT stages and drawings overlay exactly. With a `property`, segments are coloured by
#' that quantity and, by default, a colour bar naming the quantity and its unit
#' is drawn to the right of the canvas, inside the plot region, so later
#' low-level calls (`points()`, `text()`) use the same coordinates.
#'
#' Orientation is axial: a segment at 10 degrees and one at 190 degrees have
#' the same colour, and the default palette is a cyclic hue wheel over 0 to
#' 180 degrees, so 0 and 180 share a colour. Length and curvature use the
#' sequential Viridis palette over `zlim`. Mapped medial-axis scores are drawn
#' as filled mapped pixels on a 0 to 1 Viridis scale over grey contours. With
#' `property = "junctions"`, or `arms = TRUE`, junction positions are marked;
#' `arms = TRUE` also draws each arm as an arrow leaving the junction along its
#' orientation (the `arms` column of [detect_junctions()]). Junctions come from
#' the drawing's stored `junctions` (for example from
#' [compute_contour_properties()] with `rule = "mlv2023"`), its feature
#' specification, or default [detect_junctions()].
#' @param x A drawing.
#' @param property Optional segment property (`"orientation"`, `"length"`,
#'   `"curvature"`), `"junctions"`, or a mapped MAT property.
#' @param col Contour colour when no segment property is shown.
#' @param lwd Line width.
#' @param ... Arguments to the initial plot, such as `main`, `xlim` or
#'   `axes = FALSE`.
#' @param legend Draw a colour bar when colour encodes a quantity.
#' @param palette Optional palette: a function of the number of colours, or a
#'   character vector of colours. Defaults are described above.
#' @param zlim Value range mapped onto the palette for length and curvature;
#'   values outside it take the end colours. Defaults to the finite data range.
#'   Ignored for orientation (always 0 to 180 degrees) and MAT scores (0 to 1).
#' @param arms Draw junction arms.
#' @param arm_length Drawn arm length in pixels; defaults to 8% of the shorter
#'   canvas side.
#' @param junction_col Colour of junction markers and arms.
#' @param frame Outline the canvas in light grey.
#' @return The drawing invisibly.
#' @export
#' @examples
#' d <- mlv_drawing(list(rbind(c(5, 10, 35, 10)), rbind(c(20, 10, 20, 35))), 40, 40)
#' plot(d, property = "orientation")
#' plot(d, property = "junctions", arms = TRUE, axes = FALSE)
plot.mlv_drawing <- function(x,property=NULL,col="black",lwd=1,...,legend=TRUE,
                             palette=NULL,zlim=NULL,arms=FALSE,arm_length=NULL,
                             junction_col="#D55E00",frame=TRUE) {
 .check_drawing(x)
 for(flag in list(legend,arms,frame))
  if(!is.logical(flag)||length(flag)!=1L||is.na(flag))stop("legend, arms and frame must be TRUE or FALSE",call.=FALSE)
 st<-.measurement_table(x,x$feature_spec)
 colors<-rep(col,nrow(st));bar<-NULL;mat<-FALSE
 segment_property<-!is.null(property)&&property %in% c("length","orientation","curvature")
 if(!is.null(property)&&!segment_property&&property!="junctions") {
  if(is.null(x$mat_properties[[property]]))stop("unknown or unmapped plot property",call.=FALSE)
  .checked_mat_measurements(x,x$feature_spec);mat<-TRUE
 }
 if(segment_property) {
  value<-st[[property]]
  if(property=="orientation") {
   pal<-.plot_palette(palette,.orientation_colours);range<-c(0,180);value<-value%%180
   colors<-.map_colours(value,range,pal,cyclic=TRUE)
  } else {
   pal<-.plot_palette(palette,function(n)grDevices::hcl.colors(n,"Viridis"))
   range<-if(is.null(zlim)) {if(any(is.finite(value)))range(value,finite=TRUE) else c(0,1)} else zlim
   if(length(range)!=2L||any(!is.finite(range))||range[2]<range[1])stop("zlim must be an increasing finite range",call.=FALSE)
   colors<-.map_colours(value,range,pal)
  }
  bar<-list(palette=pal,range=range,label=.plot_label(x,property),cyclic=property=="orientation")
 } else if(mat) {
  pal<-.plot_palette(palette,function(n)grDevices::hcl.colors(n,"Viridis"))
  bar<-list(palette=pal,range=c(0,1),label=paste0(property,"\n(0",.range_dash(),"1)"),cyclic=FALSE)
 }
 if(!legend)bar<-NULL
 # The canvas is the raster extent: pixel centres at 1..width and 1..height,
 # so edges at 0.5 and width + 0.5, exactly as in plot.mlv_mat(), and a
 # colour bar takes the same room, so drawings and MAT stages overlay.
 side<-max(x$width,x$height)
 defaults<-list(x=NA_real_,xlim=c(.5,x$width+.5+if(is.null(bar))0 else .3*side),
                ylim=c(x$height+.5,.5),asp=1,xaxs="i",yaxs="i",xlab="",ylab="",type="n")
 arguments<-utils::modifyList(defaults,list(...))
 axes<-if(is.null(arguments$axes))TRUE else isTRUE(arguments$axes)
 arguments$axes<-FALSE;arguments$frame.plot<-FALSE
 do.call(graphics::plot,arguments)
 if(axes) {
  # Axes sit on the canvas edges, not the plot-region edges that asp = 1 moves.
  style<-list(col="grey55",col.ticks="grey55",col.axis="grey30",cex.axis=.8,tcl=-.3,mgp=c(3,.5,0))
  do.call(graphics::axis,c(list(1,at=.canvas_ticks(x$width),pos=x$height+.5),style))
  do.call(graphics::axis,c(list(2,at=.canvas_ticks(x$height),pos=.5,las=1),style))
 }
 if(frame)graphics::rect(.5,x$height+.5,x$width+.5,.5,border="grey75",lwd=.75)
 if(mat) {
  image<-x$mat_images[[property]];where<-which(is.finite(image))
  graphics::segments(st$x1,st$y1,st$x2,st$y2,col="gray70",lwd=lwd)
  # Each mapped pixel is a filled unit square on the raster extent
  # (0.5..width + 0.5), so mapped contours read as continuous strokes.
  fill<-matrix("transparent",nrow(image),ncol(image))
  fill[where]<-.map_colours(image[where],c(0,1),pal)
  graphics::rasterImage(grDevices::as.raster(fill),.5,nrow(image)+.5,ncol(image)+.5,.5,
                        interpolate=FALSE)
 } else graphics::segments(st$x1,st$y1,st$x2,st$y2,col=colors,lwd=lwd)
 if(!is.null(bar)).colour_bar(x,bar)
 if(identical(property,"junctions")||arms) {
  j<-.measurement_junctions(x,x$feature_spec)
  if(arms&&nrow(j)) {
   if(is.null(j$arms))stop("junctions lack arm directions; recompute them",call.=FALSE)
   if(is.null(arm_length))arm_length<-.08*min(x$width,x$height)
   if(length(arm_length)!=1L||!is.finite(arm_length)||arm_length<=0)stop("arm_length must be positive",call.=FALSE)
   k<-rep.int(seq_len(nrow(j)),lengths(j$arms));a<-unlist(j$arms,use.names=FALSE)*pi/180
   graphics::arrows(j$x[k],j$y[k],j$x[k]+arm_length*cos(a),j$y[k]-arm_length*sin(a),
                    length=.06,angle=25,col=junction_col,lwd=1.5*lwd)
  }
  graphics::points(j$x,j$y,pch=21,bg=junction_col,col="white",cex=1.1)
 }
 invisible(x)
}
# Cyclic hue wheel at constant chroma and luminance: hue 0 and 360 coincide,
# so axial orientations 0 and 180 degrees receive the same colour.
.orientation_colours <- function(n) grDevices::hcl(h=(15+360*(seq_len(n)-1)/n)%%360,c=70,l=55)
.plot_palette <- function(palette,default) {
 if(is.null(palette))return(default)
 if(is.function(palette))return(palette)
 if(is.character(palette)&&length(palette)>=2L)return(grDevices::colorRampPalette(palette))
 stop("palette must be a function or at least two colours",call.=FALSE)
}
# Map values onto n palette colours. Cyclic ranges wrap, so both ends of the
# range share the first colour; otherwise values are clamped to the range.
.map_colours <- function(value,range,palette,n=256L,cyclic=FALSE) {
 cols<-palette(n);out<-rep(NA_character_,length(value));ok<-is.finite(value)
 span<-diff(range)
 t<-if(span>0)(value[ok]-range[1])/span else rep(.5,sum(ok))
 index<-if(cyclic)as.integer(floor((t%%1)*n))%%n+1L else 1L+pmin(n-1L,pmax(0L,as.integer(floor(t*n))))
 out[ok]<-cols[index];out
}
.plot_label <- function(x,property) {
 unit<-if(!is.null(x$feature_spec)&&identical(x$feature_spec$settings$units,"visual_degrees"))"deg visual angle" else "px"
 switch(property,orientation="orientation\n(deg, mod 180)",length=paste0("length\n(",unit,")"),
        curvature=paste0("curvature\n(deg/",if(unit=="px")"px" else "deg visual angle",")"))
}
# An en dash for value ranges, except on the pdf() and postscript() devices,
# whose default single-byte encodings lack it and would warn and substitute.
.range_dash <- function() {
 if(names(grDevices::dev.cur())%in%c("pdf","postscript"))"-" else "\u2013"
}
.canvas_ticks <- function(n) {z<-pretty(c(1,n));unique(c(1,z[z>=1&z<=n]))}
.colour_bar <- function(x,bar,n=64L) {
 side<-max(x$width,x$height);x0<-x$width+.5+.05*side;x1<-x0+.035*side
 top<-.5+.3*x$height;bottom<-x$height+.5
 cols<-bar$palette(n)
 if(bar$cyclic)cols<-.map_colours((seq_len(n)-.5)/n*diff(bar$range)+bar$range[1],bar$range,bar$palette,cyclic=TRUE)
 # One raster, highest value in the first (top) row: adjacent rectangles
 # leave antialiasing seams on bitmap devices.
 graphics::rasterImage(grDevices::as.raster(matrix(rev(cols),ncol=1L)),x0,bottom,x1,top,
                       interpolate=FALSE,xpd=NA)
 graphics::rect(x0,top,x1,bottom,border="grey40",lwd=.5,xpd=NA)
 ticks<-if(bar$cyclic)seq(bar$range[1],bar$range[2],length.out=5L) else {
  z<-pretty(bar$range,4L);z[z>=bar$range[1]-1e-12&z<=bar$range[2]+1e-12]}
 if(diff(bar$range)>0) {
  y<-bottom-(ticks-bar$range[1])/diff(bar$range)*(bottom-top)
  graphics::segments(x1,y,x1+.012*side,y,col="grey40",lwd=.5,xpd=NA)
  graphics::text(x1+.02*side,y,format(ticks,trim=TRUE),adj=c(0,.5),cex=.7,xpd=NA)
 } else graphics::text(x1+.02*side,(top+bottom)/2,format(bar$range[1]),adj=c(0,.5),cex=.7,xpd=NA)
 graphics::text(x0,top-.04*side,bar$label,adj=c(0,0),cex=.72,xpd=NA)
 invisible(NULL)
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
