# Flow-based line extraction (FDoG with an XDoG threshold).

# Masks, drawings and feature specifications record the extraction definition.
.line_definitions <- c(base="fdog-xdog-1",repair="fdog-xdog-repair-1")
.line_providers <- c("xdog_lines","extract_lines")
.line_defaults <- list(sigma=2,k=1.6,tau=0.98,eps=0,phi=Inf,sigma_c=2,sigma_m=3,
                       smooth=5,colour=FALSE,min_component=20L)

.line_definition <- function(p) .line_definitions[[if(is.null(p$junction_gap))"base" else "repair"]]

.line_number <- function(x,name,lower,upper=Inf,strict=FALSE,lower_ok_inf=FALSE,zero_ok=FALSE) {
 if(zero_ok&&is.numeric(x)&&length(x)==1L&&!is.na(x)&&x==0)return(0)
 if(!is.numeric(x)||length(x)!=1L||is.na(x)||(!is.finite(x)&&!(upper==Inf&&x==Inf&&lower_ok_inf))||
    (if(strict)x<=lower else x<lower)||x>upper)
  stop(name," must be a single number",if(zero_ok)" equal to 0 or" else ""," ",
       if(strict)"greater than " else "at least ",lower,
       if(is.finite(upper))paste0(" and at most ",upper) else "",call.=FALSE)
 as.double(x)
}

# Canonical, validated parameter list; unknown names fail.
.line_parameters <- function(params=list()) {
 if(!is.list(params)||length(params)&&(is.null(names(params))||any(!nzchar(names(params))))||
    anyDuplicated(names(params)))stop("line extraction parameters must be uniquely named",call.=FALSE)
 unknown<-setdiff(names(params),c(names(.line_defaults),"junction_gap"))
 if(length(unknown))stop("unknown line extraction parameter(s): ",paste(unknown,collapse=", "),call.=FALSE)
 gap<-if("junction_gap"%in%names(params)).line_number(params$junction_gap,"junction_gap",0,50) else 0
 params$junction_gap<-NULL
 p<-.line_defaults;p[names(params)]<-params
 # Lower bounds keep every kernel at least a few samples wide (no 0/0).
 p$sigma<-.line_number(p$sigma,"sigma",0.1,50)
 p$k<-.line_number(p$k,"k",1,10,strict=TRUE)
 p$tau<-.line_number(p$tau,"tau",0,1)
 p$eps<-.line_number(p$eps,"eps",-1,1)
 p$phi<-.line_number(p$phi,"phi",0,Inf,strict=TRUE,lower_ok_inf=TRUE)
 p$sigma_c<-.line_number(p$sigma_c,"sigma_c",0.1,50)
 p$sigma_m<-.line_number(p$sigma_m,"sigma_m",0.1,50,zero_ok=TRUE)
 p$smooth<-.line_number(p$smooth,"smooth",0.5,50,zero_ok=TRUE)
 if(!is.logical(p$colour)||length(p$colour)!=1L||is.na(p$colour))stop("colour must be TRUE or FALSE",call.=FALSE)
 mc<-p$min_component
 if(!is.numeric(mc)||length(mc)!=1L||!is.finite(mc)||mc<1||mc!=floor(mc)||mc>.Machine$integer.max)
  stop("min_component must be a positive whole number",call.=FALSE)
 p$min_component<-as.integer(mc)
 # Canonical parameters of the original definition omit the inactive option.
 if(gap>0)p$junction_gap<-gap
 p
}
.line_method <- function(method) {
 if(!identical(method,"fdog"))stop("method must be \"fdog\", the only line extraction method",call.=FALSE)
 method
}
# Parameters passed through `...`: every one named exactly (no positional or
# partial matching), `method` checked and dropped.
.line_dots <- function(args) {
 if(length(args)&&(is.null(names(args))||any(!nzchar(names(args)))))
  stop("line extraction parameters must be named, e.g. sigma = 2",call.=FALSE)
 if("method"%in%names(args)){.line_method(args$method);args$method<-NULL}
 args
}

# sRGB in [0, 1] to CIE L*a*b* (D65), each divided by 100.
.srgb_lab <- function(r,g,b,colour) {
 lin<-function(v){o<-v/12.92;i<-v>0.04045;o[i]<-((v[i]+0.055)/1.055)^2.4;o}
 f<-function(t){o<-t/(3*(6/29)^2)+4/29;i<-t>(6/29)^3;o[i]<-t[i]^(1/3);o}
 r<-lin(r);g<-lin(g);b<-lin(b)
 # sRGB to XYZ, each row divided by its sum so that white maps exactly to the
 # D65 white point (L* = 100, a* = b* = 0).
 Y<-(0.2126729*r+0.7151522*g+0.0721750*b)/1.0000001
 fy<-f(Y);L<-(116*fy-16)/100
 if(!colour)return(list(L=L))
 X<-(0.4124564*r+0.3575761*g+0.1804375*b)/0.9504700
 Z<-(0.0193339*r+0.1191920*g+0.9503041*b)/1.0888300
 list(L=L,a=5*(f(X)-fy),b=2*(fy-f(Z)))
}

.line_image <- function(image) {
 if(is.character(image)&&is.null(dim(image))){
  if(length(image)!=1L||is.na(image))stop("image file name must be a single string",call.=FALSE)
  image<-read_mlv_image(image)
 }
 d<-dim(image)
 if(!is.numeric(image)||!(length(d)==2L||length(d)==3L&&d[3]%in%c(1L,3L,4L))||any(d[1:2]<1L))
  stop("image must be a numeric greyscale matrix or an RGB(A) array",call.=FALSE)
 if(prod(as.double(d[1:2]))>.Machine$integer.max)stop("image exceeds the supported integer index capacity",call.=FALSE)
 if(any(!is.finite(image))||any(image<0|image>1))stop("image values must be finite and in [0, 1]",call.=FALSE)
 if(length(d)==2L)return(list(r=unname(image),g=unname(image),b=unname(image)))
 channel<-function(k)matrix(image[,,k],d[1],d[2])
 if(d[3]==1L){v<-channel(1);return(list(r=v,g=v,b=v))}
 a<-if(d[3]==4L)channel(4) else 1
 # Transparent pixels are composited over white paper.
 list(r=channel(1)*a+1-a,g=channel(2)*a+1-a,b=channel(3)*a+1-a)
}

.line_pipeline <- function(image,p,stages=FALSE) {
 rgb<-.line_image(image);h<-nrow(rgb$r);w<-ncol(rgb$r)
 lab<-.srgb_lab(rgb$r,rgb$g,rgb$b,p$colour)
 .line_core(lab,h,w,p,stages)
}

# Pipeline from lightness (and optional a*, b*) channels, each h x w.
.line_core <- function(lab,h,w,p,stages=FALSE) {
 nc<-length(lab);stack<-as.double(unlist(lapply(lab,as.vector),use.names=FALSE))
 smoothed<-if(p$smooth>0)cpp_kuwahara(stack,h,w,nc,p$sigma_c,p$smooth,8,0.01) else stack
 L<-matrix(lab$L,h,w);Ls<-matrix(smoothed[seq_len(h*w)],h,w)
 fd<-cpp_fdog(Ls,smoothed,nc,p$sigma_c,p$sigma,p$k,p$tau,p$sigma_m,cos(pi/9))
 H<-fd$response
 # XDoG soft threshold; the band is where it falls below one half.
 thr<-if(is.finite(p$phi))p$eps-atanh(0.5)/p$phi else p$eps
 band<-H<thr
 mask<-cpp_band_edge(band,fd$blur,H,H-(1-p$tau)*fd$blur,thr,fd$flow$x,fd$flow$y,3*p$sigma,3*p$sigma,p$min_component)
 mask<-cpp_fill_holes(mask,2L)
 if(p$min_component>1L)mask<-cpp_remove_small(mask,p$min_component)
 if(!is.null(p$junction_gap))mask<-.repair_line_junctions(mask,L,p$junction_gap)
 attr(mask,"extraction")<-list(method=.line_definition(p),parameters=p)
 if(!stages)return(mask)
 response<-if(is.finite(p$phi))ifelse(H>=p$eps,1,1+tanh(p$phi*(H-p$eps))) else 1*(H>=p$eps)
 list(lightness=L,smoothed=Ls,flow=fd$flow,dog=fd$dog,response=H,xdog=response,band=band,
      mask=mask,extraction=attr(mask,"extraction"))
}

#' Extract a line drawing from an image
#'
#' Turns a photograph or an image of a painting into a binary ink mask whose
#' lines follow the image's edges, ready for [trace_contours()] or, through
#' [xdog_lines()], for [trace_rgb()]. The only method, `"fdog"`, is a
#' flow-based difference of Gaussians with an XDoG threshold, in four stages:
#'
#' 1. **Lightness and pre-filter.** The image is converted from sRGB to CIE
#'    `L*` (D65 white), divided by 100. An anisotropic Kuwahara filter
#'    (Kyprianidis, Kang and Doellner 2009) removes brush texture while keeping
#'    form boundaries: each pixel's neighbourhood is an ellipse along the local
#'    edge, of radius `smooth` (stretched to `smooth * (1 + A)` along the edge
#'    and shrunk to `smooth / (1 + A)` across it, `A` the anisotropy), split
#'    into 8 overlapping polynomial sectors (Kyprianidis 2010) with a fall-off
#'    `1 - 0.5426 r^2` (`r` the radius in the unit disc, close to the
#'    reference Gaussian). The output averages the sector means with weights
#'    `1 / (1 + (v / 0.01^2)^4)`, `v` a sector's lightness variance, so the
#'    most uniform sectors dominate. `smooth = 0` skips the filter. Like any
#'    Kuwahara filter it rounds sharp corners slightly (at the default radius
#'    the traced corner of a square lies up to about 1.2 pixels inside the true
#'    corner, against 0.85 pixels without the filter), and it removes dark
#'    strokes narrower than about 2 pixels unless they lie along a pixel row or
#'    column; light strokes survive from 1 pixel. For thin line art, use
#'    `smooth = 0`.
#' 2. **Edge tangent flow.** The structure tensor of the smoothed image
#'    (Sobel gradients, Gaussian scale `sigma_c`) gives the local edge
#'    tangent as its minor eigenvector (Kyprianidis and Doellner 2008).
#' 3. **Flow-based DoG.** At every pixel a 1-D difference of Gaussians
#'    `G_sigma - tau * G_(k sigma)` is taken across the flow (along the
#'    gradient, unit steps, bilinear samples), then integrated along the
#'    tangent streamline with a Gaussian of scale `sigma_m`, which joins edge
#'    fragments into coherent strokes (Kang, Lee and Chui 2007). Streamlines
#'    take unit midpoint steps for up to `ceil(2 * sigma_m)` steps each way and
#'    stop early once they have turned more than 20 degrees from their starting
#'    direction, so corners are not rounded off. `sigma_m = 0` skips the
#'    integration.
#' 4. **Threshold, then this package's placement and cleanup.** The XDoG soft threshold
#'    (Winnemoeller, Kyprianidis and Olsen 2012), `T = 1` where the response
#'    `D >= eps` and `1 + tanh(phi * (D - eps))` elsewhere, marks an ink band
#'    where `T < 1/2`; with `phi = Inf` (the default) that is `D < eps`, and a
#'    finite `phi` moves the threshold to `eps - atanh(1/2) / phi`. The band
#'    lies on the dark side of an edge, about `1.2 * sigma` deep, so its
#'    skeleton would sit beside the edge. The remaining steps are this
#'    package's own, not part of the published methods. Lines are placed
#'    between the band and its brighter side: for each band pixel and each
#'    4-neighbour outside the band that is brighter (after the across-flow
#'    `G_(k sigma)` blur, integrated along the flow), the mask keeps whichever
#'    of the two lies nearer the zero crossing of the pure difference of
#'    Gaussians (`tau = 1`, whose zero crossing lies on an isolated step edge),
#'    or nearer the threshold crossing where that does not change sign between
#'    them. Thin strokes are drawn once, along their middle: sampling the
#'    response bilinearly across the flow, a band at most `3 * sigma` pixels
#'    wide with brighter pixels on both sides is a dark stroke, and a gap of
#'    the pure difference of Gaussians at most `3 * sigma` wide with darker band
#'    on both sides (not itself a dark stroke) is a light stroke. Each is drawn
#'    as its centre line (if at least `min_component` pixels long), and edge
#'    lines beside it are dropped. Parallel dark strokes 2 to 4 pixels apart
#'    (hatching) merge into one band and are drawn as that band's two outer
#'    outlines, about 2 pixels outside the outer strokes, with no line between
#'    them; from about 6 pixels apart each gets its own centre line. Enclosed holes of one or two pixels are
#'    filled, and eight-connected components smaller than `min_component`
#'    pixels are removed.
#'
#' The decision is a sign test of a filter of lightness, so with `eps = 0`
#' and `phi = Inf` it compares edge contrast with the local lightness: with
#' `tau < 1`, an edge is drawn when its lightness step exceeds roughly
#' `(1 - tau) / 0.1` of the darker side's lightness (20 percent at the default
#' `tau`). Multiplying the lightness by a positive constant therefore leaves
#' the mask unchanged when the pre-filter is off (`smooth = 0`), but the method
#' is not invariant to scaling image intensities: the sRGB to `L*` conversion
#' is non-linear, and the Kuwahara weights use an absolute variance scale.
#' Raising `eps` widens the band; once it covers the whole image there is no
#' band boundary left and the mask is empty (XDoG's stylised output would be
#' all ink).
#'
#' FDoG can leave a gap at an occlusion T junction in the thresholded band.
#' A positive `junction_gap` enables a repair after mask cleanup. It fits
#' straight segments to the thinned contours and considers extending a free
#' endpoint to a nearby contour that crosses its projected path. It adds a
#' bridge when contrast across that path in the original CIE lightness
#' image has the same sign and at least 40 percent of the magnitude of the
#' median contrast along the visible stem. The stem must have at least 3.5 L*
#' units of median contrast. The new pixels must touch two distinct patches
#' of existing ink, one at each endpoint. The thresholded band stays as
#' computed by FDoG.
#' Validate junction measurements against annotations on representative
#' images. For line art whose junctions follow stroke centre lines, pass a
#' centre-line mask to \code{\link{trace_contours}}.
#'
#' The extraction is deterministic. With `smooth = 0` and `junction_gap = 0`,
#' the arithmetic is
#' arranged so that quarter turns and reflections of the image transform
#' every stage, and so the mask, exactly. The Kuwahara filter accumulates its
#' sums in an order that these transforms change, so with the pre-filter on the
#' intermediate images agree only to rounding (about 1e-14); the masks have
#' agreed exactly in every test, but a pixel within rounding of the threshold
#' could differ.
#'
#' With `colour = TRUE` the pre-filter and both structure tensors use all
#' three CIE `L*a*b*` channels (each divided by 100), so colour boundaries
#' steer the flow and the smoothing; the difference of Gaussians is still
#' taken on lightness, so a boundary between two colours of equal lightness is
#' not drawn.
#'
#' The defaults were chosen by a visual review of four public-domain
#' paintings 640 pixels wide (see the \href{https://bbuchsbaum.github.io/mlvtoolboxR/articles/line-extraction.html}{Line extraction} article).
#' `sigma`, `sigma_c`, `sigma_m` and `smooth` are in pixels, so images at other
#' resolutions want them scaled with the image.
#'
#' @param image Greyscale matrix or RGB(A) array with values in \[0, 1\]
#'   (sRGB-encoded, as returned by [read_mlv_image()]), or a PNG/JPEG file name.
#'   Transparent pixels are composited over white.
#' @param method Extraction method; only `"fdog"` (matched exactly).
#' @param sigma Scale of the narrow Gaussian across the flow, in pixels (0.1 to
#'   50); sets the scale of the edges that are drawn.
#' @param k Ratio of the wide to the narrow Gaussian scale (greater than 1).
#' @param tau Weight of the wide Gaussian, in \[0, 1\]; lower values draw only
#'   higher-contrast edges.
#' @param eps XDoG threshold on the response, in lightness units (\[-1, 1\]).
#' @param phi XDoG threshold sharpness (positive; `Inf` is a hard threshold).
#' @param sigma_c Gaussian scale of the structure tensor, in pixels (0.1 to 50).
#' @param sigma_m Scale of the integration along the flow, in pixels: 0 skips
#'   it, otherwise 0.1 to 50.
#' @param smooth Radius of the anisotropic Kuwahara pre-filter in pixels: 0
#'   skips it, otherwise 0.5 to 50.
#' @param colour Use CIE `L*a*b*` colour for the pre-filter and the flow.
#' @param min_component Smallest eight-connected ink component kept, in pixels.
#'   Narrow light gaps between dark marks are drawn as short centre lines, so
#'   lowering it below the default adds many short marks in textured areas.
#' @param junction_gap Optional maximum bridge length in pixels (0 to 50).
#'   Zero disables repair. The repair considers gaps longer than 1.5 pixels.
#' @return A logical ink mask (`TRUE` is ink) with the dimensions of the image
#'   and an `"extraction"` attribute giving the method and validated
#'   parameters. The original definition, `"fdog-xdog-1"`, omits
#'   `junction_gap`. With repair enabled, `"fdog-xdog-repair-1"` records it.
#' @references
#' Kang, H., Lee, S. and Chui, C. K. (2007). Coherent line drawing.
#' *Proceedings of the 5th International Symposium on Non-Photorealistic
#' Animation and Rendering (NPAR)*, 43-50. \doi{10.1145/1274871.1274878}
#'
#' Kyprianidis, J. E. and Doellner, J. (2008). Image abstraction by structure
#' adaptive filtering. *Proceedings of EG UK Theory and Practice of Computer
#' Graphics*, 51-58.
#'
#' Kyprianidis, J. E., Kang, H. and Doellner, J. (2009). Image and video
#' abstraction by anisotropic Kuwahara filtering. *Computer Graphics Forum*,
#' 28(7), 1955-1963. \doi{10.1111/j.1467-8659.2009.01574.x}
#'
#' Kyprianidis, J. E. (2010). Anisotropic Kuwahara filtering with polynomial
#' weighting functions. *Proceedings of EG UK Theory and Practice of Computer
#' Graphics*, 25-30.
#'
#' Winnemoeller, H., Kyprianidis, J. E. and Olsen, S. C. (2012). XDoG: an
#' eXtended difference-of-Gaussians compendium including advanced image
#' stylization. *Computers & Graphics*, 36(6), 740-753.
#' \doi{10.1016/j.cag.2012.03.004}
#' @seealso [line_extraction_stages()] for the intermediate images;
#'   [trace_contours()] and [trace_rgb()] to vectorise the mask.
#' @examples
#' # A dark disc on light paper.
#' xy <- expand.grid(y = 1:48, x = 1:48)
#' img <- matrix(ifelse((xy$x - 24.5)^2 + (xy$y - 24.5)^2 < 15^2, 0.2, 0.9), 48, 48)
#' ink <- extract_lines(img)
#' sum(ink)
#' attr(ink, "extraction")$method
#' d <- trace_contours(ink, simplify = 1)
#' length(d$contours)
#' @export
extract_lines <- function(image,method="fdog",sigma=2,k=1.6,tau=0.98,eps=0,phi=Inf,
                          sigma_c=2,sigma_m=3,smooth=5,colour=FALSE,min_component=20L,
                          junction_gap=0) {
 .line_method(method)
 p<-.line_parameters(list(sigma=sigma,k=k,tau=tau,eps=eps,phi=phi,sigma_c=sigma_c,sigma_m=sigma_m,
                          smooth=smooth,colour=colour,min_component=min_component,junction_gap=junction_gap))
 .line_pipeline(image,p)
}

#' Line extraction as a trace_rgb() backend
#'
#' `xdog_lines()` is [extract_lines()] with `method = "fdog"`, shaped as a
#' [trace_rgb()] backend: `trace_rgb(image, backend = xdog_lines, ...)` records
#' the provider `"xdog_lines"`, its extraction definition and the validated
#' parameters in the drawing's `source`.
#' @param image Greyscale matrix, RGB(A) array or image file name; see [extract_lines()].
#' @param ... Parameters of [extract_lines()], each named in full; `method`,
#'   if given, must be `"fdog"`.
#' @return The logical ink mask of [extract_lines()].
#' @examples
#' rgb <- array(0.9, c(40, 40, 3)); rgb[10:30, 12:28, ] <- 0.25
#' d <- trace_rgb(rgb, backend = xdog_lines, min_component = 10)
#' d$source$provider_version
#' @export
xdog_lines <- function(image,...) {
 p<-.line_parameters(.line_dots(list(...)))
 .line_pipeline(image,p)
}

#' Intermediate images of the line extractor
#'
#' Runs [extract_lines()] and returns every stage, for figures and teaching.
#' @inheritParams extract_lines
#' @param ... Parameters of [extract_lines()], each named in full; `method`,
#'   if given, must be `"fdog"`.
#' @return A list of `h x w` matrices unless noted: `lightness`
#'   (CIE `L*` divided by 100), `smoothed` (after the Kuwahara
#'   pre-filter), `flow` (a list of the unit edge-tangent components `x` and
#'   `y`, zero where the tensor is isotropic, and the `anisotropy` in \[0, 1\]),
#'   `dog` (the across-flow difference of Gaussians), `response` (after
#'   integration along the flow), `xdog` (the XDoG soft threshold in \[0, 1\]),
#'   `band` (the logical ink band), `mask` (the result of [extract_lines()]) and
#'   `extraction` (its definition and parameters).
#' @examples
#' img <- matrix(0.9, 32, 32); img[8:24, 8:24] <- 0.2
#' s <- line_extraction_stages(img, min_component = 5)
#' names(s)
#' identical(s$mask, extract_lines(img, min_component = 5))
#' @export
line_extraction_stages <- function(image,...) {
 p<-.line_parameters(.line_dots(list(...)))
 .line_pipeline(image,p,stages=TRUE)
}
