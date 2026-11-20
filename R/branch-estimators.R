.score_branch_arclength <- function(branch,property,window,power,reference_radius,
                                    fit_tolerance,units,pixels_per_degree,undefined="zero") {
 X<-branch$x;Y<-branch$y;R<-branch$radius;n<-length(X)
 if(!is.numeric(X)||!is.numeric(Y)||!is.numeric(R)||length(Y)!=n||length(R)!=n||
    any(!is.finite(c(X,Y,R)))||any(R<0))stop("branch coordinates/radii must be finite, nonnegative and equally sized",call.=FALSE)
 window<-.spec_number(window,"window",strict=TRUE)
 power<-.spec_number(power,"power",strict=TRUE)
 reference_radius<-.spec_number(reference_radius,"reference_radius",strict=TRUE)
 fit_tolerance<-.spec_number(fit_tolerance,"fit_tolerance",strict=TRUE)
 factor<-.curvature_units(units,pixels_per_degree)
 if(!n)return(numeric())
 if(n==1L)return(NA_real_)
 X<-X*factor;Y<-Y*factor;R<-R*factor
 ds<-sqrt(diff(X)^2+diff(Y)^2)
 if(any(!is.finite(ds)))stop("branch exceeds the supported numerical range",call.=FALSE)
 if(any(ds==0 & diff(R)!=0))stop("duplicate coordinates have inconsistent radii",call.=FALSE)
 original_s<-c(0,cumsum(ds));keep<-c(TRUE,ds>0)
 X<-X[keep];Y<-Y[keep];R<-R[keep];s<-original_s[keep]
 if(length(s)<2L)return(rep(NA_real_,n))
 L<-tail(s,1)
 integral_square<-function(v,w)sum(w*(head(v,-1)^2+head(v,-1)*tail(v,-1)+tail(v,-1)^2)/3)
 out<-vapply(original_s,function(at) {
  a<-max(0,at-window);b<-min(L,at+window);width<-b-a
  if(width<=0)return(NA_real_)
  knots<-c(a,s[s>a&s<b],b);w<-diff(knots)
  px<-stats::approx(s,X,xout=knots)$y;py<-stats::approx(s,Y,xout=knots)$y
  radius<-stats::approx(s,R,xout=knots)$y
  mean_radius<-sum(w*(head(radius,-1)+tail(radius,-1))/2)/width
  if(property=="separation")return(mean_radius/(mean_radius+reference_radius))
  if(property=="parallelism")return(width/sum(sqrt(w^2+diff(radius)^2)))
  if(property=="mirror") {
   dx<-tail(px,1)-px[1];dy<-tail(py,1)-py[1];chord<-sqrt(dx^2+dy^2)
   # A window whose ends coincide (a closed stretch) has no chord direction:
   # undefined, NA under branch_rule = "invariant"; 0 as before otherwise.
   if(chord==0)return(if(undefined=="na")NA_real_ else 0)
   perpendicular<-((px-px[1])*dy-(py-py[1])*dx)/chord
   error2<-integral_square(perpendicular,w)/width
   return((chord/width)/(1+error2/fit_tolerance^2))
  }
  t<-knots-(a+b)/2;centered<-radius-mean_radius
  # Exact integral of a product of affine functions, in centered coordinates.
  tr<-sum(w*(2*head(t,-1)*head(centered,-1)+head(t,-1)*tail(centered,-1)+
             tail(t,-1)*head(centered,-1)+2*tail(t,-1)*tail(centered,-1))/6)
  slope<-tr/(width^3/12)
  residual<-centered-slope*t
  error2<-integral_square(residual,w)/width
  1/sqrt(1+error2/width^2)
 },double(1))
 pmax(0,pmin(1,out))^power
}
