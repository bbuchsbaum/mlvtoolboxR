// Flow-based line extraction kernels: Bradley Buchsbaum, 2026, GPL-3.0-or-later.
// Implemented from the published definitions of the structure-tensor edge
// tangent flow (Kyprianidis & Dollner 2008), the anisotropic Kuwahara filter
// with polynomial sector weights (Kyprianidis, Kang & Dollner 2009;
// Kyprianidis 2010), the flow-based difference of Gaussians (Kang, Lee & Chui
// 2007) and the XDoG threshold (Winnemoller, Kyprianidis & Olsen 2012).
//
// Images are column-major h x w (row = y, column = x), channels stacked.
// All loops are sequential and in double precision. Arithmetic is arranged
// so that, without the Kuwahara pre-filter, the eight lattice symmetries
// (quarter turns and reflections) map every intermediate value exactly, not
// only up to rounding: symmetric kernel taps are summed in mirror pairs,
// filters accumulate deviations from the centre value (so constant regions
// stay exactly constant), orientations use the doubled-angle form of the
// structure tensor, and bilinear samples sum their four terms in diagonal
// pairs, which every lattice symmetry permutes. The Kuwahara filter sums over
// its neighbourhood in an order the symmetries change, so with it the values
// agree only to rounding (about 1e-14) and mask symmetry is empirical.
//
// Line placement and stroke centre lines (cpp_band_edge) are this package's
// own steps, not part of the published methods.
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

namespace {

inline int clampi(int v,int lo,int hi){return v<lo?lo:(v>hi?hi:v);}

// Half of a sampled Gaussian, normalized over its full support [-r, r].
std::vector<double> half_gauss(double sigma,int r){
 std::vector<double> k(r+1);
 if(sigma<=0){k.assign(r+1,0.0);k[0]=1;return k;}
 double s=0;
 for(int j=0;j<=r;++j){k[j]=std::exp(-0.5*j*j/(sigma*sigma));s+=j?2*k[j]:k[j];}
 for(double &v:k)v/=s;
 return k;
}

// One separable pass with replicated borders. out = c + sum_j k_j((x_{+j}-c)+(x_{-j}-c)).
void blur_pass(const double* in,double* out,int h,int w,const std::vector<double>& k,bool along_x){
 int r=int(k.size())-1;
 for(int x=0;x<w;++x)for(int y=0;y<h;++y){
  double c=in[y+h*x],acc=0;
  for(int j=1;j<=r;++j){
   double a,b;
   if(along_x){a=in[y+h*clampi(x+j,0,w-1)];b=in[y+h*clampi(x-j,0,w-1)];}
   else{a=in[clampi(y+j,0,h-1)+h*x];b=in[clampi(y-j,0,h-1)+h*x];}
   acc+=k[j]*((a-c)+(b-c));
  }
  out[y+h*x]=c+acc;
 }
}

// Isotropic Gaussian blur; the mean of both pass orders makes it exactly
// symmetric under transposition.
std::vector<double> blur2d(const std::vector<double>& in,int h,int w,double sigma){
 if(sigma<=0)return in;
 int r=std::max(1,int(std::ceil(3*sigma)));
 std::vector<double> k=half_gauss(sigma,r),t(in.size()),xy(in.size()),yx(in.size());
 blur_pass(in.data(),t.data(),h,w,k,true);blur_pass(t.data(),xy.data(),h,w,k,false);
 blur_pass(in.data(),t.data(),h,w,k,false);blur_pass(t.data(),yx.data(),h,w,k,true);
 for(size_t i=0;i<in.size();++i)xy[i]=0.5*(xy[i]+yx[i]);
 return xy;
}

// Bilinear sample at base pixel (bx,by) plus offset (ox,oy), replicated borders.
inline double sample(const double* img,int h,int w,int bx,int by,double ox,double oy){
 double fx=std::floor(ox),fy=std::floor(oy);
 int nx0=int(fx),ny0=int(fy);
 double wx0=1-(ox-fx),wx1=1-((fx+1)-ox),wy0=1-(oy-fy),wy1=1-((fy+1)-oy);
 int x0=clampi(bx+nx0,0,w-1),x1=clampi(bx+nx0+1,0,w-1),y0=clampi(by+ny0,0,h-1),y1=clampi(by+ny0+1,0,h-1);
 double t00=(wx0*wy0)*img[y0+h*x0],t11=(wx1*wy1)*img[y1+h*x1],
        t10=(wx1*wy0)*img[y0+h*x1],t01=(wx0*wy1)*img[y1+h*x0];
 return (t00+t11)+(t10+t01);
}

// Smoothed structure tensor of a multichannel image (Sobel gradients summed
// over channels) and its orientation: unit gradient direction (gx, gy) of the
// major eigenvector, and anisotropy (l1 - l2)/(l1 + l2). The edge tangent is
// (-gy, gx). Isotropic pixels get a zero direction and zero anisotropy.
void tensor_flow(const std::vector<double>& img,int nc,int h,int w,double sigma_c,
                 std::vector<double>& gx,std::vector<double>& gy,std::vector<double>& an){
 size_t n=size_t(h)*w;
 std::vector<double> E(n,0.0),F(n,0.0),G(n,0.0);
 for(int ch=0;ch<nc;++ch){
  const double* I=img.data()+n*ch;
  for(int x=0;x<w;++x)for(int y=0;y<h;++y){
   int xm=clampi(x-1,0,w-1),xp=clampi(x+1,0,w-1),ym=clampi(y-1,0,h-1),yp=clampi(y+1,0,h-1);
   double dx=((I[ym+h*xp]+I[yp+h*xp])+2*I[y+h*xp])-((I[ym+h*xm]+I[yp+h*xm])+2*I[y+h*xm]);
   double dy=((I[yp+h*xm]+I[yp+h*xp])+2*I[yp+h*x])-((I[ym+h*xm]+I[ym+h*xp])+2*I[ym+h*x]);
   dx*=0.125;dy*=0.125;size_t i=y+size_t(h)*x;
   E[i]+=dx*dx;F[i]+=dx*dy;G[i]+=dy*dy;
  }
 }
 E=blur2d(E,h,w,sigma_c);F=blur2d(F,h,w,sigma_c);G=blur2d(G,h,w,sigma_c);
 gx.assign(n,0.0);gy.assign(n,0.0);an.assign(n,0.0);
 for(size_t i=0;i<n;++i){
  double d=E[i]-G[i],f2=2*F[i],r=std::sqrt(d*d+f2*f2),tr=E[i]+G[i];
  if(!(r>0)||!(tr>0))continue;
  double c=d/r,s=f2/r,sg=s<0?-1.0:1.0;  // the direction is axial: s = 0 takes +
  gx[i]=std::sqrt(std::max(0.0,(1+c)*0.5));gy[i]=sg*std::sqrt(std::max(0.0,(1-c)*0.5));
  an[i]=std::min(1.0,r/tr);
 }
}

// Polynomial sector weights for 8 sectors (Kyprianidis 2010) at unit-disc
// coordinates (u, v); u runs along the edge tangent. Sector k points at angle
// k * 45 degrees. Returns the weights normalized to sum to one.
inline void sector_weights(double u,double v,double zeta,double eta,double* wk){
 double s=0,z;
 double uu=zeta-eta*u*u,vv=zeta-eta*v*v;
 z=std::max(0.0, u+vv);wk[0]=z*z;
 z=std::max(0.0, v+uu);wk[2]=z*z;
 z=std::max(0.0,-u+vv);wk[4]=z*z;
 z=std::max(0.0,-v+uu);wk[6]=z*z;
 const double h=0.70710678118654752440;
 double p=h*(u-v),q=h*(u+v);
 uu=zeta-eta*p*p;vv=zeta-eta*q*q;
 z=std::max(0.0, p+vv);wk[7]=z*z;
 z=std::max(0.0, q+uu);wk[1]=z*z;
 z=std::max(0.0,-p+vv);wk[3]=z*z;
 z=std::max(0.0,-q+uu);wk[5]=z*z;
 for(int k=0;k<8;++k)s+=wk[k];
 if(s>0){double is=1/s;for(int k=0;k<8;++k)wk[k]*=is;}
}

// Entry-point guard: a positive h x w image with nc channels stored in len values.
void check_image(int h,int w,int nc,R_xlen_t len){
 if(h<=0||w<=0||nc<=0)stop("image dimensions and channel count must be positive");
 if(double(h)*double(w)*double(nc)!=double(len))stop("image buffer does not match its dimensions");
}
void check_plane(const NumericMatrix& m,int h,int w,const char* name){
 if(m.nrow()!=h||m.ncol()!=w)stop("%s must match the image dimensions",name);
}

}  // namespace

// [[Rcpp::export]]
List cpp_line_flow(NumericVector image,int h,int w,int nc,double sigma_c){
 check_image(h,w,nc,image.size());
 std::vector<double> img(image.begin(),image.end()),gx,gy,an;
 tensor_flow(img,nc,h,w,sigma_c,gx,gy,an);
 NumericMatrix tx(h,w),ty(h,w),A(h,w);
 for(size_t i=0;i<gx.size();++i){tx[i]=-gy[i];ty[i]=gx[i];A[i]=an[i];}
 return List::create(_["x"]=tx,_["y"]=ty,_["anisotropy"]=A);
}

// Anisotropic Kuwahara filter. Each pixel's neighbourhood is an ellipse
// aligned with the edge tangent, semi-axes radius*(1+A) along and
// radius/(1+A) across the flow, split into 8 overlapping polynomial sectors
// with a fall-off 1 - 0.5426 r^2 (r the unit-disc radius; it matches the
// reference implementation's Gaussian exp(-0.78125 r^2) at the centre and the
// rim, within 0.03 between); the output is the mean of the sector means
// weighted by 1/(1 + (s_k^2/s0^2)^(q/2)), s_k^2 the sector variance summed
// over channels. Isotropic pixels use a disc with axis-aligned sectors.
// [[Rcpp::export]]
NumericVector cpp_kuwahara(NumericVector image,int h,int w,int nc,double sigma_c,
                           double radius,double q,double s0){
 check_image(h,w,nc,image.size());
 size_t n=size_t(h)*w;
 std::vector<double> img(image.begin(),image.end()),gx,gy,an;
 NumericVector out(image.size());
 if(radius<=0){std::copy(img.begin(),img.end(),out.begin());return out;}
 tensor_flow(img,nc,h,w,sigma_c,gx,gy,an);
 const double pi=3.14159265358979323846,zeta=0.33;
 const double eta=(zeta+std::cos(pi/8))/std::pow(std::sin(pi/8),2);
 std::vector<double> W(8),M(8*nc),S(8),wk(8),c0(nc),dv(nc),dvm(nc),num(nc);
 double s02=s0*s0;bool q8=q==8;
 for(int x=0;x<w;++x){
  if((x&31)==0)checkUserInterrupt();
  for(int y=0;y<h;++y){
   size_t i=y+size_t(h)*x;
   double A=an[i],gxi=gx[i],gyi=gy[i];
   double tx=-gyi,ty=gxi;
   if(gxi==0&&gyi==0){tx=1;ty=0;gxi=0;gyi=1;}
   double a=radius*(1+A),b=radius/(1+A);
   // Quadratic form of the ellipse in pixel offsets, for the visited range of
   // dy at each dx; inclusion itself is decided by u^2 + v^2 <= 1 below.
   double qa=ty*ty/(a*a)+gyi*gyi/(b*b),qb=tx*ty/(a*a)+gxi*gyi/(b*b),qc=tx*tx/(a*a)+gxi*gxi/(b*b);
   int ex=int(std::ceil(std::sqrt(a*a*tx*tx+b*b*gxi*gxi))),ey=int(std::ceil(std::sqrt(a*a*ty*ty+b*b*gyi*gyi)));
   std::fill(W.begin(),W.end(),0.0);std::fill(M.begin(),M.end(),0.0);std::fill(S.begin(),S.end(),0.0);
   for(int ch=0;ch<nc;++ch)c0[ch]=img[i+n*ch];
   bool inside=x-ex>=0&&x+ex<w&&y-ey>=0&&y+ey<h;
   // Offsets d and -d are visited together: their unit-disc coordinates are
   // negatives, so the sector weights of -d are those of d rotated by four
   // sectors. Only dx > 0, or dx = 0 and dy > 0, is enumerated; the centre
   // pixel adds equally to every sector.
   {
    double g=1.0/8;
    for(int k=0;k<8;++k){W[k]+=g;}
   }
   for(int dx=0;dx<=ex;++dx){
    double disc=qb*qb*dx*dx-qa*(qc*dx*dx-1);
    if(disc<0)disc=0;
    double rt=std::sqrt(disc);
    int lo=int(std::floor((-qb*dx-rt)/qa))-1,hi=int(std::ceil((-qb*dx+rt)/qa))+1;
    if(dx==0)lo=std::max(lo,1);
    int xp=inside?x+dx:clampi(x+dx,0,w-1),xm=inside?x-dx:clampi(x-dx,0,w-1);
    for(int dy=lo;dy<=hi;++dy){
     double u=(dx*tx+dy*ty)/a,v=(dx*gxi+dy*gyi)/b,rr=u*u+v*v;
     if(rr>1)continue;
     sector_weights(u,v,zeta,eta,wk.data());
     double g=1-0.5426*rr;
     size_t jp=(inside?y+dy:clampi(y+dy,0,h-1))+size_t(h)*xp;
     size_t jm=(inside?y-dy:clampi(y-dy,0,h-1))+size_t(h)*xm;
     if(nc==1){
      double dp=img[jp]-c0[0],dm=img[jm]-c0[0],dp2=dp*dp,dm2=dm*dm;
      for(int k=0;k<8;++k){
       int kk=(k+4)&7;double wg=wk[k]*g;
       W[k]+=wg;M[k]+=wg*dp;S[k]+=wg*dp2;
       W[kk]+=wg;M[kk]+=wg*dm;S[kk]+=wg*dm2;
      }
     }else{
      double sp=0,sm=0;
      for(int ch=0;ch<nc;++ch){
       dv[ch]=img[jp+n*ch]-c0[ch];sp+=dv[ch]*dv[ch];
       dvm[ch]=img[jm+n*ch]-c0[ch];sm+=dvm[ch]*dvm[ch];
      }
      for(int k=0;k<8;++k){
       int kk=(k+4)&7;double wg=wk[k]*g;
       W[k]+=wg;S[k]+=wg*sp;W[kk]+=wg;S[kk]+=wg*sm;
       for(int ch=0;ch<nc;++ch){M[k+8*ch]+=wg*dv[ch];M[kk+8*ch]+=wg*dvm[ch];}
      }
     }
    }
   }
   double aw=0;std::fill(num.begin(),num.end(),0.0);
   for(int k=0;k<8;++k){
    if(!(W[k]>0))continue;
    double var=S[k]/W[k];
    for(int ch=0;ch<nc;++ch){double m=M[k+8*ch]/W[k];var-=m*m;}
    var=std::max(0.0,var);
    double t=var/s02,p;
    if(q8){p=t*t;p*=p;}else p=std::pow(t,0.5*q);
    double alpha=1/(1+p);
    aw+=alpha;
    for(int ch=0;ch<nc;++ch)num[ch]+=alpha*(M[k+8*ch]/W[k]);
   }
   for(int ch=0;ch<nc;++ch)out[i+n*ch]=c0[ch]+(aw>0?num[ch]/aw:0.0);
  }
 }
 return out;
}

// Flow-based difference of Gaussians. dog: 1-D DoG G_sigma - tau G_{k sigma}
// sampled across the flow (along the gradient direction, unit steps,
// bilinear); blur: the 1-D G_{k sigma} on the same samples. Both are then
// integrated with a Gaussian of scale sigma_m along edge-tangent streamlines
// (unit midpoint steps, nearest-pixel tangent, half-length ceil(2 sigma_m),
// ending early where the streamline has turned by more than acos(cos_turn)
// from its starting direction), giving response and blur_m.
// [[Rcpp::export]]
List cpp_fdog(NumericMatrix L,NumericVector smooth,int nc,double sigma_c,double sigma,
              double k,double tau,double sigma_m,double cos_turn=-2){
 int h=L.nrow(),w=L.ncol();size_t n=size_t(h)*w;
 check_image(h,w,nc,smooth.size());
 std::vector<double> img(smooth.begin(),smooth.end()),gx,gy,an;
 tensor_flow(img,nc,h,w,sigma_c,gx,gy,an);
 int J=std::max(1,int(std::ceil(3*k*sigma)));
 std::vector<double> gs=half_gauss(sigma,J),gk=half_gauss(k*sigma,J),f(J+1);
 for(int j=1;j<=J;++j)f[j]=gs[j]-tau*gk[j];
 NumericMatrix D(h,w),B(h,w),H(h,w),Bm(h,w);
 const double* I=L.begin();
 for(int x=0;x<w;++x){
  if((x&63)==0)checkUserInterrupt();
  for(int y=0;y<h;++y){
   size_t i=y+size_t(h)*x;double c=I[i],dd=0,bb=0;
   for(int j=1;j<=J;++j){
    double a=sample(I,h,w,x,y,j*gx[i],j*gy[i])-c,b=sample(I,h,w,x,y,-j*gx[i],-j*gy[i])-c;
    dd+=f[j]*(a+b);bb+=gk[j]*(a+b);
   }
   D[i]=(1-tau)*c+dd;B[i]=c+bb;
  }
 }
 if(!(sigma_m>0)){std::copy(D.begin(),D.end(),H.begin());std::copy(B.begin(),B.end(),Bm.begin());}
 else{
  int S=std::max(1,int(std::ceil(2*sigma_m)));
  std::vector<double> wm(S+1);for(int s=0;s<=S;++s)wm[s]=std::exp(-0.5*s*s/(sigma_m*sigma_m));
  const double* Dp=D.begin();const double* Bp=B.begin();
  std::vector<double> fd(S+1),bd(S+1),fb(S+1),bb(S+1);std::vector<char> hf(S+1),hb(S+1);
  for(int x=0;x<w;++x){
   if((x&63)==0)checkUserInterrupt();
   for(int y=0;y<h;++y){
    size_t i=y+size_t(h)*x;double c=Dp[i],cb=Bp[i];
    for(int dir=0;dir<2;++dir){
     std::vector<double>& vd=dir?bd:fd;std::vector<double>& vb=dir?bb:fb;std::vector<char>& has=dir?hb:hf;
     std::fill(has.begin(),has.end(),0);
     double sgn=dir?-1.0:1.0,dxs=sgn*(-gy[i]),dys=sgn*gx[i],ox=0,oy=0;
     if(dxs==0&&dys==0)continue;
     const double t0x=dxs,t0y=dys;
     // Tangent at the pixel nearest offset (ux, uy), signed to agree with
     // (dxs, dys); the current direction where the tensor is isotropic.
     auto tangent=[&](double ux,double uy,double& tx,double& ty)->bool{
      int px=x+int(std::round(ux)),py=y+int(std::round(uy));
      if(px<0||px>=w||py<0||py>=h)return false;
      size_t j=py+size_t(h)*px;tx=-gy[j];ty=gx[j];
      if(tx==0&&ty==0){tx=dxs;ty=dys;return true;}
      if(tx*dxs+ty*dys<0){tx=-tx;ty=-ty;}
      return true;
     };
     for(int s=1;s<=S;++s){
      // Midpoint (second-order) step of unit length along the flow.
      double mx,my;
      if(!tangent(ox+0.5*dxs,oy+0.5*dys,mx,my))break;
      ox+=mx;oy+=my;dxs=mx;dys=my;
      int px=x+int(std::round(ox)),py=y+int(std::round(oy));
      if(px<0||px>=w||py<0||py>=h)break;
      vd[s]=sample(Dp,h,w,x,y,ox,oy)-c;vb[s]=sample(Bp,h,w,x,y,ox,oy)-cb;has[s]=1;
      double tx,ty;tangent(ox,oy,tx,ty);dxs=tx;dys=ty;
      if(dxs*t0x+dys*t0y<cos_turn)break;
     }
    }
    double ad=0,ab=0,wt=wm[0];
    for(int s=1;s<=S;++s){
     if(hf[s]&&hb[s]){ad+=wm[s]*(fd[s]+bd[s]);ab+=wm[s]*(fb[s]+bb[s]);wt+=2*wm[s];}
     else if(hf[s]){ad+=wm[s]*fd[s];ab+=wm[s]*fb[s];wt+=wm[s];}
     else if(hb[s]){ad+=wm[s]*bd[s];ab+=wm[s]*bb[s];wt+=wm[s];}
    }
    H[i]=c+ad/wt;Bm[i]=cb+ab/wt;
   }
  }
 }
 NumericMatrix tx(h,w),ty(h,w),A(h,w);
 for(size_t i=0;i<n;++i){tx[i]=-gy[i];ty[i]=gx[i];A[i]=an[i];}
 return List::create(_["dog"]=D,_["blur"]=Bm,_["response"]=H,
                     _["flow"]=List::create(_["x"]=tx,_["y"]=ty,_["anisotropy"]=A));
}

// Lines from the ink band (this package's placement step; the band itself is
// the XDoG threshold of the flow-based DoG). The band lies on the dark side of
// an edge, about 1.2 sigma deep, and covers a thin dark stroke.
//
// Edges: for each band pixel p and each 4-neighbour q outside the band that
// is brighter in the blur, the edge lies between them. Its position is taken
// from the zero crossing of `zero` (the response of the pure DoG, tau = 1,
// whose zero crossing lies on an isolated step edge) when that changes sign
// from p to q, and from the threshold crossing of the response otherwise, by
// linear interpolation; the pixel nearer the crossing is kept (both at the
// midpoint).
//
// Strokes: from each pixel the response is sampled bilinearly at unit steps
// along +/- the gradient direction (across the flow) until it crosses the
// threshold, located by linear interpolation between samples. A band pixel
// whose band is left within max_width pixels (total width) on both sides into
// brighter blur is in a dark stroke. A non-band pixel whose non-band run is
// left within max_width on both sides into darker blur, at band pixels that
// are not in dark strokes, is in a light stroke (a light line on a dark
// ground; the gap between two dark strokes is not one). A stroke pixel is a
// centre pixel when its distance from the stroke's midpoint, measured across
// the flow, is at most max(|gx|, |gy|) / 2 + 0.05, which makes the centre line
// one or two pixels wide and 8-connected at any orientation. Centre pixels in 8-connected
// groups of at least min_size pixels are drawn, and edge pixels within
// floor(max_width / 2) + 1 pixels (8-neighbour steps) of them are dropped, so
// a stroke is drawn once, along its middle. Fragmentary stroke detections
// leave the edges in place.
// [[Rcpp::export]]
LogicalMatrix cpp_band_edge(LogicalMatrix band,NumericMatrix blur,NumericMatrix response,
                            NumericMatrix zero,double thr,NumericMatrix tx,NumericMatrix ty,
                            double max_width,double light_width,int min_size){
 int h=band.nrow(),w=band.ncol();size_t n=size_t(h)*w;
 check_image(h,w,1,band.size());
 check_plane(blur,h,w,"blur");check_plane(response,h,w,"response");check_plane(zero,h,w,"zero");
 check_plane(tx,h,w,"tx");check_plane(ty,h,w,"ty");
 const double* H=response.begin();const double* B=blur.begin();const double* Z=zero.begin();
 LogicalMatrix centre(h,w),edge(h,w);
 std::vector<char> dark(n,0);
 // Walk from pixel (x, y) along sg * (gxu, gyu) while the band state equals
 // `inside`; returns the crossing distance (0 if none within M steps or the
 // walk leaves the image) and the blur and nearest pixel just past it.
 auto walk=[&](const double* F,double level,int steps,int x,int y,double gxu,double gyu,double sg,bool inside,
               double& dist,double& exit_blur,long& exit_pixel){
  size_t i=y+size_t(h)*x;double prev=F[i];dist=0;
  for(int st=1;st<=steps;++st){
   double ox=sg*st*gxu,oy=sg*st*gyu;
   int px=x+int(std::round(ox)),py=y+int(std::round(oy));
   if(px<0||px>=w||py<0||py>=h)return;
   double v=sample(F,h,w,x,y,ox,oy);
   if((v<level)!=inside){
    double t=(level-prev)/(v-prev);if(!(t>=0&&t<=1))t=0.5;
    dist=(st-1)+t;exit_blur=sample(B,h,w,x,y,ox,oy);exit_pixel=long(py+size_t(h)*px);
    return;
   }
   prev=v;
  }
 };
 auto is_centre=[&](double d0,double d1,double gxu,double gyu){
  // An even-width stroke has two middle pixels exactly +/- 1/2 from the
  // midpoint; the 0.05 px allowance keeps both despite interpolation error.
  return std::abs(d0-d1)*0.5<=0.5*std::max(std::abs(gxu),std::abs(gyu))+0.05;
 };
 for(int pass=0;pass<2;++pass)for(int x=0;x<w;++x)for(int y=0;y<h;++y){
  size_t i=y+size_t(h)*x;bool in=band[i];
  if(pass==0&&!in)continue;
  if(pass==1&&(in||Z[i]<0))continue;
  double gxu=ty[i],gyu=-tx[i];
  if(gxu==0&&gyu==0)continue;
  double d[2],eb[2];long ep[2]={-1,-1};
  double limit=in?max_width:light_width;
  int steps=int(std::floor(limit))+1;
  for(int side=0;side<2;++side){
   if(in)walk(H,thr,steps,x,y,gxu,gyu,side?-1.0:1.0,true,d[side],eb[side],ep[side]);
   else walk(Z,0.0,steps,x,y,gxu,gyu,side?-1.0:1.0,false,d[side],eb[side],ep[side]);
  }
  if(!(d[0]>0&&d[1]>0&&d[0]+d[1]<=limit))continue;
  if(in){
   if(!(eb[0]>B[i]&&eb[1]>B[i]))continue;
   dark[i]=1;
  }else{
   if(!(eb[0]<B[i]&&eb[1]<B[i]))continue;
   if(dark[ep[0]]||dark[ep[1]])continue;
  }
  if(is_centre(d[0],d[1],gxu,gyu))centre[i]=true;
 }
 for(int x=0;x<w;++x)for(int y=0;y<h;++y){
  size_t i=y+size_t(h)*x;if(!band[i])continue;
  for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy){
   if((dx!=0)==(dy!=0))continue;int xx=x+dx,yy=y+dy;
   if(xx<0||xx>=w||yy<0||yy>=h)continue;
   size_t j=yy+size_t(h)*xx;
   if(band[j]||!(blur[j]>blur[i]))continue;
   double a,b;
   if(zero[i]<0&&zero[j]>=0){a=-zero[i];b=zero[j];}
   else{a=thr-response[i];b=response[j]-thr;}
   if(a<=b)edge[i]=true;
   if(a>=b)edge[j]=true;
  }
 }
 // Keep centre-line groups of at least min_size pixels.
 std::vector<char> seen(n,0),keep(n,0);std::vector<int> queue;
 for(size_t s=0;s<n;++s){
  if(!centre[s]||seen[s])continue;
  queue.clear();queue.push_back(int(s));seen[s]=1;
  for(size_t qi=0;qi<queue.size();++qi){
   int p=queue[qi],r=p%h,c=p/h;
   for(int dc=-1;dc<=1;++dc)for(int dr=-1;dr<=1;++dr){
    int rr=r+dr,cc=c+dc;if(rr<0||rr>=h||cc<0||cc>=w)continue;
    int j=rr+h*cc;if(centre[j]&&!seen[j]){seen[j]=1;queue.push_back(j);}
   }
  }
  if(int(queue.size())>=min_size)for(int p:queue)keep[p]=1;
 }
 int reach=int(std::floor(max_width/2))+1;
 // Square (2 reach + 1)^2 dilation of keep, as two separable sliding-window
 // passes: O(h w) whatever the reach.
 std::vector<char> across(n,0),near(n,0);
 std::vector<int> run;
 for(int y=0;y<h;++y){
  if((y&255)==0)checkUserInterrupt();
  run.assign(size_t(w)+1,0);
  for(int x=0;x<w;++x)run[x+1]=run[x]+(keep[y+size_t(h)*x]?1:0);
  for(int x=0;x<w;++x){
   int lo=std::max(0,x-reach),hi=std::min(w-1,x+reach);
   across[y+size_t(h)*x]=run[hi+1]-run[lo]>0;
  }
 }
 for(int x=0;x<w;++x){
  if((x&255)==0)checkUserInterrupt();
  run.assign(size_t(h)+1,0);
  for(int y=0;y<h;++y)run[y+1]=run[y]+(across[y+size_t(h)*x]?1:0);
  for(int y=0;y<h;++y){
   int lo=std::max(0,y-reach),hi=std::min(h-1,y+reach);
   near[y+size_t(h)*x]=run[hi+1]-run[lo]>0;
  }
 }
 LogicalMatrix out(h,w);
 for(size_t i=0;i<n;++i)out[i]=keep[i]||(edge[i]&&!near[i]);
 return out;
}

// Fill holes in an ink mask: 4-connected background components of at most
// max_size pixels that do not touch the image border. Such pinholes between
// the two pixel rows of a line would otherwise thin into small loops.
// [[Rcpp::export]]
LogicalMatrix cpp_fill_holes(LogicalMatrix mask,int max_size){
 int h=mask.nrow(),w=mask.ncol();LogicalMatrix out=clone(mask);
 std::vector<char> seen(size_t(h)*w,0);std::vector<int> queue;
 for(int s=0;s<h*w;++s){
  if((s&65535)==0)checkUserInterrupt();
  if(mask[s]||seen[s])continue;
  queue.clear();queue.push_back(s);seen[s]=1;bool border=false,big=false;
  for(size_t qi=0;qi<queue.size();++qi){
   int p=queue[qi],r=p%h,c=p/h;
   if(r==0||c==0||r==h-1||c==w-1)border=true;
   const int dr[4]={1,-1,0,0},dc[4]={0,0,1,-1};
   for(int k=0;k<4;++k){
    int rr=r+dr[k],cc=c+dc[k];if(rr<0||rr>=h||cc<0||cc>=w)continue;
    int j=rr+h*cc;if(!mask[j]&&!seen[j]){seen[j]=1;queue.push_back(j);}
   }
   if(int(queue.size())>max_size)big=true;
  }
  if(!border&&!big)for(int p:queue)out[p]=true;
 }
 return out;
}
