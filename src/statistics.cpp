// Exact-order histogram kernels, Bradley Buchsbaum, 2026, GPL-3.0-or-later.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {
using Accum = long double;
double finish_sum(Accum value) {
 const double limit=std::numeric_limits<double>::max();
 if(value>limit)return R_PosInf;
 if(value < -limit)return R_NegInf;
 return static_cast<double>(value);
}
void dimensions(int groups,int bins) {
 if(groups<0 || bins<1 || double(groups)*bins>std::numeric_limits<int>::max())
  stop("histogram exceeds index capacity");
}
// Euclidean remainder, using extended intermediates as R does for %%.
// Inputs here are bounded angles or bin numbers; no large-quotient case.
double remainder_positive(double x,double period) {
 // These bounded nonnegative cases are exact; retain the general arithmetic
 // for negative/subnormal values where division can underflow.
 if(x>=0 && x<period)return x;
 if(x==period)return 0;
 if(x>period && x<2*period)return x-period;
 Accum residual=static_cast<Accum>(x)-std::floor(x/period)*static_cast<Accum>(period);
 return static_cast<double>(residual-std::floor(residual/period)*period);
}
// Bin-edge snapping, identical to .snap_to_edges() in R/statistics.R: a value
// within 1e-9 * max(1, |edge|) of a finite edge is replaced by that edge, so
// values that differ from an edge only by rounding share its (upper) bin.
double edge_tolerance(double edge) {return 1e-9*std::max(1.0,std::abs(edge));}
bool near_edge(double value,double edge) {
 return std::isfinite(edge) && std::abs(value-edge)<=edge_tolerance(edge);
}
double snap_to_edges(double value,const NumericVector& edges) {
 if(std::isnan(value))return value;
 const R_xlen_t k=edges.size();
 // Number of edges <= value, as findInterval(); candidates are the edge at or
 // below the value, then the edge above it (which wins if both are near).
 const R_xlen_t i=std::upper_bound(edges.begin(),edges.end(),value)-edges.begin();
 const R_xlen_t below=std::max<R_xlen_t>(i,1)-1,above=std::min<R_xlen_t>(i+1,k)-1;
 if(near_edge(value,edges[below]))value=edges[below];
 if(near_edge(value,edges[above]))value=edges[above];
 return value;
}
int interval(double value,const NumericVector& edges) {
 // findInterval(snapped, ..., all.inside=TRUE): interior ties enter the upper bin.
 value=snap_to_edges(value,edges);
 int index=std::upper_bound(edges.begin(),edges.end(),value)-edges.begin()-1;
 return std::max(0,std::min(index,static_cast<int>(edges.size())-2));
}
void check_edges(const NumericVector& edges) {
 if(edges.size()<2 || edges.size()>std::numeric_limits<int>::max())
  stop("invalid histogram edges");
 for(R_xlen_t j=0;j<edges.size();++j)
  if(std::isnan(edges[j]) || (j && edges[j]<edges[j-1]))stop("invalid histogram edges");
}
// Two-term sum of separately rounded products, accumulated as R's sum() does.
double sum2(double a1,double b1,double a2,double b2) {
 volatile double p=a1*b1,q=a2*b2;
 Accum total=p;total+=q;return static_cast<double>(total);
}
double difference_of_products(double a1,double b1,double a2,double b2) {
 volatile double p=a1*b1,q=a2*b2;return p-q;
}
struct Piece {double value,weight;};
// Vertex-turning curvature pieces of one contour: mirrors .canonical_curve(
// reversal = TRUE) and the turning branch of .contour_profile() in
// R/curvature.R. Each vertex's absolute turn (180 for an exact reversal) is
// spread uniformly over its dual cell, from the midpoint of the segment before
// it to the midpoint of the segment after it; the end halves of an open
// contour have zero curvature. Values equal the R profile's; weights are whole
// cells, where the R profile splits each cell at its vertex.
void turning_pieces(const double* z,int n,std::vector<Piece>& out) {
 out.clear();
 const R_xlen_t st=n;const double tol=64*std::numeric_limits<double>::epsilon();
 std::vector<double> x(n+1),y(n+1);
 x[0]=z[0];y[0]=z[st];
 for(int i=0;i<n;++i){x[i+1]=z[i+2*st];y[i+1]=z[i+3*st];}
 const bool closed=x[0]==x[n] && y[0]==y[n];
 if(closed){x.pop_back();y.pop_back();}
 for(;;) {
  const int m=x.size();std::vector<char> drop(m,0);bool any=false;
  const int first=closed?0:1,last=closed?m:(m>2?m-1:1);
  for(int i=first;i<last;++i) {
   const int ip=i==0?m-1:i-1,in=i==m-1?0:i+1;
   double a1=x[i]-x[ip],a2=y[i]-y[ip],b1=x[in]-x[i],b2=y[in]-y[i];
   const double na=std::sqrt(sum2(a1,a1,a2,a2)),nb=std::sqrt(sum2(b1,b1,b2,b2));
   a1/=na;a2/=na;b1/=nb;b2/=nb;
   const double cross=difference_of_products(a1,b2,a2,b1),dot=sum2(a1,b1,a2,b2);
   if(std::abs(cross)<=tol && dot>0){drop[i]=1;any=true;}
  }
  if(!any)break;
  std::vector<double> kx,ky;
  for(int i=0;i<m;++i)if(!drop[i]){kx.push_back(x[i]);ky.push_back(y[i]);}
  x.swap(kx);y.swap(ky);
 }
 if(closed){x.push_back(x[0]);y.push_back(y[0]);}
 const int N=x.size()-1;
 std::vector<double> dx(N),dy(N),len(N),s(N+1);
 Accum run=0;s[0]=0;
 for(int k=0;k<N;++k) {
  dx[k]=x[k+1]-x[k];dy[k]=y[k+1]-y[k];
  len[k]=std::sqrt(sum2(dx[k],dx[k],dy[k],dy[k]));
  run+=len[k];s[k+1]=static_cast<double>(run);
 }
 const double L=s[N];
 if(!std::isfinite(L) || !(L>0))stop("contour length exceeds the supported numerical range");
 // Vertex v (1-based, as in R) joins segment `previous` to segment v.
 const int v0=closed?1:2;
 double lo_first=0,hi_last=0;bool any_vertex=false;
 for(int v=v0;v<=N;++v) {
  const int prev=v==1?N-1:v-2,cur=v-1;
  const double a1=dx[prev]/len[prev],a2=dy[prev]/len[prev],b1=dx[cur]/len[cur],b2=dy[cur]/len[cur];
  const double dot=sum2(a1,b1,a2,b2);
  double absolute;
  if(std::abs(difference_of_products(a1,b2,a2,b1))<=tol && dot<0)absolute=180;
  else absolute=std::abs(std::atan2(difference_of_products(a2,b1,a1,b2),dot)*180/M_PI);
  const double position=s[v-1];
  const double lo=position-len[prev]/2,hi=position+len[cur]/2;
  out.push_back({absolute/(hi-lo),hi-lo});
  if(!any_vertex)lo_first=lo;
  hi_last=hi;any_vertex=true;
 }
 if(!closed) {
  if(!any_vertex)out.push_back({0,L});
  else {out.push_back({0,lo_first});out.push_back({0,L-hi_last});}
 }
}
double segment_length(double dx,double dy,bool row_sum) {
 // Keep R's separate squared terms; do not contract multiply-add.
 volatile double xx=dx*dx,yy=dy*dy;
 if(row_sum){Accum squared=xx;squared+=yy;return std::sqrt(static_cast<double>(squared));}
 const double squared=xx+yy;
 return std::sqrt(squared);
}
}

// [[Rcpp::export]]
NumericMatrix cpp_grouped_sums(NumericVector weights,IntegerVector groups,
                               IntegerVector bins,int ng,int nb) {
 dimensions(ng,nb);
 if(weights.size()!=groups.size() || weights.size()!=bins.size())
  stop("inconsistent grouping vectors");
 NumericMatrix out(ng,nb);std::vector<Accum> sums(out.size(),0);
 for(R_xlen_t i=0;i<weights.size();++i) {
  if((i&16383)==0)checkUserInterrupt();
  // split() omits missing group labels; sum(..., na.rm=TRUE) omits missing weights.
  if(groups[i]==NA_INTEGER || bins[i]==NA_INTEGER)continue;
  if(groups[i]<1 || groups[i]>ng || bins[i]<1 || bins[i]>nb)
   stop("invalid group/bin index");
  if(!std::isnan(weights[i]))sums[(groups[i]-1)+ng*(bins[i]-1)]+=weights[i];
 }
 for(R_xlen_t i=0;i<out.size();++i)out[i]=finish_sum(sums[i]);
 return out;
}

// [[Rcpp::export]]
NumericVector cpp_contour_lengths(List contours) {
 NumericVector out(contours.size());
 for(R_xlen_t c=0;c<contours.size();++c) {
  if((c&1023)==0)checkUserInterrupt();
  SEXP cell=contours[c];
  if(!Rf_isMatrix(cell) || Rf_ncols(cell)!=4)stop("expected four contour columns");
  if(TYPEOF(cell)!=REALSXP && TYPEOF(cell)!=INTSXP)stop("expected numeric contour matrix");
  const int n=Rf_nrows(cell);
  Shield<SEXP> storage(TYPEOF(cell)==REALSXP?cell:Rf_coerceVector(cell,REALSXP));
  const double* z=REAL(storage);
  const R_xlen_t stride=n;
  Accum sum=0;
  for(int i=0;i<n;++i) {
   if((i&16383)==0)checkUserInterrupt();
   sum+=segment_length(z[i+2*stride]-z[i],z[i+3*stride]-z[i+stride],true);
  }
  out[c]=finish_sum(sum);
 }
 return out;
}

// Codes: orientation=1, contour length=2, published curvature=3,
// horizontal/vertical=4, vertex-turning curvature=5 (edges in deg/px, not
// log10; pixel units only). Only requested columns are accumulated.
// [[Rcpp::export]]
List cpp_geometry_histograms(List contours,IntegerVector codes,List edge_list,
                              bool by_drawing,double scale) {
 if(codes.size()!=edge_list.size())stop("inconsistent property definitions");
 if(!(scale>0))stop("invalid measurement scale");
 if(contours.size()>std::numeric_limits<int>::max())stop("too many contours");
 const int nc=contours.size(),ng=by_drawing?1:nc,np=codes.size();
 std::vector<NumericVector> edges;std::vector<int> widths;
 std::vector<std::vector<Accum>> sums;
 bool angles=false,segments=false;
 for(int p=0;p<np;++p) {
  if(codes[p]<1 || codes[p]>5)stop("invalid geometry property");
  if(codes[p]==5 && scale!=1)stop("turning histograms are computed in pixels");
  NumericVector edge=edge_list[p];check_edges(edge);int k=edge.size()-1;
  dimensions(ng,k);edges.push_back(edge);widths.push_back(k);
  sums.emplace_back(static_cast<size_t>(ng)*k,0);
  const bool per_segment=codes[p]==1 || codes[p]==3 || codes[p]==4;
  angles=angles || per_segment;segments=segments || per_segment;
 }
 NumericVector lengths(nc);std::vector<Piece> pieces;
 for(int c=0;c<nc;++c) {
  if((c&255)==0)checkUserInterrupt();
  SEXP cell=contours[c];
  if(!Rf_isMatrix(cell) || Rf_ncols(cell)!=4)stop("expected four contour columns");
  if(TYPEOF(cell)!=REALSXP && TYPEOF(cell)!=INTSXP)stop("expected numeric contour matrix");
  const int n=Rf_nrows(cell),group=by_drawing?0:c;
  Shield<SEXP> storage(TYPEOF(cell)==REALSXP?cell:Rf_coerceVector(cell,REALSXP));
  const double* z=REAL(storage);
  const R_xlen_t stride=n;
  std::vector<double> lens(segments?n:0),ori(angles?n:0);
  Accum total=0;
  for(int i=0;i<n;++i) {
   if((i&16383)==0)checkUserInterrupt();
   const double dx=z[i+2*stride]-z[i],dy=z[i+3*stride]-z[i+stride];
   const double raw_length=segment_length(dx,dy,false);
   if(!std::isfinite(raw_length) || raw_length<=0)
    stop("contour geometry must have finite positive segment lengths");
   if(std::numeric_limits<Accum>::digits==std::numeric_limits<double>::digits)
    total+=raw_length;
   else total+=segment_length(dx,dy,true);
   if(segments)lens[i]=raw_length;
   if(angles)ori[i]=remainder_positive(std::atan2(-dy,dx)*180/M_PI,360);
  }
  lengths[c]=finish_sum(total)*scale;
  for(int p=0;p<np;++p)if(codes[p]==2) {
   const int bin=interval(std::log10(lengths[c]+1),edges[p]);
   if(!std::isnan(lengths[c]))sums[p][group+ng*bin]+=lengths[c];
  }
  bool have_pieces=false;
  for(int p=0;p<np;++p)if(codes[p]==5) {
   if(!have_pieces){turning_pieces(z,n,pieces);have_pieces=true;}
   for(const Piece& piece:pieces)sums[p][group+ng*interval(piece.value,edges[p])]+=piece.weight;
  }
  for(int i=0;segments && i<n;++i) {
   if((i&16383)==0)checkUserInterrupt();
   const double weight=lens[i]*scale;
   const int neighbor=i+1<n?i+1:(n==1?i:i-1);
   const double delta=std::abs(ori[i]-ori[neighbor]);
   const double curvature=(std::min(delta,360-delta)/lens[i])/scale;
   for(int p=0;p<np;++p) {
    int bin;
    switch(codes[p]) {
     case 1: {
      // Edges at (m - 1/2) * 180 / bins; same snapping as .orientation_index().
      const double axial=remainder_positive(ori[i],180);
      const double scaled=axial*widths[p]/180+.5;
      const double nearest=std::floor(scaled+.5),edge=(nearest-.5)*180/widths[p];
      const double rounded=near_edge(axial,edge)?nearest:std::floor(scaled);
      bin=static_cast<int>(remainder_positive(rounded,widths[p]));break;
     }
     case 3: bin=interval(std::log10(curvature+1),edges[p]);break;
     case 4: {
      const double angle=ori[i]*M_PI/180;
      const double value=std::max(-1.0,std::min(1.0,std::abs(std::cos(angle))-std::abs(std::sin(angle))));
      bin=interval(value,edges[p]);break;
     }
     default:continue;
    }
    if(!std::isnan(weight))sums[p][group+ng*bin]+=weight;
   }
  }
 }
 List counts(np);
 for(int p=0;p<np;++p) {
  NumericMatrix out(ng,widths[p]);
  for(R_xlen_t i=0;i<out.size();++i)out[i]=finish_sum(sums[p][i]);
  counts[p]=out;
 }
 return List::create(_["lengths"]=lengths,_["counts"]=counts);
}
