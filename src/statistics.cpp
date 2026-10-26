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
int interval(double value,const NumericVector& edges) {
 // findInterval(..., all.inside=TRUE): interior ties enter the upper bin.
 int index=std::upper_bound(edges.begin(),edges.end(),value)-edges.begin()-1;
 return std::max(0,std::min(index,static_cast<int>(edges.size())-2));
}
void check_edges(const NumericVector& edges) {
 if(edges.size()<2 || edges.size()>std::numeric_limits<int>::max())
  stop("invalid histogram edges");
 for(R_xlen_t j=0;j<edges.size();++j)
  if(std::isnan(edges[j]) || (j && edges[j]<edges[j-1]))stop("invalid histogram edges");
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
// horizontal/vertical=4. Only requested columns are accumulated.
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
  if(codes[p]<1 || codes[p]>4)stop("invalid geometry property");
  NumericVector edge=edge_list[p];check_edges(edge);int k=edge.size()-1;
  dimensions(ng,k);edges.push_back(edge);widths.push_back(k);
  sums.emplace_back(static_cast<size_t>(ng)*k,0);
  angles=angles || codes[p]!=2;segments=segments || codes[p]!=2;
 }
 NumericVector lengths(nc);
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
      const double axial=remainder_positive(ori[i],180);
      const double rounded=std::floor(axial*widths[p]/180+.5);
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
