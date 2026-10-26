// Isolated performance experiment; not linked into mlvtoolbox.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix grouped_sums(NumericVector weights, IntegerVector groups,
                          IntegerVector bins, int ng, int nb) {
 if(ng<0 || nb<1 || double(ng)*nb>std::numeric_limits<int>::max())
   stop("invalid output dimensions");
 if(weights.size()!=groups.size() || weights.size()!=bins.size())
   stop("inconsistent grouping vectors");
 NumericMatrix out(ng,nb);
 std::vector<long double> sums(out.size(),0);
 for(R_xlen_t i=0;i<weights.size();++i) {
  if((i&16383)==0)checkUserInterrupt();
  if(groups[i]<1 || groups[i]>ng || bins[i]<1 || bins[i]>nb)
    stop("invalid group/bin index");
  if(!std::isnan(weights[i])) sums[(groups[i]-1)+ng*(bins[i]-1)]+=weights[i];
 }
 for(R_xlen_t i=0;i<out.size();++i)out[i]=static_cast<double>(sums[i]);
 return out;
}

// [[Rcpp::export]]
NumericVector contour_lengths(List contours) {
 NumericVector out(contours.size());
 for(R_xlen_t c=0;c<contours.size();++c) {
  if((c&1023)==0)checkUserInterrupt();
  NumericMatrix z=contours[c];
  if(z.ncol()!=4)stop("expected four columns");
  long double sum=0;
  for(int i=0;i<z.nrow();++i) {
   const double dx=z(i,2)-z(i,0),dy=z(i,3)-z(i,1);
   // Match separate R powers/rowSums; do not contract multiply-add.
   volatile double xx=dx*dx,yy=dy*dy;
   long double squared=xx;squared+=yy;
   sum+=std::sqrt(static_cast<double>(squared));
  }
  out[c]=static_cast<double>(sum);
 }
 return out;
}

// [[Rcpp::export]]
int accumulator_bytes(){return sizeof(long double);}
