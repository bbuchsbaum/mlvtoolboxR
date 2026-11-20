// Geometry kernels: Bradley Buchsbaum, 2026, GPL-3.0-or-later.
// AOF formulation: Morteza Rezanejad (2019), MLV Toolbox.
// Pratt lookup tables are credited separately in pratt_lut.h.
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <set>
#include <map>
#include <array>
#include "pratt_lut.h"
using namespace Rcpp;
static double cross(double x,double y,double a,double b){return x*b-y*a;}
// Segment extent used to convert the absolute gap into a parameter fraction:
// the larger coordinate difference (upstream MLV; one raster step, so up to
// sqrt(2) px along diagonals) or the Euclidean length (rotation invariant).
static double gap_norm(double x,double y,bool euclidean){
 return euclidean?std::hypot(x,y):std::max(std::abs(x),std::abs(y));
}

// [[Rcpp::export]]
NumericMatrix cpp_intersections(NumericMatrix z,IntegerVector contour,
                                double ae,double re,bool sweep,bool euclidean=false) {
 int n=z.nrow(); std::vector<double> lo(n),hi(n),ly(n),hy(n);
 std::vector<int> order(n); std::iota(order.begin(),order.end(),0);
 for(int i=0;i<n;++i){
  double dx=z(i,2)-z(i,0),dy=z(i,3)-z(i,1);
  double t=std::min(re,ae/gap_norm(dx,dy,euclidean));
  // The candidate box is only a filter; a relative pad keeps it conservative
  // against rounding of the extension, so sweep finds every pair exhaustive
  // search accepts.
  double pad=1e-9*(1+std::max(std::max(std::abs(z(i,0)),std::abs(z(i,2))),
                                std::max(std::abs(z(i,1)),std::abs(z(i,3)))));
  lo[i]=std::min(z(i,0),z(i,2))-std::abs(dx)*t-pad;
  hi[i]=std::max(z(i,0),z(i,2))+std::abs(dx)*t+pad;
  ly[i]=std::min(z(i,1),z(i,3))-std::abs(dy)*t-pad;
  hy[i]=std::max(z(i,1),z(i,3))+std::abs(dy)*t+pad;
 }
 if(sweep)std::stable_sort(order.begin(),order.end(),[&](int i,int j){return lo[i]<lo[j];});
 std::vector<std::array<double,4>> result; std::vector<int> active;
 // Canonical frame of a segment: its lexicographically smaller endpoint Q,
 // the other endpoint R and direction R - Q. It does not depend on the stored
 // direction, and pairs are evaluated in geometric order, so the accepted
 // set and the reported point are identical under reversal, renumbering and
 // either search method (MLV-R015).
 auto lexless=[&](int k){return z(k,0)<z(k,2)||(z(k,0)==z(k,2)&&z(k,1)<=z(k,3));};
 auto check=[&](int i,int j){
  if(contour[i]==contour[j])return;
  if(sweep && (hy[i]<ly[j]||hy[j]<ly[i]))return;
  double f[2][4];
  for(int s2=0;s2<2;++s2){int k=s2?j:i;bool l=lexless(k);
   f[s2][0]=l?z(k,0):z(k,2);f[s2][1]=l?z(k,1):z(k,3);f[s2][2]=l?z(k,2):z(k,0);f[s2][3]=l?z(k,3):z(k,1);}
  int A=0,B=1;
  for(int c=0;c<4;++c){if(f[0][c]<f[1][c])break;if(f[0][c]>f[1][c]){A=1;B=0;break;}if(c==3&&j<i){A=1;B=0;}}
  const double *qa=f[A],*qb=f[B];
  double ax=qa[2]-qa[0],ay=qa[3]-qa[1],bx=qb[2]-qb[0],by=qb[3]-qb[1];
  double den=cross(ax,ay,bx,by);
  double scale=std::hypot(ax,ay)*std::hypot(bx,by);
  if(!std::isfinite(den)||std::abs(den)<=32*std::numeric_limits<double>::epsilon()*scale)return;
  // Parameters measured from each endpoint toward the other, so an
  // extension beyond either end is tested by the same arithmetic; an
  // extension exactly equal to the limit is accepted.
  double ta=std::min(re,ae/gap_norm(ax,ay,euclidean));
  double tb=std::min(re,ae/gap_norm(bx,by,euclidean));
  double t0=cross(qb[0]-qa[0],qb[1]-qa[1],bx,by)/den;
  double t1=cross(qb[0]-qa[2],qb[1]-qa[3],bx,by)/(-den);
  double u0=cross(qb[0]-qa[0],qb[1]-qa[1],ax,ay)/den;
  double u1=cross(qa[0]-qb[2],qa[1]-qb[3],ax,ay)/den;
  if(!(t0>=-ta && t1>=-ta && u0>=-tb && u1>=-tb))return;
  int a=std::min(i,j),b=std::max(i,j);
  result.push_back({double(a+1),double(b+1),qa[0]+t0*ax,qa[1]+t0*ay});
 };
 for(int q=0;q<n;++q){
  if((q&1023)==0)checkUserInterrupt();
  int i=order[q];
  if(sweep){
   active.erase(std::remove_if(active.begin(),active.end(),[&](int j){return hi[j]<lo[i];}),active.end());
   for(int j:active)check(i,j);
   active.push_back(i);
  }else for(int j=0;j<i;++j)check(i,j);
 }
 std::sort(result.begin(),result.end(),[](auto a,auto b){return a[0]<b[0]||(a[0]==b[0]&&a[1]<b[1]);});
 NumericMatrix out(result.size(),4);
 for(size_t i=0;i<result.size();++i)for(int k=0;k<4;++k)out(i,k)=result[i][k];
 return out;
}
// [[Rcpp::export]]
IntegerVector cpp_cluster_points(NumericMatrix p,double radius){
 int n=p.nrow(); std::vector<int> parent(n),order(n);
 std::iota(parent.begin(),parent.end(),0);std::iota(order.begin(),order.end(),0);
 auto find=[&](int x){int r=x;while(parent[r]!=r)r=parent[r];while(parent[x]!=x){int q=parent[x];parent[x]=r;x=q;}return r;};
 std::stable_sort(order.begin(),order.end(),[&](int i,int j){return p(i,0)<p(j,0);});
 for(int q=0;q<n;++q){if((q&1023)==0)checkUserInterrupt();int i=order[q];
  for(int k=q-1;k>=0;--k){int j=order[k];if(p(i,0)-p(j,0)>radius)break;
   if(std::hypot(p(i,0)-p(j,0),p(i,1)-p(j,1))<=radius){int a=find(i),b=find(j);parent[std::max(a,b)]=std::min(a,b);}
  }
 }
 IntegerVector out(n);for(int i=0;i<n;++i)out[i]=find(i)+1;return out;
}
// [[Rcpp::export]]
List cpp_distance_transform(LogicalMatrix features){
 int h=features.nrow(),w=features.ncol();double inf=R_PosInf;
 NumericMatrix first(h,w),d(h,w);IntegerMatrix nearest(h,w),idx(h,w);
 for(int c=0;c<w;++c){int last=-1;
  for(int r=0;r<h;++r){if(features(r,c))last=r;
   first(r,c)=last<0?inf:double(r-last)*(r-last);nearest(r,c)=last;}
  last=-1;
  for(int r=h-1;r>=0;--r){if(features(r,c))last=r;
   if(last>=0 && double(last-r)*(last-r)<first(r,c)){first(r,c)=double(last-r)*(last-r);nearest(r,c)=last;}
  }
 }
 std::vector<int> v(w);std::vector<double> breaks(w+1);
 for(int r=0;r<h;++r){
  if((r&255)==0)checkUserInterrupt();int k=-1;
  for(int q=0;q<w;++q){
   if(!std::isfinite(first(r,q)))continue;
   double s=-inf;
   while(k>=0){int p=v[k];s=((first(r,q)+double(q)*q)-(first(r,p)+double(p)*p))/(2.0*(q-p));
    if(s>breaks[k])break;--k;
   }
   ++k;v[k]=q;breaks[k]=(k==0?-inf:s);breaks[k+1]=inf;
  }
  if(k<0){for(int c=0;c<w;++c){d(r,c)=inf;idx(r,c)=NA_INTEGER;}continue;}
  int t=0;for(int c=0;c<w;++c){
   while(t<k && breaks[t+1]<c)++t;
   int q=v[t];double dd=double(c-q)*(c-q)+first(r,q);
   d(r,c)=std::sqrt(dd);idx(r,c)=nearest(r,q)+h*q+1;
  }
 }
 return List::create(_["distance"]=d,_["index"]=idx);
}
// [[Rcpp::export]]
NumericMatrix cpp_aof(NumericMatrix distance,IntegerMatrix idx,int samples,double epsilon){
 int h=distance.nrow(),w=distance.ncol();NumericMatrix out(h,w);
 std::vector<double> nx(samples),ny(samples);
 const double alpha=(2*M_PI)/samples;
 for(int s=0;s<samples;++s){nx[s]=std::cos(alpha*s);ny[s]=std::sin(alpha*s);}
 for(int c=0;c<w;++c){if((c&31)==0)checkUserInterrupt();for(int r=0;r<h;++r){
  if(!std::isfinite(distance(r,c))||distance(r,c)<=-1.5 ||
     r+1<=epsilon || c+1<=epsilon || r+1>=h-epsilon || c+1>=w-epsilon)continue;
  double flux=0;
  for(int s=0;s<samples;++s){
   double px=(r+1)+epsilon*nx[s]+0.5,py=(c+1)+epsilon*ny[s]+0.5;
   int rr=int(std::floor(px))-1,cc=int(std::floor(py))-1;
   if(rr<0||rr>=h||cc<0||cc>=w)continue;
   int label=idx(rr,cc);if(label==NA_INTEGER||label<1||label>h*w)continue;
   int b=label-1;double dx=(b%h+1)-px,dy=(b/h+1)-py,len=std::hypot(dx,dy);
   if(len>0)flux+=(dx*nx[s]+dy*ny[s])/len;
  }
  out(r,c)=flux;
 }}return out;
}
static std::vector<unsigned char> lut_apply(const std::vector<unsigned char>& x,int h,int w,const unsigned char* lut){
 std::vector<unsigned char> out(x.size());
 for(int c=0;c<w;++c)for(int r=0;r<h;++r){
  int code=0,bit=8;
  for(int dc=-1;dc<=1;++dc)for(int dr=-1;dr<=1;++dr,--bit){
   int rr=r+dr,cc=c+dc;if(rr>=0&&rr<h&&cc>=0&&cc<w&&x[rr+h*cc])code|=1<<bit;
  }
  out[r+h*c]=lut[code];
 }
 return out;
}
// [[Rcpp::export]]
LogicalMatrix cpp_skeletonize(LogicalMatrix mask){
 int h=mask.nrow(),w=mask.ncol();std::vector<unsigned char> x(mask.size());
 for(R_xlen_t i=0;i<mask.size();++i)x[i]=mask[i]?1:0;
 for(int iteration=0;;++iteration){
  if((iteration&15)==0)checkUserInterrupt();
  auto a=lut_apply(x,h,w,pratt_conditional);auto b=lut_apply(a,h,w,pratt_unconditional);
  bool change=false;for(size_t i=0;i<x.size();++i){unsigned char next=x[i]&b[i];change|=next!=x[i];x[i]=next;}
  if(!change)break;
 }
 x=lut_apply(x,h,w,pratt_bridge);LogicalMatrix out(h,w);
 for(size_t i=0;i<x.size();++i)out[i]=x[i];return out;
}
// [[Rcpp::export]]
LogicalMatrix cpp_remove_small(LogicalMatrix mask,int minimum){
 int h=mask.nrow(),w=mask.ncol();LogicalMatrix out=clone(mask);
 std::vector<unsigned char> seen(mask.size());std::vector<int> queue;
 for(int i=0;i<mask.size();++i){
  if((i&65535)==0)checkUserInterrupt();
  if(!mask[i]||seen[i])continue;queue.clear();queue.push_back(i);seen[i]=1;
  for(size_t q=0;q<queue.size();++q){int p=queue[q],r=p%h,c=p/h;
   for(int dc=-1;dc<=1;++dc)for(int dr=-1;dr<=1;++dr){
    int rr=r+dr,cc=c+dc;if(rr<0||rr>=h||cc<0||cc>=w)continue;
    int j=rr+h*cc;if(mask[j]&&!seen[j]){seen[j]=1;queue.push_back(j);}
   }
  }if(int(queue.size())<minimum)for(int j:queue)out[j]=false;
 }return out;
}
// [[Rcpp::export]]
IntegerVector cpp_raster_indices(NumericMatrix seg,int h,int w,double width){
 std::vector<int> hits;
 auto plot=[&](int x,int y){if(x>=0&&x<w&&y>=0&&y<h)hits.push_back(y+h*x+1);};
 for(int i=0;i<seg.nrow();++i){
  if((i&1023)==0)checkUserInterrupt();
  double x0=seg(i,0)-1,y0=seg(i,1)-1,x1=seg(i,2)-1,y1=seg(i,3)-1;
  if(width>1){
   int minx=int(std::max(0.0,std::ceil(std::min(x0,x1)-width/2))),maxx=int(std::min(double(w-1),std::floor(std::max(x0,x1)+width/2)));
   int miny=int(std::max(0.0,std::ceil(std::min(y0,y1)-width/2))),maxy=int(std::min(double(h-1),std::floor(std::max(y0,y1)+width/2)));
   double dx=x1-x0,dy=y1-y0,den=dx*dx+dy*dy;
   for(int x=minx;x<=maxx;++x)for(int y=miny;y<=maxy;++y){double t=den?std::max(0.0,std::min(1.0,((x-x0)*dx+(y-y0)*dy)/den)):0;
    if(std::hypot(x-(x0+t*dx),y-(y0+t*dy))<=width/2)plot(x,y);}
  }else{
   double dx=x1-x0,dy=y1-y0,enter=0,leave=1;
   bool ok=true;double pp[4]={-dx,dx,-dy,dy},qq[4]={x0,double(w-1)-x0,y0,double(h-1)-y0};
   for(int k=0;k<4;++k){if(pp[k]==0){if(qq[k]<0)ok=false;}else{double t=qq[k]/pp[k];if(pp[k]<0)enter=std::max(enter,t);else leave=std::min(leave,t);}}
   if(!ok||enter>leave)continue;
   int ax=int(std::round(x0+enter*dx)),ay=int(std::round(y0+enter*dy)),bx=int(std::round(x0+leave*dx)),by=int(std::round(y0+leave*dy));
   // Major-axis DDA between the rounded endpoints. The minor coordinate is the
   // nearest pixel to the exact line; an exact half-pixel tie goes to the larger
   // minor coordinate for positive slope and the smaller for negative slope.
   // This is invariant to segment direction (Bresenham's error-term ties are
   // not) and reproduces MATLAB insertShape('Line',LineWidth=1,SmoothEdges=false)
   // on the Forest_DBW reference raster (MLV-R010).
   long long mx=bx-ax,my=by-ay;bool up=(mx>0&&my>0)||(mx<0&&my<0);
   auto minor=[&](long long k,long long num,long long den,int base){
    // nearest integer to base+k*num/den with den>0; exact ties resolved by slope sign
    long long p=k*num,q=den;long long fl=p>=0?p/q:-((-p+q-1)/q);long long r2=2*(p-fl*q);
    if(r2<q)return base+int(fl); if(r2>q)return base+int(fl)+1; return base+int(up?fl+1:fl);
   };
   if(std::llabs(mx)>=std::llabs(my)){
    long long n=std::llabs(mx),s=mx>=0?1:-1;
    for(long long k=0;k<=n;++k)plot(ax+int(s*k),n?minor(k,my,n,ay):ay);
   }else{
    long long n=std::llabs(my),s=my>=0?1:-1;
    for(long long k=0;k<=n;++k)plot(minor(k,mx,n,ax),ay+int(s*k));
   }
  }
 }
 std::sort(hits.begin(),hits.end());hits.erase(std::unique(hits.begin(),hits.end()),hits.end());return wrap(hits);
}
// [[Rcpp::export]]
List cpp_trace_graph(LogicalMatrix mask,int minimum){
 int h=mask.nrow(),w=mask.ncol(),n=mask.size();
 std::vector<std::vector<int>> adjacent(n);
 for(int i=0;i<n;++i)if(mask[i]){int r=i%h,c=i/h;
  for(int dc=-1;dc<=1;++dc)for(int dr=-1;dr<=1;++dr){
   if(!dr&&!dc)continue;int rr=r+dr,cc=c+dc;if(rr<0||rr>=h||cc<0||cc>=w||!mask[rr+h*cc])continue;
   if(dr&&dc && (mask[r+h*cc]||mask[rr+h*c]))continue;
   adjacent[i].push_back(rr+h*cc);
  }
  std::sort(adjacent[i].begin(),adjacent[i].end());
 }
 std::set<std::pair<int,int>> visited;std::vector<std::vector<int>> paths;
 auto edge=[](int a,int b){return std::make_pair(std::min(a,b),std::max(a,b));};
 auto trace=[&](int a,int b){
  std::vector<int> path{a};int previous=a,current=b;visited.insert(edge(a,b));
  for(;;){path.push_back(current);if(adjacent[current].size()!=2)break;
   int next=adjacent[current][0]==previous?adjacent[current][1]:adjacent[current][0];
   if(visited.count(edge(current,next)))break;
   visited.insert(edge(current,next));previous=current;current=next;
  }
  if(int(path.size())>=minimum){for(int &p:path)++p;paths.push_back(path);}
 };
 for(int i=0;i<n;++i)if(mask[i] && adjacent[i].size()!=2)for(int j:adjacent[i])if(!visited.count(edge(i,j)))trace(i,j);
 for(int i=0;i<n;++i)if(mask[i]){if((i&65535)==0)checkUserInterrupt();for(int j:adjacent[i])if(!visited.count(edge(i,j)))trace(i,j);}
 List out(paths.size());for(size_t i=0;i<paths.size();++i)out[i]=wrap(paths[i]);return out;
}
struct KDNode{int point,left,right,axis;};
class KDTree {
 const NumericMatrix& p;std::vector<int> order;std::vector<KDNode> nodes;
 int build(int a,int b,int depth){
  if(a>=b)return -1;int axis=depth%2,m=(a+b)/2;
  std::nth_element(order.begin()+a,order.begin()+m,order.begin()+b,[&](int i,int j){return p(i,axis)<p(j,axis)||(p(i,axis)==p(j,axis)&&i<j);});
  int node=nodes.size();nodes.push_back({order[m],-1,-1,axis});
  int l=build(a,m,depth+1),r=build(m+1,b,depth+1);nodes[node].left=l;nodes[node].right=r;return node;
 }
 void search(int node,double x,double y,int& best,double& dist)const{
  if(node<0)return;const KDNode& nd=nodes[node];int i=nd.point;
  double dx=x-p(i,0),dy=y-p(i,1),d=dx*dx+dy*dy;
  if(d<dist||(d==dist&&(best<0||i<best))){dist=d;best=i;}
  double delta=nd.axis?dy:dx;int first=delta<0?nd.left:nd.right,second=delta<0?nd.right:nd.left;
  search(first,x,y,best,dist);if(delta*delta<=dist)search(second,x,y,best,dist);
 }
 void range(int node,double x,double y,double r2,std::vector<int>& hits)const{
  if(node<0)return;const KDNode& nd=nodes[node];int i=nd.point;
  double dx=x-p(i,0),dy=y-p(i,1);if(dx*dx+dy*dy<=r2)hits.push_back(i);
  double delta=nd.axis?dy:dx;int first=delta<0?nd.left:nd.right,second=delta<0?nd.right:nd.left;
  range(first,x,y,r2,hits);if(delta*delta<=r2)range(second,x,y,r2,hits);
 }
 public:
 KDTree(const NumericMatrix& input):p(input),order(input.nrow()){std::iota(order.begin(),order.end(),0);build(0,order.size(),0);}
 void within(double x,double y,double r,std::vector<int>& hits)const{hits.clear();if(!nodes.empty())range(0,x,y,r*r,hits);}
 void nearest(double x,double y,int& best,double& d)const{best=-1;d=R_PosInf;if(!nodes.empty())search(0,x,y,best,d);}
};
// [[Rcpp::export]]
List cpp_nearest(NumericMatrix points,NumericMatrix query){
 KDTree tree(points);IntegerVector idx(query.nrow());NumericVector distance(query.nrow());
 for(int i=0;i<query.nrow();++i){if((i&4095)==0)checkUserInterrupt();int best;double d;tree.nearest(query(i,0),query(i,1),best,d);idx[i]=best<0?NA_INTEGER:best+1;distance[i]=std::sqrt(d);}
 return List::create(_["index"]=idx,_["distance"]=distance);
}
// [[Rcpp::export]]
List cpp_fit_lines(NumericMatrix xy,double tolerance){
 int n=xy.nrow();std::vector<std::array<int,2>> spans;NumericVector scores(n);
 int start=0;
 while(start<n-1){
  int end=start+1;
  for(int t=start+2;t<n;++t){
   if((t&1023)==0)checkUserInterrupt();
   double dx=xy(t,0)-xy(start,0),dy=xy(t,1)-xy(start,1),len=std::hypot(dx,dy),sum=0;
   for(int j=start+1;j<t;++j)sum+=len>0?std::abs(cross(xy(j,0)-xy(start,0),xy(j,1)-xy(start,1),dx,dy))/len:std::hypot(xy(j,0)-xy(start,0),xy(j,1)-xy(start,1));
   if(sum/(t-start-1)>tolerance)break;end=t;
  }
  spans.push_back({start+1,end+1});for(int j=start;j<=end;++j)scores[j]=1.0-1.0/(end-start+1);start=end;
 }
 IntegerMatrix out(spans.size(),2);for(size_t i=0;i<spans.size();++i){out(i,0)=spans[i][0];out(i,1)=spans[i][1];}
 return List::create(_["spans"]=out,_["scores"]=scores);
}
// Sequential eight-direction hit-or-miss thinning (Serra; Gonzalez & Woods,
// "thinning" with rotated edge/corner structuring elements). Each direction
// deletes, in parallel, every pixel whose 3x3 neighbourhood matches its element;
// cycles repeat until no pixel changes. The result is one pixel wide and
// homotopic. Direction order (background side): S, SW, E, SE, W, NW, NE, N.
// This order was selected by exhaustive search over all 8! orders: it reproduces
// MATLAB bwmorph(...,'skel',Inf) exactly on the Forest_DBW reference flux mask
// (24,904 pixels) and on its ink-thinning stage. MATLAB's algorithm is
// undocumented; equivalence beyond this evidence is not claimed. Replaces the
// non-thin Octave Pratt default (MLV-R013).
// [[Rcpp::export]]
LogicalMatrix cpp_skeletonize_hmt(LogicalMatrix mask){
 int h=mask.nrow(),w=mask.ncol();std::vector<unsigned char> x(mask.size());
 for(R_xlen_t i=0;i<mask.size();++i)x[i]=mask[i]?1:0;
 // 1 foreground, 0 background, -1 don't care; [dr+1][dc+1]
 int se[8][3][3]={{{0,0,0},{-1,1,-1},{1,1,1}},{{-1,0,0},{1,1,0},{-1,1,-1}}};
 for(int k=2;k<8;++k)for(int i=0;i<3;++i)for(int j=0;j<3;++j)se[k][i][j]=se[k-2][2-j][i];
 const int order[8]={4,5,2,3,6,7,1,0};
 std::vector<int> del;
 for(int iteration=0;;++iteration){
  if((iteration&15)==0)checkUserInterrupt();
  bool change=false;
  for(int q=0;q<8;++q){const int (*b)[3]=se[order[q]];del.clear();
   for(int c=0;c<w;++c)for(int r=0;r<h;++r){if(!x[r+h*c])continue;bool hit=true;
    for(int dr=-1;dr<=1&&hit;++dr)for(int dc=-1;dc<=1;++dc){int want=b[dr+1][dc+1];if(want<0)continue;
     int rr=r+dr,cc=c+dc;int v=(rr>=0&&rr<h&&cc>=0&&cc<w)?x[rr+h*cc]:0;if(v!=want){hit=false;break;}}
    if(hit)del.push_back(r+h*c);}
   for(int i:del)x[i]=0;if(!del.empty())change=true;}
  if(!change)break;
 }
 LogicalMatrix out(h,w);for(size_t i=0;i<x.size();++i)out[i]=x[i];return out;
}
// Tie-invariant, pixel-centred AOF. Each sample point p=(i,j)+epsilon*(cos,sin)
// takes the mean unit vector toward every equidistant nearest boundary pixel of
// the nearest pixel cell (of both cells when p lies on a cell boundary). The
// result is equivariant under the eight lattice symmetries, unlike the MLV
// sampling (circle centred at (i+0.5,j+0.5), labels from floor(p), one label; MLV-M027).
// [[Rcpp::export]]
NumericMatrix cpp_aof_centered(NumericMatrix distance,LogicalMatrix ink,LogicalMatrix outer,int samples,double epsilon){
 int h=distance.nrow(),w=distance.ncol();R_xlen_t n=(R_xlen_t)h*w;
 std::vector<R_xlen_t> ptr(n+1);std::vector<int> lab;
 for(R_xlen_t i=0;i<n;++i){
  if((i&65535)==0)checkUserInterrupt();
  ptr[i]=lab.size();int r=i%h,c=i/h;
  if(outer[i]){lab.push_back(int(i));continue;}
  double D=distance[i];if(!std::isfinite(D))continue;
  long long S=llround(D*D);int d=int(std::floor(std::sqrt(double(S))))+1;bool on_ink=ink[i];
  std::vector<int> hits;
  for(int dr=-d;dr<=d;++dr){long long rem=S-(long long)dr*dr;if(rem<0)continue;
   long long dc=llround(std::sqrt(double(rem)));if(dc*dc!=rem)continue;
   for(int sg=-1;sg<=1;sg+=2){if(dc==0&&sg>0)continue;int rr=r+dr,cc=c+sg*int(dc);
    if(rr<0||rr>=h||cc<0||cc>=w)continue;R_xlen_t j=rr+(R_xlen_t)h*cc;
    bool feature=on_ink?(!ink[j]||outer[j]):bool(ink[j]);if(feature)hits.push_back(int(j));}}
  std::sort(hits.begin(),hits.end());lab.insert(lab.end(),hits.begin(),hits.end());
 }
 ptr[n]=lab.size();
 NumericMatrix out(h,w);std::vector<double> nx(samples),ny(samples);
 const double alpha=(2*M_PI)/samples,tol=1e-9;
 for(int s=0;s<samples;++s){nx[s]=std::cos(alpha*s);ny[s]=std::sin(alpha*s);}
 std::vector<int> labels;
 for(int c=0;c<w;++c){if((c&31)==0)checkUserInterrupt();for(int r=0;r<h;++r){
  // Symmetric border rule: every sample point must lie within [1,h]x[1,w] (MLV uses
  // i > epsilon and i < h - epsilon, which skips one more row at the bottom than at the top).
  if(!std::isfinite(distance(r,c))||distance(r,c)<=-1.5||r+1<1+epsilon||c+1<1+epsilon||r+1>h-epsilon||c+1>w-epsilon)continue;
  double flux=0;
  for(int s=0;s<samples;++s){
   double px=(r+1)+epsilon*nx[s],py=(c+1)+epsilon*ny[s];int ci[2],cj[2],ni=1,nj=1;
   double fx=std::floor(px),fy=std::floor(py);
   if(std::abs(px-fx-.5)<tol){ci[0]=int(fx);ci[1]=int(fx)+1;ni=2;}else ci[0]=int(std::floor(px+.5));
   if(std::abs(py-fy-.5)<tol){cj[0]=int(fy);cj[1]=int(fy)+1;nj=2;}else cj[0]=int(std::floor(py+.5));
   labels.clear();
   for(int a=0;a<ni;++a)for(int b=0;b<nj;++b){int rr=ci[a]-1,cc=cj[b]-1;if(rr<0||rr>=h||cc<0||cc>=w)continue;
    R_xlen_t k=rr+(R_xlen_t)h*cc;for(R_xlen_t t=ptr[k];t<ptr[k+1];++t)labels.push_back(lab[t]);}
   std::sort(labels.begin(),labels.end());labels.erase(std::unique(labels.begin(),labels.end()),labels.end());
   if(labels.empty())continue;double vx=0,vy=0;
   for(int L:labels){double dx=(L%h+1)-px,dy=(L/h+1)-py,len=std::hypot(dx,dy);if(len>0){vx+=dx/len;vy+=dy/len;}}
   flux+=(vx*nx[s]+vy*ny[s])/labels.size();
  }
  out(r,c)=flux;
 }}
 return out;
}
// Published-compatible branch tracing: a line-by-line port of Kovesi's edgelink
// (trackedge/nextpoint/findendsjunctions) as called by MLV traceSkeleton/GetConSeg,
// before cleanedgelist (which is ported in R). The input must already be thinned
// with MATLAB-equivalent morphology. Faithful details include the raster scan in
// row-major order, the first edge labelled -0 (= background), junction lookup
// that ignores labels of other edges, the "remember last" fallback, and loop
// closure that repeats the first point. Returns 1-based (row, col) point lists.
// [[Rcpp::export]]
List cpp_edgelink_trace(LogicalMatrix mask){
 const int R=mask.nrow(),C=mask.ncol();
 std::vector<int> E((size_t)R*C,0);
 auto at=[&](int r,int c)->int&{return E[(size_t)r+(size_t)R*c];};
 // bwmorph(...,'clean'): remove isolated pixels
 for(int c=0;c<C;++c)for(int r=0;r<R;++r)if(mask(r,c)){
  bool nb=false;for(int dc=-1;dc<=1&&!nb;++dc)for(int dr=-1;dr<=1;++dr){if(!dr&&!dc)continue;int rr=r+dr,cc=c+dc;
   if(rr>=0&&rr<R&&cc>=0&&cc<C&&mask(rr,cc)){nb=true;break;}}
  at(r,c)=nb?1:0;}
 // junctions: crossings >= 6 around the 8-neighbourhood (zero padded)
 std::vector<unsigned char> J((size_t)R*C,0);
 const int ar[8]={-1,0,1,1,1,0,-1,-1},ac[8]={-1,-1,-1,0,1,1,1,0};
 for(int c=0;c<C;++c)for(int r=0;r<R;++r)if(at(r,c)==1){
  int v[8];for(int k=0;k<8;++k){int rr=r+ar[k],cc=c+ac[k];v[k]=(rr>=0&&rr<R&&cc>=0&&cc<C&&at(rr,cc)==1)?1:0;}
  int cr=0;for(int k=0;k<8;++k)cr+=std::abs(v[k]-v[(k+1)%8]);
  if(cr>=6)J[(size_t)r+(size_t)R*c]=1;}
 auto isj=[&](int r,int c){return J[(size_t)r+(size_t)R*c]!=0;};
 const int nr[8]={-1,0,1,0,-1,-1,1,1},nc[8]={0,1,0,-1,-1,1,1,-1};
 const int br[8]={-1,0,1,1,1,0,-1,-1},bc[8]={-1,-1,-1,0,1,1,1,0};
 auto nextpoint=[&](int rp,int cp,int edgeNo,int& outr,int& outc)->int{
  for(int k=0;k<8;++k){int r=rp+nr[k],c=cp+nc[k];if(r<0||r>=R||c<0||c>=C)continue;
   if(isj(r,c)&&at(r,c)!=-edgeNo){outr=r;outc=c;return 2;}}
  bool flag=false;int mr=0,mc=0;
  for(int k=0;k<8;++k){int r=rp+nr[k],c=cp+nc[k];if(r<0||r>=R||c<0||c>=C)continue;
   if(at(r,c)==1){int cnt=0;
    for(int q=0;q<8;++q){int rr=r+br[q],cc=c+bc[q];int val=(rr>=0&&rr<R&&cc>=0&&cc<C)?at(rr,cc):0;if(val==-edgeNo)++cnt;}
    if(cnt<2){outr=r;outc=c;return 1;}
    flag=true;mr=r;mc=c;}}
  if(flag){outr=mr;outc=mc;return 1;}
  return 0;};
 std::vector<IntegerMatrix> edges;int edgeNo=0;
 for(int r=0;r<R;++r){if((r&63)==0)checkUserInterrupt();for(int c=0;c<C;++c)if(at(r,c)==1){
  std::vector<std::array<int,2>> pts{{r,c}};at(r,c)=-edgeNo;int qr,qc;
  int st=nextpoint(r,c,edgeNo,qr,qc);
  while(st!=0){pts.push_back({qr,qc});at(qr,qc)=-edgeNo;if(st==2)st=0;else st=nextpoint(qr,qc,edgeNo,qr,qc);}
  if(!isj(r,c)){std::reverse(pts.begin(),pts.end());st=nextpoint(r,c,edgeNo,qr,qc);
   while(st!=0){pts.push_back({qr,qc});at(qr,qc)=-edgeNo;if(st==2)st=0;else st=nextpoint(qr,qc,edgeNo,qr,qc);}}
  if(pts.size()>=4&&std::abs(pts.front()[0]-pts.back()[0])<=1&&std::abs(pts.front()[1]-pts.back()[1])<=1)pts.push_back(pts.front());
  IntegerMatrix m(pts.size(),2);for(size_t i=0;i<pts.size();++i){m(i,0)=pts[i][0]+1;m(i,1)=pts[i][1]+1;}
  edges.push_back(m);++edgeNo;}}
 List out(edges.size());for(size_t i=0;i<edges.size();++i)out[i]=edges[i];return out;
}
// Line-by-line port of Kovesi's cleanedgelist (MLV vendored copy); see R/mat-compat.R.
// The node adjacency A and edge adjacency B of the original are symmetric 0/1
// matrices; they are held as ordered adjacency sets (memory O(nodes + links)
// instead of (2E)^2 bytes), with the original's row/column copy, OR-merge and
// removal semantics, its "Nconnections = sum(A)" recomputation after a merge and
// its subtract-then-zero update in removeedge. Initial node matching uses a hash
// of node coordinates.
// [[Rcpp::export]]
List cpp_cleanedgelist(List edgelist,double minlength){
 int Ne=edgelist.size(),Nn=2*Ne;
 std::vector<std::vector<std::array<int,2>>> E(Ne);std::vector<char> alive(Ne,1);
 for(int n=0;n<Ne;++n){IntegerMatrix m=edgelist[n];for(int i=0;i<m.nrow();++i)E[n].push_back({m(i,0),m(i,1)});}
 std::vector<std::set<int>> A(Nn),B(Ne);
 auto nodept=[&](int k){const auto& e=E[k/2];return k%2==0?e.front():e.back();};
 {std::map<std::pair<int,int>,std::vector<int>> at;
  for(int k=0;k<Nn;++k){auto p=nodept(k);at[{p[0],p[1]}].push_back(k);}
  for(auto& kv:at){auto& v=kv.second;for(size_t i=0;i<v.size();++i)for(size_t j=0;j<v.size();++j)if(i!=j){
   A[v[i]].insert(v[j]);B[v[i]/2].insert(v[j]/2);}}}
 std::vector<double> Nc(Nn),Ec(Ne);
 auto recount=[&](){for(int j=0;j<Nn;++j)Nc[j]=A[j].size();for(int j=0;j<Ne;++j)Ec[j]=B[j].size();};
 recount();
 auto bself=[&](int n){return B[n].count(n)>0;};
 auto empty=[&](int n){return !alive[n];};
 auto len=[&](int n){double s=0;const auto& e=E[n];for(size_t i=1;i<e.size();++i)s+=std::sqrt(double((e[i][0]-e[i-1][0])*(e[i][0]-e[i-1][0])+(e[i][1]-e[i-1][1])*(e[i][1]-e[i-1][1])));return s;};
 auto removeedge=[&](int n){
  alive[n]=0;E[n].clear();
  for(int j:B[n])Ec[j]-=1;Ec[n]=0;
  std::vector<int> bn(B[n].begin(),B[n].end());for(int j:bn)B[j].erase(n);B[n].clear();
  for(int d:{2*n,2*n+1})for(int j:A[d])Nc[j]-=1;
  for(int d:{2*n,2*n+1}){std::vector<int> an(A[d].begin(),A[d].end());for(int j:an)A[j].erase(d);A[d].clear();}
 };
 auto copyRowCol=[&](int dst,int src){ // A(dst,:)=A(src,:); A(:,dst)=A(:,src)
  std::vector<int> old(A[dst].begin(),A[dst].end());for(int j:old)A[j].erase(dst);A[dst].clear();
  A[src].erase(dst);
  for(int j:A[src])if(j!=dst){A[dst].insert(j);A[j].insert(dst);}
 };
 auto mergenodes=[&](int n1,int n2){
  int e1=n1/2,e2=n2/2;if(e1==e2)return;
  int s1=2*e1,f1=2*e1+1,s2=2*e2,f2=2*e2+1;
  if(!A[n1].count(n2))stop("Attempt to merge nodes that are not connected");
  bool flip1=(n1%2==0),flip2=(n2%2==1);
  auto x=E[e1],y=E[e2];
  if(flip1)std::reverse(x.begin(),x.end());if(flip2)std::reverse(y.begin(),y.end());
  x.insert(x.end(),y.begin(),y.end());E[e1]=x;
  if(!flip1&&!flip2){copyRowCol(f1,f2);}
  else if(!flip1&&flip2){copyRowCol(f1,s2);}
  else if(flip1&&!flip2){copyRowCol(s1,f1);copyRowCol(f1,f2);}
  else {copyRowCol(s1,f1);copyRowCol(f1,s2);}
  std::vector<int> b2(B[e2].begin(),B[e2].end());
  for(int j:b2){B[e1].insert(j);B[j].insert(e1);}
  B[e1].erase(e1);
  recount();
  removeedge(e2);
 };
 struct Info{double sd;int sp,sn,en;double sc,ec;};
 auto info=[&](int n)->Info{Info z{0,-1,-1,-1,0,0};if(empty(n))return z;
  z.sn=2*n;z.en=2*n+1;z.sc=Nc[z.sn];z.ec=Nc[z.en];
  if(z.sc==0&&z.ec>=1){z.sd=z.ec;z.sp=z.en;}else if(z.sc>=1&&z.ec==0){z.sd=z.sc;z.sp=z.sn;}
  return z;};
 auto findrow=[&](int node){return std::vector<int>(A[node].begin(),A[node].end());};
 for(int n=0;n<Ne;++n){if((n&63)==0)checkUserInterrupt();
  if(!bself(n)&&!empty(n)){Info z=info(n);
   if(z.sc==1){auto r=findrow(z.sn);if(r.size()!=1)stop("unexpected connection count");mergenodes(r[0],z.sn);}
   if(!empty(n)&&z.ec==1){auto r=findrow(z.en);if(r.size()!=1)stop("unexpected connection count");mergenodes(r[0],z.en);}
  }}
 if(minlength>0){
  for(int n=0;n<Ne;++n){Info z=info(n);
   if(!empty(n)&&len(n)<minlength){
    if(!Ec[n]||(Ec[n]==1&&bself(n)))removeedge(n);
    else if(z.sd==2){std::vector<int> linking(B[n].begin(),B[n].end());
     if(linking.size()==1)removeedge(n);
     else{std::vector<int> spurs{n};std::vector<double> lens{len(n)};
      for(int L:linking)if(info(L).sd){spurs.push_back(L);lens.push_back(len(L));}
      size_t im=0;for(size_t i=1;i<lens.size();++i)if(lens[i]<lens[im])im=i;
      int del=spurs[im];int sp=info(del).sp;auto m2=findrow(sp);
      if(m2.size()!=2)stop("attempt to merge other than 2 nodes");
      removeedge(del);mergenodes(m2[0],m2[1]);}}
    else if(z.sd==3)removeedge(n);}}
  for(int n=0;n<Ne;++n)if(!empty(n)&&len(n)<minlength)if(!Ec[n]||(Ec[n]==1&&bself(n)))removeedge(n);
 }
 std::vector<IntegerMatrix> out;for(int n=0;n<Ne;++n)if(alive[n]){IntegerMatrix m(E[n].size(),2);for(size_t i=0;i<E[n].size();++i){m(i,0)=E[n][i][0];m(i,1)=E[n][i][1];}out.push_back(m);}
 List o(out.size());for(size_t i=0;i<out.size();++i)o[i]=out[i];return o;
}
// Published-compatible mirror fit: port of MLV fitLineSegments.m and
// getDistanceFromLineSegment.m, including the forced terminal span (MLV-M016).
// [[Rcpp::export]]
NumericVector cpp_fit_lines_mlv(NumericMatrix xy){
 int n=xy.nrow();NumericVector scores(n);
 if(n<2)return scores;
 auto dist=[&](int s,int t){ // rows s..t inclusive, 0-based
  int m=t-s+1;if(m<=2)return 0.0;
  double x1=xy(s,0),y1=xy(s,1),x2=xy(t,0),y2=xy(t,1),mx=x1-x2,my=y1-y2,total=0;
  if(my!=0){double a=1,b=-mx/my,c=-(x1+b*y1);
   for(int i=s+1;i<t;++i)total+=std::abs(a*xy(i,0)+b*xy(i,1)+c)/std::sqrt(a*a+b*b);}
  else for(int i=s+1;i<t;++i)total+=std::abs(xy(i,1)-y1);
  return total/(m-2);};
 std::vector<std::array<int,2>> segs;int start=1,t=start+2; // 1-based as in MATLAB
 while(t<=n){
  t=start+2;bool still=true;
  while(still){
   double d=dist(start-1,t-1);
   if(d>1||t>=n){still=false;if(n-t<=1)t=n+1;segs.push_back({start,t-1});start=t-1;}
   else ++t;
  }
 }
 int counter=1;
 for(auto& sg:segs){int N=sg[1]-sg[0]+1;for(int i=counter;i<=counter+N-1;++i)scores[i-1]=1.0-1.0/N;counter+=N-1;}
 return scores;
}

// Order-invariant nearest-sample mapping: for each query, the minimum distance to
// the sample points and, per value column, the maximum finite value over every
// sample within tol of that minimum (exact ties up to floating rounding). NA when
// no tied sample has a finite value.
// [[Rcpp::export]]
List cpp_nearest_max(NumericMatrix points,NumericMatrix query,NumericMatrix values,double tol){
 KDTree tree(points);int nq=query.nrow(),P=values.ncol();
 NumericVector distance(nq);NumericMatrix out(nq,P);std::fill(out.begin(),out.end(),NA_REAL);
 std::vector<int> hits;
 for(int i=0;i<nq;++i){if((i&4095)==0)checkUserInterrupt();int best;double d;tree.nearest(query(i,0),query(i,1),best,d);
  if(best<0){distance[i]=R_PosInf;continue;}double dm=std::sqrt(d);distance[i]=dm;
  tree.within(query(i,0),query(i,1),dm+tol,hits);
  for(int q=0;q<P;++q){double m=NA_REAL;bool any=false;for(int h:hits){double v=values(h,q);if(std::isfinite(v)){if(!any||v>m)m=v;any=true;}}out(i,q)=m;}
 }
 return List::create(_["distance"]=distance,_["value"]=out);
}
// Minimal skeleton: hmt8 thinning followed by sequential removal (column-major
// scan, repeated to a fixed point) of every 8-simple pixel (Yokoi 8-connectivity
// number 1) that has at least two foreground neighbours. No such pixel remains;
// a 2x2 block survives only where all four pixels are non-simple. Topology is
// preserved (only simple pixels are removed), end points are kept; idempotent.
static int yokoi8(const std::vector<unsigned char>& a,int h,int w,int r,int c,int& nb){
 auto at=[&](int rr,int cc){return (rr>=0&&rr<h&&cc>=0&&cc<w)?int(a[rr+h*cc]):0;};
 int x[10];x[1]=at(r,c+1);x[2]=at(r-1,c+1);x[3]=at(r-1,c);x[4]=at(r-1,c-1);x[5]=at(r,c-1);
 x[6]=at(r+1,c-1);x[7]=at(r+1,c);x[8]=at(r+1,c+1);x[9]=x[1];
 nb=0;for(int i=1;i<=8;++i)nb+=x[i];
 int b[12];for(int i=1;i<=9;++i)b[i]=1-x[i];b[10]=b[2];b[11]=b[3];
 int N=0;for(int q=1;q<=7;q+=2)N+=b[q]-b[q]*b[q+1]*b[q+2];
 return N;
}
// [[Rcpp::export]]
LogicalMatrix cpp_skeletonize_minimal(LogicalMatrix mask){
 int h=mask.nrow(),w=mask.ncol();
 LogicalMatrix t=cpp_skeletonize_hmt(mask);
 std::vector<unsigned char> a(t.size());for(R_xlen_t i=0;i<t.size();++i)a[i]=t[i];
 for(bool change=true;change;){change=false;checkUserInterrupt();
  for(int c=0;c<w;++c)for(int r=0;r<h;++r)if(a[r+h*c]){int nb;int N=yokoi8(a,h,w,r,c,nb);
   if(N==1&&nb>=2){a[r+h*c]=0;change=true;}}}
 LogicalMatrix out(h,w);for(size_t i=0;i<a.size();++i)out[i]=a[i];return out;
}
// [[Rcpp::export]]
int cpp_redundant_pixels(LogicalMatrix m){
 int h=m.nrow(),w=m.ncol(),k=0;std::vector<unsigned char> a(m.size());for(R_xlen_t i=0;i<m.size();++i)a[i]=m[i];
 for(int c=0;c<w;++c)for(int r=0;r<h;++r)if(a[r+h*c]){int nb;if(yokoi8(a,h,w,r,c,nb)==1&&nb>=2)++k;}
 return k;
}
// Mirror scores of an isolated closed loop for the direction-symmetric rule: the
// mlv2023 sample estimator (greedy line fit with mean perpendicular tolerance,
// moving mean of width three, exponent) averaged over every start point and both
// directions of the closed path (points 1..m, repeating point 1 at the end).
// [[Rcpp::export]]
NumericVector cpp_mirror_loop(NumericMatrix xy,double power,double tolerance){
 int m=xy.nrow();NumericVector acc(m);std::vector<double> comp(m,0.0); // Kahan compensation
 NumericMatrix p(m+1,2);
 for(int dir=0;dir<2;++dir)for(int k=0;k<m;++k){
  std::vector<int> ord(m);for(int i=0;i<m;++i)ord[i]=dir==0?(k+i)%m:((k-i)%m+m)%m;
  for(int i=0;i<=m;++i){int o=ord[i%m];p(i,0)=xy(o,0);p(i,1)=xy(o,1);}
  NumericVector sc=as<NumericVector>(cpp_fit_lines(p,tolerance)["scores"]);
  int n=m+1;std::vector<double> sm(n);
  for(int i=0;i<n;++i){int lo=std::max(0,i-1),hi=std::min(n-1,i+1);double s=0;for(int j=lo;j<=hi;++j)s+=sc[j];sm[i]=std::pow(s/(hi-lo+1),power);}
  sm[0]=(sm[0]+sm[m])/2;
  for(int i=0;i<m;++i){int o=ord[i];double y=sm[i]-comp[o],t=acc[o]+y;comp[o]=(t-acc[o])-y;acc[o]=t;}
 }
 for(int i=0;i<m;++i)acc[i]/=2.0*m;
 return acc;
}
