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
#include <array>
#include "pratt_lut.h"
using namespace Rcpp;
static double cross(double x,double y,double a,double b){return x*b-y*a;}

// [[Rcpp::export]]
NumericMatrix cpp_intersections(NumericMatrix z,IntegerVector contour,
                                double ae,double re,bool sweep) {
 int n=z.nrow(); std::vector<double> lo(n),hi(n),ly(n),hy(n);
 std::vector<int> order(n); std::iota(order.begin(),order.end(),0);
 for(int i=0;i<n;++i){
  double dx=z(i,2)-z(i,0),dy=z(i,3)-z(i,1);
  double t=std::min(re,ae/std::max(std::abs(dx),std::abs(dy)));
  lo[i]=std::min(z(i,0),z(i,2))-std::abs(dx)*t;
  hi[i]=std::max(z(i,0),z(i,2))+std::abs(dx)*t;
  ly[i]=std::min(z(i,1),z(i,3))-std::abs(dy)*t;
  hy[i]=std::max(z(i,1),z(i,3))+std::abs(dy)*t;
 }
 if(sweep)std::stable_sort(order.begin(),order.end(),[&](int i,int j){return lo[i]<lo[j];});
 std::vector<std::array<double,4>> result; std::vector<int> active;
 auto check=[&](int i,int j){
  if(contour[i]==contour[j])return;
  if(sweep && (hy[i]<ly[j]||hy[j]<ly[i]))return;
  double ax=z(i,2)-z(i,0),ay=z(i,3)-z(i,1),bx=z(j,2)-z(j,0),by=z(j,3)-z(j,1);
  double den=cross(ax,ay,bx,by);
  double scale=std::hypot(ax,ay)*std::hypot(bx,by);
  if(!std::isfinite(den)||std::abs(den)<=32*std::numeric_limits<double>::epsilon()*scale)return;
  double cx=z(j,0)-z(i,0),cy=z(j,1)-z(i,1);
  double t=cross(cx,cy,bx,by)/den,u=cross(cx,cy,ax,ay)/den;
  double ta=std::min(re,ae/std::max(std::abs(ax),std::abs(ay)));
  double tb=std::min(re,ae/std::max(std::abs(bx),std::abs(by)));
  if(t < -ta || t>1+ta || u < -tb || u>1+tb)return;
  int a=std::min(i,j),b=std::max(i,j);
  result.push_back({double(a+1),double(b+1),z(i,0)+t*ax,z(i,1)+t*ay});
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
  }out[r+h*c]=lut[code];
 }return out;
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
   int ddx=std::abs(bx-ax),ddy=-std::abs(by-ay),sx=ax<bx?1:-1,sy=ay<by?1:-1,error=ddx+ddy;
   for(;;){plot(ax,ay);if(ax==bx&&ay==by)break;int twice=2*error;if(twice>=ddy){error+=ddy;ax+=sx;}if(twice<=ddx){error+=ddx;ay+=sy;}}
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
  }std::sort(adjacent[i].begin(),adjacent[i].end());
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
 public:
 KDTree(const NumericMatrix& input):p(input),order(input.nrow()){std::iota(order.begin(),order.end(),0);build(0,order.size(),0);}
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
