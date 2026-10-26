// Research-only Boost.Polygon probe. Not linked into the R package.
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/segment_data.hpp>
#include <boost/polygon/point_data.hpp>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc,char** argv) {
  using boost::polygon::point_data;
  using boost::polygon::segment_data;
  using Diagram=boost::polygon::voronoi_diagram<double>;
  const int repeats=argc>1?std::stoi(argv[1]):1000;
  int n;if(!(std::cin>>n)||n<1||n>20000||repeats<1||repeats>100000)return 2;
  std::vector<segment_data<int>> input;input.reserve(n);
  for(int i=0;i<n;++i) {
    int a,b,c,d;if(!(std::cin>>a>>b>>c>>d)||(a==c&&b==d))return 3;
    input.emplace_back(point_data<int>(a,b),point_data<int>(c,d));
  }
  Diagram diagram;
  boost::polygon::construct_voronoi(input.begin(),input.end(),&diagram);
  auto start=std::chrono::steady_clock::now();
  for(int i=0;i<repeats;++i) {
    Diagram probe;
    boost::polygon::construct_voronoi(input.begin(),input.end(),&probe);
    if(probe.num_edges()!=diagram.num_edges())return 4;
  }
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repeats;
  std::cerr<<std::setprecision(17)<<"construct_seconds="<<seconds<<" repeats="<<repeats
           <<" cells="<<diagram.num_cells()<<" vertices="<<diagram.num_vertices()<<"\n";
  std::cout<<"site_a,site_b,category_a,category_b,primary,linear,finite,x1,y1,x2,y2\n";
  std::cout<<std::setprecision(17);
  for(const auto& edge:diagram.edges()) {
    if(!std::less<const Diagram::edge_type*>()(&edge,edge.twin()))continue;
    std::cout<<edge.cell()->source_index()+1<<','<<edge.twin()->cell()->source_index()+1<<','
             <<edge.cell()->source_category()<<','<<edge.twin()->cell()->source_category()<<','
             <<edge.is_primary()<<','<<edge.is_linear()<<','<<edge.is_finite();
    if(edge.is_finite())std::cout<<','<<edge.vertex0()->x()<<','<<edge.vertex0()->y()<<','
                                <<edge.vertex1()->x()<<','<<edge.vertex1()->y();
    else std::cout<<",NA,NA,NA,NA";
    std::cout<<'\n';
  }
}
