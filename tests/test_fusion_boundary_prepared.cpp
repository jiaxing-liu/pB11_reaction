#include "fusion_boundary_segment.h"
#include <atomic>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>
#include <stdexcept>
void require(bool x){if(!x)throw std::runtime_error("prepared parity check");}
bool eq(double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;}
bool eq(const fusion_flux_boundary_value_v1&a,const fusion_flux_boundary_value_v1&b){return eq(a.polygon_distance_m,b.polygon_distance_m)&&eq(a.uncertainty_m,b.uncertainty_m)&&eq(a.angular_margin_m2,b.angular_margin_m2)&&a.classification==b.classification;}
bool eq(const fusion_boundary_segment_value_v1&a,const fusion_boundary_segment_value_v1&b){return eq(a.fraction_lo,b.fraction_lo)&&eq(a.fraction_hi,b.fraction_hi)&&eq(a.boundary_band_m,b.boundary_band_m)&&a.result==b.result&&a.tolerance_met==b.tolerance_met&&a.termination_reason==b.termination_reason&&a.intervals_visited==b.intervals_visited;}
int main(){try{
 long queries=0,segments=0;
 for(auto coefficients:std::vector<std::vector<double>>{{2,.4,.4},{2,.4,.65},{2,.4,.02,.004,.001,.0001,.6,-.01,.003,.0004,.0001}}){
  int m=(coefficients.size()-1)/2;
  for(int N:{64,512,2048}){
   auto a=coefficients;fusion_flux_boundary_prepared*b=nullptr;
   require(fusion_c_flux_boundary_prepare(m,a.data(),N,1e-12,1e-12,&b)==0&&b);
   for(int j=0;j<100;++j){
    double R=2+1.2*((j%10)/9.-.5),Z=1.5*((j/10)/9.-.5);
    fusion_flux_boundary_value_v1 old{},now{};
    int rc=fusion_c_flux_boundary(m,a.data(),N,R,Z,1e-12,1e-12,&old);
    require(fusion_c_flux_boundary_prepared_point(b,R,Z,&now)==rc&&eq(old,now));++queries;
   }
   for(double t:{0.,.31,1.57,3.14})for(double offset:{-1e-13,0.,1e-13}){
    double R=a[0],Z=0;for(int j=1;j<=m;++j){R+=a[j]*std::cos(j*t);Z+=a[m+j]*std::sin(j*t);}R+=offset;
    fusion_flux_boundary_value_v1 old{},now{};int rc=fusion_c_flux_boundary(m,a.data(),N,R,Z,1e-12,1e-12,&old);
    require(fusion_c_flux_boundary_prepared_point(b,R,Z,&now)==rc&&eq(old,now));++queries;
   }
   double start[3]={2,0,0};
   for(double endR:{2.1,2.4,2.8})for(int budget:{1,50,512}){
    double end[3]={endR,0,0};fusion_boundary_segment_value_v1 old{},now{};
    int rc=fusion_c_boundary_segment(m,a.data(),N,start,end,1e-12,1e-12,1e-5,budget,&old);
    require(fusion_c_boundary_prepared_segment(b,start,end,1e-5,budget,&now)==rc&&eq(old,now));++segments;
   }
   fusion_flux_boundary_value_v1 expected{};require(fusion_c_flux_boundary_prepared_point(b,2,0,&expected)==0);
   a.assign(a.size(),std::numeric_limits<double>::quiet_NaN());
   std::atomic<bool> okay{true};std::vector<std::thread> threads;
   for(int j=0;j<4;++j)threads.emplace_back([&]{for(int i=0;i<10;++i){fusion_flux_boundary_value_v1 v{};if(fusion_c_flux_boundary_prepared_point(b,2,0,&v)||!eq(v,expected))okay=false;}});
   for(auto&t:threads)t.join();require(okay);fusion_c_flux_boundary_destroy(b);
  }
 }
 double valid[]={2,.4,.4};fusion_flux_boundary_prepared*b=reinterpret_cast<fusion_flux_boundary_prepared*>(1);
 require(fusion_c_flux_boundary_prepare(0,valid,512,1e-12,1e-12,&b)!=0&&b==nullptr);
 require(fusion_c_flux_boundary_prepare(1,valid,15,1e-12,1e-12,&b)!=0&&b==nullptr);
 require(fusion_c_flux_boundary_prepare(1,valid,512,0,1e-12,&b)!=0&&b==nullptr);
 double bad[]={.1,.4,.4};require(fusion_c_flux_boundary_prepare(1,bad,512,1e-12,1e-12,&b)==FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY&&b==nullptr);
 valid[1]=std::numeric_limits<double>::quiet_NaN();require(fusion_c_flux_boundary_prepare(1,valid,512,1e-12,1e-12,&b)!=0&&b==nullptr);
 fusion_flux_boundary_value_v1 out{1,2,3,4};require(fusion_c_flux_boundary_prepared_point(nullptr,2,0,&out)!=0&&eq(out,fusion_flux_boundary_value_v1{}));
 double max_coeff[]={2,.4,.4};
 require(fusion_c_flux_boundary_prepare(1,max_coeff,65536,1e-12,1e-12,&b)==0);
 fusion_flux_boundary_value_v1 legacy{},cached{};
 require(fusion_c_flux_boundary(1,max_coeff,65536,2,0,1e-12,1e-12,&legacy)==0);
 require(fusion_c_flux_boundary_prepared_point(b,2,0,&cached)==0&&eq(legacy,cached));
 cached={1,2,3,4};require(fusion_c_flux_boundary_prepared_point(b,std::numeric_limits<double>::infinity(),0,&cached)!=0&&eq(cached,fusion_flux_boundary_value_v1{}));
 require(fusion_c_flux_boundary_prepared_point(b,2,0,nullptr)==PB11_STATUS_NULL_OUTPUT);
 fusion_c_flux_boundary_destroy(b);
 double reversed[]={2,.4,-.4};require(fusion_c_flux_boundary_prepare(1,reversed,512,1e-12,1e-12,&b)==FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY&&b==nullptr);
 fusion_c_flux_boundary_destroy(nullptr);
 std::cout<<"PASS point="<<queries<<" segment="<<segments<<" ownership/threads/errors\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
