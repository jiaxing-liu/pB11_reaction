#include "fusion_thermal_markers.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
using M=fusion_thermal_marker_v1;
int main(){try{
 M a[4]{};double q[]={0,std::nextafter(.25,0.),.25,std::nextafter(.25,1.),std::nextafter(1.,0.)};M out[5]{};
 a[1].number_weight=1;a[3].number_weight=3;
 for(int i=0;i<4;++i){a[i].kinetic_energy_J=2;a[i].initial.x_m[0]=i;a[i].initial.u_m_s[1]=i+10;}
 M saved[4];std::memcpy(saved,a,sizeof a);
 auto run=[&](){return fusion_c_thermal_marker_resample(4,a,5,q,out);};
 check(run()==0,"weighted endpoints");int expected[]={1,1,3,3,3};long double n=0,e=0;
 for(int i=0;i<5;++i){check(out[i].initial.x_m[0]==expected[i]&&out[i].initial.u_m_s[1]==expected[i]+10,"copied state");check(out[i].number_weight==.8&&out[i].kinetic_energy_J==2,"equal weight and K");n+=out[i].number_weight;e+=out[i].number_weight*out[i].kinetic_energy_J;}
 check(std::abs(n-4)<1e-15L&&std::abs(e-8)<2e-15L,"N/E closure");check(!std::memcmp(a,saved,sizeof a),"unchanged input");M replay[5];std::memcpy(replay,out,sizeof out);check(run()==0&&!std::memcmp(replay,out,sizeof out),"replay");
 auto clear=[&](){for(const auto&m:out){M zero{};check(!std::memcmp(&m,&zero,sizeof m),"failure clear");}};
 a[3].kinetic_energy_J=3;check(run()!=0,"mono K rejection");clear();std::memcpy(a,saved,sizeof a);
 for(double bad:{-1.,1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){q[4]=bad;check(run()!=0,"uniform rejection");clear();}q[4]=.9;
 for(int slot=0;slot<8;++slot){std::memcpy(a,saved,sizeof a);double*fields[]={&a[0].initial.x_m[0],&a[0].initial.x_m[1],&a[0].initial.x_m[2],&a[0].initial.u_m_s[0],&a[0].initial.u_m_s[1],&a[0].initial.u_m_s[2],&a[0].number_weight,&a[0].kinetic_energy_J};*fields[slot]=std::numeric_limits<double>::quiet_NaN();check(run()!=0,"nonfinite zero-weight input");clear();*fields[slot]=std::numeric_limits<double>::infinity();check(run()!=0,"infinite input");clear();}
 std::memcpy(a,saved,sizeof a);a[1].number_weight=-1;check(run()!=0,"negative N");clear();
 for(auto&m:a){m.number_weight=0;}check(run()!=0,"zero total");clear();
 a[0].number_weight=std::numeric_limits<double>::max();a[1].number_weight=std::numeric_limits<double>::max();check(run()!=0,"total overflow");clear();
 for(auto&m:a){m.number_weight=0;m.kinetic_energy_J=0;}a[0].number_weight=1;a[1].number_weight=std::numeric_limits<double>::denorm_min();check(run()!=0,"positive CDF disappearance");clear();
 M tiny{};tiny.number_weight=std::numeric_limits<double>::denorm_min();double u[2]={0,.5};M t[2]{};
 check(fusion_c_thermal_marker_resample(1,&tiny,2,u,t)!=0,"output underflow");
 check(fusion_c_thermal_marker_resample(1,&tiny,1,u,t)==0&&t[0].number_weight==tiny.number_weight,"representable subnormal");
 tiny.number_weight=3*std::numeric_limits<double>::denorm_min();check(fusion_c_thermal_marker_resample(1,&tiny,2,u,t)!=0,"subnormal closure rejection");
 tiny.number_weight=std::numeric_limits<double>::denorm_min();tiny.kinetic_energy_J=.25;check(fusion_c_thermal_marker_resample(1,&tiny,1,u,t)!=0,"energy underflow");
 tiny.number_weight=std::numeric_limits<double>::max();tiny.kinetic_energy_J=2;check(fusion_c_thermal_marker_resample(1,&tiny,1,u,t)!=0,"energy overflow");
 check(fusion_c_thermal_marker_resample(0,a,5,q,out)!=0,"invalid count");clear();check(fusion_c_thermal_marker_resample(1000001,a,5,q,out)!=0,"count cap");clear();check(fusion_c_thermal_marker_resample(4,nullptr,5,q,out)!=0,"null input");clear();check(fusion_c_thermal_marker_resample(4,a,5,nullptr,out)!=0,"null uniforms");clear();check(fusion_c_thermal_marker_resample(4,a,0,q,out)!=0,"invalid draws");check(fusion_c_thermal_marker_resample(4,a,1000001,q,out)!=0,"draw cap");check(fusion_c_thermal_marker_resample(4,a,5,q,nullptr)!=0,"null output");
 std::memcpy(a,saved,sizeof a);std::vector<double> grid(100);std::vector<M> res(100);for(int i=0;i<100;++i)grid[i]=(i+.5)/100.;check(fusion_c_thermal_marker_resample(4,a,100,grid.data(),res.data())==0,"q grid");int first=0;for(auto&m:res)first+=m.initial.x_m[0]==1;check(first==25,"categorical counts (not iid variance)");
 std::cout<<"thermal marker resample PASS\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
