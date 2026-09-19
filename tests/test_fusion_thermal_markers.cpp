#include "fusion_thermal_markers.h"
#include "fusion_prompt_reduce.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
void check(bool ok,const char*s){if(!ok)throw std::runtime_error(s);}
int main(){try{
 fusion_birth_spatial_node_v1 x[2]={{{1,2,3},1},{{4,5,6},3}};
 fusion_birth_direction_v1 a[6]{};
 for(int k=0;k<3;++k){a[2*k].direction[k]=1;a[2*k+1].direction[k]=-1;}
 for(auto&v:a)v.probability_weight=1./6;
 fusion_thermal_marker_v1 m[12]{},saved[12]{};
 constexpr double mass=1.67262192369e-27,K=1.602176634e-15,c=299792458.;
 auto run=[&](double N=120,double energy=1.602176634e-15){return fusion_c_thermal_markers(1,N,energy,mass,4,1e-12,2,x,6,a,12,m);};
 check(run()==0,"valid tensor");std::memcpy(saved,m,sizeof(m));
 fusion_prompt_weight_v1 items[12]{};
 long double n=0,e=0,mean[3]{};
 for(int i=0;i<12;++i){
  check(std::abs(m[i].number_weight-(i<6?5:15))<1e-13,"unequal spatial volumes");
  double u2=0;for(int k=0;k<3;++k){check(m[i].initial.x_m[k]==x[i/6].x_m[k],"spatial order");u2+=m[i].initial.u_m_s[k]*m[i].initial.u_m_s[k];mean[k]+=m[i].number_weight*m[i].initial.u_m_s[k];}
  double reconstructed=mass*u2/(std::sqrt(1+u2/(c*c))+1);
  check(std::abs(reconstructed/K-1)<1e-14,"proper velocity reconstructs K");
  n+=m[i].number_weight;e+=(long double)m[i].number_weight*K;
  items[i]={m[i].number_weight,K,0,i%3};
 }
 check(std::abs(n/120-1)<1e-14&&std::abs(e/(120*K)-1)<1e-14,"N/E weights");
 for(auto v:mean)check(v==0,"antipodal momentum");
 fusion_prompt_totals_v1 totals{};check(fusion_c_prompt_reduce(12,items,&totals)==0,"reduce generated weights");
 check(std::abs(totals.total_number/120-1)<1e-14&&std::abs(totals.number_fraction_low-1./3)<1e-14,"classified reconstruction");
 check(run()==0&&std::memcmp(saved,m,sizeof(m))==0,"deterministic replay");
 auto cleared=[&]{fusion_thermal_marker_v1 zero[12]{};return std::memcmp(zero,m,sizeof(m))==0;};
 x[1].volume_weight_m3=2;check(run()!=0&&cleared(),"bad shell volume atomic");x[1].volume_weight_m3=3;
 a[0].direction[0]=2;check(run()!=0&&cleared(),"nonunit direction");a[0].direction[0]=1;
 for(auto&v:a){v.direction[0]=1;v.direction[1]=v.direction[2]=0;}
 check(run()!=0&&cleared(),"anisotropic rule");
 for(int k=0;k<3;++k){a[2*k].direction[0]=a[2*k+1].direction[0]=0;a[2*k].direction[k]=1;a[2*k+1].direction[k]=-1;}
 check(run(std::numeric_limits<double>::denorm_min())!=0&&cleared(),"positive split underflow rejected");
 check(run(0)==0,"zero amount");for(auto v:m)check(v.number_weight==0,"zero amount weights");
 check(run(120,0)==0,"zero K");for(auto v:m)for(double q:v.initial.u_m_s)check(q==0,"zero K proper velocity");
 check(fusion_c_thermal_markers(1,120,K,mass,4,1e-12,1000000,x,1000000,a,12,m)!=0&&cleared(),"product overflow before access");
 // A rotated cube rule and a mixed rule exercise non-axis directions and
 // nonuniform angular weights independently of the six-axis reference rule.
 std::vector<fusion_birth_direction_v1> cube;
 for(int sx:{-1,1})for(int sy:{-1,1})for(int sz:{-1,1}){
  double q=1/std::sqrt(3.),angle=.371;
  cube.push_back({{q*(sx*std::cos(angle)-sy*std::sin(angle)),
                   q*(sx*std::sin(angle)+sy*std::cos(angle)),q*sz},1./8});
 }
 std::vector<fusion_birth_direction_v1> mixed(a,a+6);
 for(auto&v:mixed)v.probability_weight*=.5;
 for(auto v:cube){v.probability_weight*=.5;mixed.push_back(v);}
 for(const auto&rule:{cube,mixed})for(double energy:{1e-30,K,1e-10,1e-7}){
  std::vector<fusion_thermal_marker_v1> out(2*rule.size());
  check(fusion_c_thermal_markers(1,120,energy,mass,4,1e-12,2,x,int(rule.size()),rule.data(),int(out.size()),out.data())==0,"rotated rule scales");
  long double sum=0;
  for(auto&r:out){
   long double u2=0;for(double u:r.initial.u_m_s)u2+=(long double)u*u;
   long double reconstructed=mass*u2/(std::sqrt(1+u2/((long double)c*c))+1);
   check(std::abs(reconstructed/energy-1)<2e-14L,"nonrelativistic to ultrarelativistic energy reconstruction");
   sum+=r.number_weight;
  }
  check(std::abs(sum/120-1)<1e-14L,"mixed angular weight conservation");
 }
 // Neither caller measure nor directions may be silently repaired.
 a[0].probability_weight*=1.1;check(run()!=0&&cleared(),"bad angular measure");a[0].probability_weight/=1.1;
 x[0].x_m[0]=std::numeric_limits<double>::infinity();check(run()!=0&&cleared(),"nonfinite position");x[0].x_m[0]=1;
 check(run(-1)!=0&&cleared(),"negative amount");
 check(run(120,-1)!=0&&cleared(),"negative kinetic energy");
 check(run(1e308,1e308)!=0&&cleared(),"unrepresentable total energy");
 std::cout<<"PASS thermal tensor weights, proper velocity, reduction, replay and atomic rejection\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
