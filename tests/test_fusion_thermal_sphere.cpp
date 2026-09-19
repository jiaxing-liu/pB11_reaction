#include "fusion_thermal_markers.h"
#include "fusion_prompt_reduce.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
using L=long double;
constexpr double pi=3.1415926535897932384626433832795;
void ck(bool b,const char*s){if(!b)throw std::runtime_error(s);}
double norm(const double*x){return std::hypot(x[0],x[1],x[2]);}
int point(void*,const double*x,int*k){double r=norm(x);*k=r<1?FUSION_BOUNDARY_INSIDE:r>1?FUSION_BOUNDARY_OUTSIDE:FUSION_BOUNDARY_AMBIGUOUS;return 0;}
int field(void*,const double*x,double*b){if(norm(x)>=1)return 999;b[0]=b[1]=b[2]=0;return 0;}
int segment(void*,const double*a,const double*b,double tol,fusion_boundary_segment_value_v1*r){
 *r={};if(norm(a)>=1)return 998;
 // A sphere is convex, so inside endpoints prove the entire chord inside.
 if(norm(b)<1){r->result=FUSION_SEGMENT_CLEAR;return 0;}
 double lo=0,hi=1;
 for(int i=0;i<64&&hi-lo>tol*.25;++i){double m=(lo+hi)/2,x[3];if(m==lo||m==hi)break;for(int k=0;k<3;++k)x[k]=a[k]+m*(b[k]-a[k]);if(norm(x)<1)lo=m;else hi=m;}
 hi=std::min(1.,hi+tol*.125);double x[3];for(int k=0;k<3;++k)x[k]=a[k]+hi*(b[k]-a[k]);
 if(norm(x)<=1)return 0;
 r->result=FUSION_SEGMENT_GEOMETRIC_CANDIDATE;r->fraction_lo=lo;r->fraction_hi=hi;r->tolerance_met=hi-lo<=tol;return 0;
}
std::vector<std::pair<double,double>> gauss(int n){
 std::vector<std::pair<double,double>> r(n);
 for(int i=0;i<(n+1)/2;++i){L z=std::cos(pi*(i+.75L)/(n+.5L)),dp=0;
  for(int it=0;it<60;++it){L p=1,q=0;for(int j=1;j<=n;++j){L old=p;p=((2*j-1)*z*p-(j-1)*q)/j;q=old;}dp=n*(z*p-q)/(z*z-1);L dz=p/dp;z-=dz;if(std::abs(dz)<1e-19L)break;}
  L w=2/((1-z*z)*dp*dp);r[i]={double(-z),double(w)};r[n-i-1]={double(z),double(w)};
 }return r;
}
std::vector<fusion_birth_direction_v1> angles(int n){std::vector<fusion_birth_direction_v1>a;for(auto [mu,w]:gauss(n))for(int j=0;j<2*n;++j){double phi=2*pi*(j+.5)/(2*n),v=std::sqrt(1-mu*mu);a.push_back({{v*std::cos(phi),v*std::sin(phi),mu},w/(4*n)});}return a;}
std::vector<fusion_birth_spatial_node_v1> space(int n){std::vector<fusion_birth_spatial_node_v1>x;auto a=angles(2*n);for(auto [q,w]:gauss(n)){double radius=std::cbrt((1+q)/2);for(auto d:a)x.push_back({{radius*d.direction[0],radius*d.direction[1],radius*d.direction[2]},(4*pi/3)*(w/2)*d.probability_weight});}return x;}
void run(const char*label,int ns,int na){
 auto a=angles(na);auto x=ns?space(ns):std::vector<fusion_birth_spatial_node_v1>{{{0,0,.5},4*pi/3}};
 L volume=0,r2=0,first[3]{},second[3]{};for(auto p:x){volume+=p.volume_weight_m3;for(int k=0;k<3;++k){first[k]+=L(p.volume_weight_m3)*p.x_m[k];second[k]+=L(p.volume_weight_m3)*p.x_m[k]*p.x_m[k];}r2+=L(p.volume_weight_m3)*std::pow(norm(p.x_m),2);}
 ck(std::abs(volume/(4*pi/3)-1)<1e-12,"volume quadrature");if(ns)for(int k=0;k<3;++k){ck(std::abs(first[k]/volume)<1e-12,"volume first moment");ck(std::abs(second[k]/volume-r2/(3*volume))<1e-12,"isotropic volume second moments");}
 constexpr double K=1.602176634e-15,mass=1.67262192369e-27,c=299792458.;
 double u=std::sqrt(K/mass)*std::sqrt(2+K/(mass*c*c)),v=u/std::hypot(1.,u/c),T=1/v;
 std::vector<fusion_thermal_marker_v1> markers(x.size()*a.size());
 ck(fusion_c_thermal_markers(1,1,K,mass,4*pi/3,1e-11,int(x.size()),x.data(),int(a.size()),a.data(),int(markers.size()),markers.data())==0,"tensor construction");
 fusion_prompt_options_v1 o{T,T/2,T*1e-9,T*1e-7,1e-7,1e-8,1e-10,2,4,64};
 std::vector<fusion_prompt_weight_v1> w;w.reserve(markers.size());L exact=0;double max_time=0;int mismatches=0;
 for(auto m:markers){L dot=0,rr=0;for(int k=0;k<3;++k){dot+=L(m.initial.x_m[k])*m.initial.u_m_s[k]/u;rr+=L(m.initial.x_m[k])*m.initial.x_m[k];}
  // Positive root of |x+l*d|^2=1, independent of segment bisection/driver.
  L length=-dot+std::sqrt(dot*dot+1-rr);bool exits=length<=1;
  if(exits)exact+=m.number_weight;
  fusion_prompt_result_v1 r{};int rc=fusion_c_prompt_orbit(mass,1.602176634e-19,&m.initial,&o,field,point,segment,nullptr,&r);
  ck(rc==0,"orbit API error invalidates ensemble");
  if(r.outcome!=FUSION_PROMPT_UNRESOLVED&&((r.outcome==FUSION_PROMPT_EVENT)!=exits))++mismatches;
  if(r.outcome==FUSION_PROMPT_EVENT)max_time=std::max(max_time,std::abs((r.event_time_lo_s+r.event_time_hi_s)/2-double(length)/v)/T);
  w.push_back({m.number_weight,K,rc,r.outcome});
 }
 fusion_prompt_totals_v1 t{};ck(fusion_c_prompt_reduce(int(w.size()),w.data(),&t)==0,"ensemble reduce");
 ck(mismatches==0,"discrete analytical classification mismatch");
 ck(exact>=t.number_fraction_low-1e-12&&exact<=t.number_fraction_high+1e-12,"discrete exact weight outside classification interval");
 ck(std::abs(t.total_number-1)<1e-11&&std::abs(t.total_energy_J/K-1)<1e-11,"ensemble N/E closure");
 ck(std::abs(t.number_fraction_low-t.energy_fraction_low)<1e-12,"monoenergy N/E fraction");
 ck(max_time<1e-7,"straight first-exit time accuracy");
 double target=ns?11./16:5./8;
 std::cout<<label<<','<<ns<<','<<na<<','<<markers.size()<<','<<double(r2/volume)<<','<<double(exact)<<','<<target<<','<<double(exact-target)<<','<<t.number[0]<<','<<max_time<<'\n';
}
int main(){try{std::cout<<std::setprecision(16)<<"family,radial_order,angle_order,markers,mean_r2,discrete_exact,continuous_exact,quadrature_error,unresolved_weight,max_relative_time_error\n";
 for(int a:{4,8,16,32,64})run("angle_only",0,a);
 for(int s:{4,6,8})run("spatial",s,8);
 for(int a:{4,6})run("angular",8,a);
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
