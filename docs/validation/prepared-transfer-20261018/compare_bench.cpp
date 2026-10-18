#include "fusion_transfer_prepared_internal.h"
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <random>
#include <vector>
#include <cstdint>
#include <cstdlib>
extern "C" int reference_coulomb_transfer_rate(double,double,double,const fusion_maxwellian_bath_v1*,double*);
uint64_t bits(double x){uint64_t u;std::memcpy(&u,&x,8);return u;}
size_t count=0,prepared_count=0,zero=0,subnormal=0,overflow=0;
void check(double energy,double mass,double charge,const fusion_maxwellian_bath_v1* b,bool nullout=false){
 double a=-17,z=-17;int sa=reference_coulomb_transfer_rate(energy,mass,charge,b,nullout?nullptr:&a),sz=fusion_c_coulomb_transfer_rate(energy,mass,charge,b,nullout?nullptr:&z);++count;
 if(sa!=sz||bits(a)!=bits(z)){std::printf("FAIL scalar #%zu status %d/%d bits %llx/%llx\n",count,sa,sz,(unsigned long long)bits(a),(unsigned long long)bits(z));std::exit(1);}
 if(sa==PB11_STATUS_NUMERICAL_FAILURE)++overflow;
 if(sa==0&&!nullout){if(a==0)++zero;if(std::fpclassify(a)==FP_SUBNORMAL)++subnormal;}
 bool valid=b&&std::isfinite(energy)&&energy>=0&&std::isfinite(mass)&&mass>0&&std::isfinite(charge)&&std::isfinite(b->density_m3)&&b->density_m3>=0&&std::isfinite(b->mass_kg)&&b->mass_kg>0&&std::isfinite(b->mean_charge_squared)&&b->mean_charge_squared>=0&&std::isfinite(b->kT_J)&&b->kT_J>0&&std::isfinite(b->coulomb_log)&&b->coulomb_log>0;
 if(valid&&!nullout){fusion_detail::prepared_transfer p;fusion_detail::prepare_transfer_validated(mass,charge,*b,p);double q=-17;int sq=fusion_detail::evaluate_transfer_prepared(energy,p,&q);++prepared_count;if(sa!=sq||bits(a)!=bits(q)){std::printf("FAIL prepared #%zu status %d/%d bits %llx/%llx\n",count,sa,sq,(unsigned long long)bits(a),(unsigned long long)bits(q));std::exit(1);}}
}
volatile double sink=0;
int main(){
 constexpr double kev=1.602176634e-16,amu=1.66053906892e-27;
 fusion_maxwellian_bath_v1 base{5e19,2*amu,1,10*kev,15};
 const double ma=4*amu,max=std::numeric_limits<double>::max(),dn=std::numeric_limits<double>::denorm_min(),mn=std::numeric_limits<double>::min();
 std::vector<double> ex={-INFINITY,-max,-1.,-dn,-0.,0.,dn,mn,1e-300,1e-100,1e-27,1.,1e100,1e300,max,INFINITY,NAN};
 for(int field=0;field<8;++field)for(double v:ex){double e=kev,m=ma,c=2;auto b=base;if(field==0)e=v;if(field==1)m=v;if(field==2)c=v;if(field==3)b.density_m3=v;if(field==4)b.mass_kg=v;if(field==5)b.mean_charge_squared=v;if(field==6)b.kT_J=v;if(field==7)b.coulomb_log=v;check(e,m,c,&b);check(e,m,c,&b,true);}
 check(0,ma,2,nullptr);check(NAN,-1,INFINITY,nullptr,true);
 // Extreme Cartesian products across dynamic scales; finite positive domain.
 std::vector<double> scales={dn,mn,1e-300,1e-100,1e-27,1.,1e100,1e300,max};
 for(double t:scales)for(double mb:scales)for(double n:scales)for(double e:scales){auto b=base;b.kT_J=t;b.mass_kg=mb;b.density_m3=n;check(e,ma,2,&b);}
 // Dedicated y=.01 neighbors plus tails near double normal/subnormal/zero.
 for(double y2:{0.,std::nextafter(.0001,0.),.0001,std::nextafter(.0001,1.),600.,700.,710.,720.,730.,740.,750.,760.,800.}){check(y2*base.kT_J*(ma/base.mass_kg),ma,2,&base);}
 std::mt19937_64 gen(885);std::uniform_real_distribution<double> logscale(-323,308),normal(-1,1);
 for(int i=0;i<40000;++i){auto pos=[&](){return std::pow(10.,logscale(gen));};auto b=base;b.density_m3=pos();b.mass_kg=pos();b.mean_charge_squared=pos();b.kT_J=pos();b.coulomb_log=pos();check(pos(),pos(),normal(gen)*pos(),&b);}
 // Thousands of invalid variants check status precedence and clearing.
 for(int i=0;i<6000;++i){auto b=base;double en=kev,ms=ma,ch=2;switch(i%10){case 0:en=NAN;break;case 1:ms=0;break;case 2:ch=INFINITY;break;case 3:b.density_m3=-1;break;case 4:b.mass_kg=0;break;case 5:b.mean_charge_squared=-1;break;case 6:b.kT_J=-1;break;case 7:b.coulomb_log=0;break;case 8:en=-1;b.density_m3=0;break;case 9:b.kT_J=NAN;ch=0;break;}check(en,ms,ch,&b);}
 std::ifstream input("edges.txt");std::vector<double> edges;double x;while(input>>x)edges.push_back(x);if(edges.size()!=509)return 2;
 std::vector<double> centers;for(int j=0;j<508;++j)centers.push_back(double(((long double)edges[j]+edges[j+1])/2));
 std::array<double,6> masses={amu,2*amu,3*amu,3*amu,4*amu,11*amu},charges={1,1,1,2,2,5};
 std::array<fusion_maxwellian_bath_v1,7> baths;for(int j=0;j<6;++j)baths[j]={5e19,masses[j],charges[j]*charges[j],.1*kev,10};baths[6]={7e16,12*amu,36,.1*kev,10};
 for(int i=0;i<6;++i)for(auto& b:baths)for(double e:centers)check(e,masses[i],charges[i],&b);
 std::printf("comparison count=%zu prepared=%zu mismatch=0 zero=%zu subnormal=%zu numerical_failure=%zu\n",count,prepared_count,zero,subnormal,overflow);
 // Include preparation once per projectile/bath/repetition, matching coupled.
 auto bench=[&](bool prepared){auto start=std::chrono::steady_clock::now();double total=0;for(int rep=0;rep<100;++rep)for(int i=0;i<6;++i)for(auto& b:baths){fusion_detail::prepared_transfer p;if(prepared)fusion_detail::prepare_transfer_validated(masses[i],charges[i],b,p);for(double e:centers){double value;int st=prepared?fusion_detail::evaluate_transfer_prepared(e,p,&value):reference_coulomb_transfer_rate(e,masses[i],charges[i],&b,&value);if(st)std::exit(3);total+=value;}}sink=total;return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
 bench(false);bench(true);for(int pair=0;pair<6;++pair){double scalar,prep;if(pair%2){prep=bench(true);scalar=bench(false);}else{scalar=bench(false);prep=bench(true);}std::printf("pair=%d scalar=%.9f prepared=%.9f speedup=%.6f checksum=%.17g evaluations=%d\n",pair,scalar,prep,scalar/prep,(double)sink,100*6*7*508);std::fflush(stdout);}
}
