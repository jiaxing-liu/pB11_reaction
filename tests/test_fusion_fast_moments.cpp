#include "fusion_fast_moments.h"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
void check(bool p){if(!p)throw std::runtime_error("fast moment check failed");}
void near(double x,double y){check(y==0?x==0:std::abs(x-y)<=2e-14*std::abs(y));}
int main(){try{
 double edges[]={0,2e-15,6e-15};std::array<double,12>s{},t{};
 // Hand-counted six species, distinct S/T and nonuniform energy cells.
 for(int i=0;i<6;++i){s[2*i]=i+1;t[2*i+1]=2*(i+1);}
 fusion_fast_moments_v1 m{};check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)==0);
 for(int i=0;i<6;++i){near(m.number_m3[i],3*(i+1));near(m.energy_J_m3[i],9e-15*(i+1));}
 // Charges p,D,T,He3,He4,B11 = 1,1,1,2,2,5.
 near(m.charge_number_m3,162);near(m.charge_squared_number_m3,576);near(m.pressure_Pa,126e-15);
 auto original=m;auto total=s;for(int i=0;i<12;++i)total[i]+=t[i];std::array<double,12>zero{};
 check(fusion_c_fast_moments(2,edges,total.data(),zero.data(),&m)==0);
 near(m.pressure_Pa,original.pressure_Pa);near(m.charge_number_m3,original.charge_number_m3);
 auto cleared=[&](){for(double v:m.number_m3)check(v==0);for(double v:m.energy_J_m3)check(v==0);check(m.pressure_Pa==0&&m.charge_number_m3==0&&m.charge_squared_number_m3==0);};
 s[0]=-1;check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)!=0);cleared();s[0]=1;
 edges[1]=edges[2];check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)!=0);cleared();edges[1]=2e-15;
 s[10]=std::numeric_limits<double>::max();check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)!=0);cleared();
 s.fill(0);t.fill(0);s[0]=std::numeric_limits<double>::denorm_min();check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)!=0);cleared();
 s.fill(0);check(fusion_c_fast_moments(2,edges,s.data(),t.data(),&m)==0);cleared();
 check(fusion_c_fast_moments(0,edges,s.data(),t.data(),&m)!=0);cleared();
 std::cout<<"PASS canonical charge moments, nonuniform S/T energy, NR pressure and rejection clearing\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
