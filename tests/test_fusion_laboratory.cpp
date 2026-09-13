#include "fusion_laboratory.h"
#include "fusion_products.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace {
using R=long double;
constexpr R c=299792458.L,c2=c*c,kev=1.602176634e-16L;
void check(bool b,const char *s){if(!b) throw std::runtime_error(s);}
void near(R a,R b,R scale,R tol=1e-11L){check(std::abs(a-b)<=tol*scale,"scaled comparison");}
void shell(const fusion_particle_four_vector_v1 &p){
 R pc2=0;for(double v:p.momentum_kg_m_s) pc2+=static_cast<R>(v)*v*c2;
 R t=p.kinetic_energy_J,r=static_cast<R>(p.mass_kg)*c2;
 near(t*(t+2*r),pc2,std::max(pc2,t*(t+2*r)));
}
void two_body(){
 const double d[3]={.6,0,.8};
 for(double ma:{1.67e-27,6.64e-27}) for(double mb:{1.67e-27,5.00e-27})
 for(double energy:{0.,double(1*kev),double(18000*kev)}) {
  fusion_particle_four_vector_v1 p[2];
  check(fusion_c_two_body_cm(ma,mb,energy,d,p)==0,"two-body status");
  shell(p[0]);shell(p[1]);
  near(static_cast<R>(p[0].kinetic_energy_J)+p[1].kinetic_energy_J,energy,energy);
  for(int j=0;j<3;++j) near(static_cast<R>(p[0].momentum_kg_m_s[j])+p[1].momentum_kg_m_s[j],0,std::sqrt(2*ma*energy));
  // Independent total-energy ratio from the two-body invariant masses.
  R ra=static_cast<R>(ma)*c2,rb=static_cast<R>(mb)*c2,W=ra+rb+energy;
  R expected_energy_a=(W*W+ra*ra-rb*rb)/(2*W);
  near(ra+p[0].kinetic_energy_J,expected_energy_a,W);
 }
 fusion_particle_four_vector_v1 p[2];const double bad[3]={1,1,0};
 check(fusion_c_two_body_cm(1,1,1,bad,p)==PB11_STATUS_OUT_OF_RANGE && p[0].mass_kg==0,"reject nonunit direction");
 check(fusion_c_two_body_cm(0,1,1,d,p)==PB11_STATUS_OUT_OF_RANGE,"reject mass");
}
void boosts(){
 const double d[3]={.6,0,.8};fusion_particle_four_vector_v1 start[2];
 check(fusion_c_two_body_cm(6.644657345e-27,1.67492750056e-27,double(17589*kev),d,start)==0,"DT-like pair");
 int cases=0;
 for(double speed:{0.,1e3,1e6,.1*double(c),.8*double(c)})
 for(int axis=0;axis<3;++axis){
  double v[3]={};v[axis]=speed;fusion_particle_four_vector_v1 lab[2],back[2];fusion_boost_ledger_v1 l;
  check(fusion_c_boost_particles(2,v,start,lab,&l)==0,"boost status");
  R before_e=0,after_e=0,after_p[3]={};
  for(int i=0;i<2;++i){shell(lab[i]);before_e+=static_cast<R>(start[i].mass_kg)*c2+start[i].kinetic_energy_J;
   after_e+=static_cast<R>(lab[i].mass_kg)*c2+lab[i].kinetic_energy_J;
   for(int j=0;j<3;++j)after_p[j]+=lab[i].momentum_kg_m_s[j];}
  R invariant=after_e*after_e;for(R x:after_p)invariant-=x*x*c2;
  near(invariant,before_e*before_e,before_e*before_e);
  // Since initial total momentum is zero, total lab energy is gamma*Ecm.
  R gamma=1/std::sqrt(1-static_cast<R>(speed)*speed/c2);
  near(after_e,gamma*before_e,after_e);
  v[axis]=-speed;
  check(fusion_c_boost_particles(2,v,lab,back,&l)==0,"inverse boost");
  for(int i=0;i<2;++i){near(back[i].kinetic_energy_J,start[i].kinetic_energy_J,start[i].kinetic_energy_J);
   for(int j=0;j<3;++j)near(back[i].momentum_kg_m_s[j],start[i].momentum_kg_m_s[j],std::sqrt(2*start[i].mass_kg*start[i].kinetic_energy_J));}
  ++cases;
 }
 // Non-collinear velocity and arbitrary three-alpha plane orientation.
 fusion_three_body_cm_v1 event;check(fusion_c_three_equal_sequential_cm(6.644657345e-27,double(8680*kev),double(2700*kev),.37,&event)==0,"three-body seed");
 fusion_particle_four_vector_v1 alpha[3]{},lab[3],back[3];
 for(int i=0;i<3;++i){alpha[i].mass_kg=6.644657345e-27;alpha[i].kinetic_energy_J=event.kinetic_energy_J[i];
  alpha[i].momentum_kg_m_s[0]=.6*event.momentum_x_kg_m_s[i];alpha[i].momentum_kg_m_s[1]=.8*event.momentum_x_kg_m_s[i];alpha[i].momentum_kg_m_s[2]=event.momentum_z_kg_m_s[i];}
 double v[3]={1e6,2e6,-3e6};fusion_boost_ledger_v1 l;
 check(fusion_c_boost_particles(3,v,alpha,lab,&l)==0,"rotated three-body laboratory boost");
 for(double &x:v)x=-x;
 check(fusion_c_boost_particles(3,v,lab,back,&l)==0,"three-body inverse");
 for(int i=0;i<3;++i)near(back[i].kinetic_energy_J,alpha[i].kinetic_energy_J,alpha[i].kinetic_energy_J);
 // At-rest input becomes the expected classical limit at small speed.
 fusion_particle_four_vector_v1 rest{3e-27,0,{0,0,0}},moving;
 double small[3]={1e3,0,0};
 check(fusion_c_boost_particles(1,small,&rest,&moving,&l)==0,"rest acceleration");
 near(moving.kinetic_energy_J,.5L*rest.mass_kg*1e6L,.5L*rest.mass_kg*1e6L,1e-10L);
 small[0]=-small[0];
 check(fusion_c_boost_particles(1,small,&moving,&rest,&l)==0,"stopping deboost");
 check(rest.kinetic_energy_J<1e-25*moving.kinetic_energy_J,"stopped product cancellation");
 small[0]=double(c);
 check(fusion_c_boost_particles(1,small,&moving,&rest,&l)==PB11_STATUS_OUT_OF_RANGE && rest.mass_kg==0,"reject luminal boost");
 small[0]=0;moving.kinetic_energy_J*=2;
 check(fusion_c_boost_particles(1,small,&moving,&rest,&l)==PB11_STATUS_OUT_OF_RANGE,"reject off-shell input");
 std::cout<<"PASS: "<<cases<<" pair boosts, inverse and total invariants; rotated three-alpha boost; stopped limit\n";
}
}
int main(){try{two_body();boosts();return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
