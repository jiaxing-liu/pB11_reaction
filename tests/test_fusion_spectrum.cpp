#include "fusion_nuclear_coulomb.h"
#include "fusion_alpha_amplitudes.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "nuclear_coulomb_reference.inc"
namespace {
constexpr double mev=1.602176634e-13,m=3727.3794118;
void require(bool x,const char *label){if(!x)throw std::runtime_error(label);}
bool near(double a,double b,double tol){return std::abs(a-b)<=tol*std::max(std::abs(a),std::abs(b));}
double norm(const fusion_alpha_amplitudes_v1 &a,bool sym){double s=0;for(int i=0;i<5;++i){double re=sym?a.sym_real[i]:a.unsym_real[i],im=sym?a.sym_imag[i]:a.unsym_imag[i];s+=re*re+im*im;}return s;}
void reference(){
 fusion_nuclear_coulomb_v1 r;
 for(const auto &v:nuclear_references){
  require(fusion_c_nuclear_coulomb(v.channel,v.E*mev,&r)==0,"independent reference evaluation");
  require(std::abs(r.log_penetrability-v.logP)<2e-7,"off-node log penetrability reference");
  require(std::abs(r.shift-v.S)<2e-8*(1+std::abs(v.S)),"off-node shift reference");
  require(std::hypot(r.phase_real-v.re,r.phase_imag-v.im)<2e-7,"off-node Coulomb phase reference");
 }

 require(fusion_c_nuclear_coulomb(4,3.129*mev,&r)==0,"boundary coefficient status");
 require(std::abs(r.log_penetrability-std::log(.9110014643199444962))<2e-7,"independent high-precision penetrability");
 require(std::abs(r.shift+.923134650192041157)<4e-8,"independent high-precision shift");
 for(int ch=0;ch<5;++ch)for(double E:{.001,.01,.1,1.,3.129,6.,12.}){
  require(fusion_c_nuclear_coulomb(ch,E*mev,&r)==0,"Coulomb domain including endpoints");
  require(std::abs(std::hypot(r.phase_real,r.phase_imag)-1)<3e-7,"Coulomb phase norm");}
 require(fusion_c_nuclear_coulomb(4,.0009*mev,&r)==PB11_STATUS_OUT_OF_RANGE && r.rho==0,"no unvalidated extrapolation");
}
void shapes(){
 int checked=0;
 for(int l=1;l<=3;++l)for(double q:{1.,3.129,6.})for(double cosine:{-.6,.1,.7}){
  const double A=9.3;
  fusion_alpha_amplitudes_v1 a,b;
  require(fusion_c_alpha_amplitudes(l,A*mev,q*mev,cosine,&a)==0,"amplitude status");
  require(norm(a,true)>0 && std::isfinite(norm(a,true)),"finite positive coherent weight");
  require(fusion_c_alpha_amplitudes(l,A*mev,q*mev,-cosine,&b)==0,"secondary reflection status");
  require(near(norm(a,true),norm(b,true),3e-6),"identical-secondary exchange");
  const double p0=std::sqrt(4*m*(A-q)/3),star=std::sqrt(m*q),px=star*std::sqrt(1-cosine*cosine);
  const double p[3][3]={{0,0,p0},{px,0,-p0/2-star*cosine},{-px,0,-p0/2+star*cosine}};
  // Relabel another alpha as primary and reconstruct Jacobi variables;
  // the summed spin norm must be invariant under the accompanying rotation.
  for(int i=1;i<3;++i){int j=(i+1)%3,k=(i+2)%3;double qnew=0,pnorm=0,dot=0;
   for(int d=0;d<3;++d){double relative=(p[j][d]-p[k][d])/2;qnew+=relative*relative/m;pnorm+=p[i][d]*p[i][d];dot+=p[i][d]*relative;}
   double cnew=-dot/std::sqrt(pnorm*m*qnew);
   require(fusion_c_alpha_amplitudes(l,A*mev,qnew*mev,cnew,&b)==0,"cyclic event status");
   require(near(norm(a,true),norm(b,true),3e-6),"coherent cyclic and rotational invariance");}
  ++checked;
 }
 // Independently derived J=2,l=2,Jb=2 angular algebra for ONE permutation.
 // This tests the explicit CG amplitude; Laursen Eq2's separate printed
 // angle fit remains a documented discrepancy, not silently substituted.
 fusion_alpha_amplitudes_v1 base,x;
 require(fusion_c_alpha_amplitudes(2,9.3*mev,3.129*mev,0,&base)==0,"unsym base");
 for(double c:{.2,.5,.7071067811865475,.9,1.}){
  require(fusion_c_alpha_amplitudes(2,9.3*mev,3.129*mev,c,&x)==0,"unsym angle");
  double expected=1-2.25*c*c*(1-c*c);
  require(near(norm(x,false)/norm(base,false),expected,2e-12),"independent unsym quadrupole angular polynomial");}
 require(fusion_c_alpha_amplitudes(2,9.3*mev,0,0,&x)==PB11_STATUS_OUT_OF_RANGE && x.phase_space_J==0,"invalid event clears output");
 std::cout<<"PASS: "<<checked<<" coherent events, cyclic/secondary/rotation checks and independent angular polynomial\n";
}
}
int main(){try{reference();shapes();return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
