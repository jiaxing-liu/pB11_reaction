#include "fusion_reaction_event.h"
#include "fusion_products.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace {
using R=long double;
constexpr R c=299792458.L,c2=c*c;
bool put(R x,double&y){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;y=double(x);return x==0||y!=0;}
R energy(double mass,const double*p){R pc2=0;for(int j=0;j<3;++j)pc2+=R(p[j])*p[j]*c2;R rest=R(mass)*c2;return pc2/(std::sqrt(rest*rest+pc2)+rest);}
void cross(const R*a,const R*b,R*out){for(int j=0;j<3;++j)out[j]=a[(j+1)%3]*b[(j+2)%3]-a[(j+2)%3]*b[(j+1)%3];}
}
extern "C" int fusion_c_reaction_parent(int ch,int convention,const double*pa,const double*pb,fusion_reaction_parent_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};if(!pa||!pb||(convention!=0&&convention!=1))return PB11_STATUS_INVALID_ARGUMENT;
 fusion_nuclear_channel_v1 reaction{};int status=fusion_c_nuclear_channel(ch,&reaction);if(status)return status;
 fusion_nuclear_mass_v1 mass[2]{},product{};
 for(int i=0;i<2;++i){status=fusion_c_nuclear_mass(reaction.reactant_ids[i],&mass[i]);if(status)return status;}
 R restout=0;for(int i=0;i<reaction.product_count;++i){status=fusion_c_nuclear_mass(reaction.product_ids[i],&product);if(status)return status;restout+=R(product.mass_kg)*c2;}
 fusion_reaction_parent_v1 result{};R chosen=0,correction=0,P[3]{},p2=0,relv2=0;
 const double*input[2]={pa,pb};
 for(int i=0;i<2;++i){
  R momentum2=0;for(int j=0;j<3;++j){double x=input[i][j];if(!std::isfinite(x))return PB11_STATUS_INVALID_ARGUMENT;momentum2+=R(x)*x;P[j]+=x;}
  R nr=momentum2/(2*R(mass[i].mass_kg)),shell=energy(mass[i].mass_kg,input[i]);
  if(!put(nr,result.classical_kinetic_J[i])||!put(shell,result.on_shell_kinetic_J[i]))return PB11_STATUS_NUMERICAL_FAILURE;
  R selected=convention? shell:nr;
  if(!put(selected,result.reactant_kinetic_J[i]))return PB11_STATUS_NUMERICAL_FAILURE;
  chosen+=selected;correction+=shell*shell/(2*R(mass[i].mass_kg)*c2);
 }
 for(int j=0;j<3;++j){p2+=P[j]*P[j];R v=R(pa[j])/mass[0].mass_kg-R(pb[j])/mass[1].mass_kg;relv2+=v*v;}
 R mu=R(mass[0].mass_kg)/(1+R(mass[0].mass_kg)/mass[1].mass_kg);
 const R L=chosen+reaction.q_J,total=restout+L;
 // Stable invariant excess: subtract only kinetic-scale quantities, never
 // two nearly equal nuclear rest energies to obtain the MeV-scale release.
 R delta=L*(2*restout+L)-p2*c2;
 if(!std::isfinite(delta)||delta<0)return PB11_STATUS_NUMERICAL_FAILURE;
 R invariant=std::sqrt(restout*restout+delta),available=delta/(invariant+restout);
 if(!put(mu*relv2/2,result.relative_classical_energy_J)||!put(available,result.available_cm_energy_J)||
    !put(L,result.expected_product_lab_kinetic_J)||!put(correction,result.classical_minus_on_shell_J))return PB11_STATUS_NUMERICAL_FAILURE;
 R v2=0;for(int j=0;j<3;++j){if(!put(P[j],result.momentum_sum_kg_m_s[j])||!put(c2*P[j]/total,result.boost_velocity_m_s[j]))return PB11_STATUS_NUMERICAL_FAILURE;v2+=R(result.boost_velocity_m_s[j])*result.boost_velocity_m_s[j];}
 if(v2>=c2)return PB11_STATUS_NUMERICAL_FAILURE;
 result.product_count=reaction.product_count;result.convention=convention;*out=result;return PB11_STATUS_OK;
}
extern "C" int fusion_c_reaction_lab_event(int ch,int convention,const double*pa,const double*pb,
 const double*direction,double q,double cosine,double azimuth,fusion_particle_four_vector_v1*out,
 fusion_reaction_parent_v1*parent,fusion_boost_ledger_v1*ledger){
 if(out)std::fill(out,out+3,fusion_particle_four_vector_v1{});
 if(parent)*parent={};
 if(ledger)*ledger={};
 if(!out||!parent||!ledger)return PB11_STATUS_NULL_OUTPUT;
 if(!direction||!std::isfinite(q)||!std::isfinite(cosine)||!std::isfinite(azimuth))return PB11_STATUS_INVALID_ARGUMENT;
 fusion_reaction_parent_v1 p{};int status=fusion_c_reaction_parent(ch,convention,pa,pb,&p);if(status)return status;
 R n[3],norm2=0;for(int j=0;j<3;++j){if(!std::isfinite(direction[j]))return PB11_STATUS_INVALID_ARGUMENT;n[j]=direction[j];norm2+=n[j]*n[j];}
 if(std::abs(norm2-1)>2e-12L)return PB11_STATUS_OUT_OF_RANGE;
 R norm=std::sqrt(norm2);for(auto&x:n)x/=norm;
 fusion_nuclear_channel_v1 reaction{};status=fusion_c_nuclear_channel(ch,&reaction);if(status)return status;
 fusion_nuclear_mass_v1 masses[3]{};for(int i=0;i<p.product_count;++i){status=fusion_c_nuclear_mass(reaction.product_ids[i],&masses[i]);if(status)return status;}
 fusion_particle_four_vector_v1 cm[3]{},lab[3]{};
 if(p.product_count==2){
  if(q!=0||cosine!=0||azimuth!=0)return PB11_STATUS_OUT_OF_RANGE;
  status=fusion_c_two_body_cm(masses[0].mass_kg,masses[1].mass_kg,p.available_cm_energy_J,direction,cm);if(status)return status;
 }else{
  if(q<0||q>p.available_cm_energy_J||std::abs(cosine)>1)return PB11_STATUS_OUT_OF_RANGE;
  fusion_three_body_cm_v1 seq{};status=fusion_c_three_equal_sequential_cm(masses[0].mass_kg,p.available_cm_energy_J,q,cosine,&seq);if(status)return status;
  R reference[3]={0,0,1},t[3],b[3];if(std::abs(n[2])>=.9L){reference[0]=1;reference[2]=0;}
  cross(reference,n,t);R tnorm=0;for(auto x:t)tnorm+=x*x;tnorm=std::sqrt(tnorm);for(auto&x:t)x/=tnorm;cross(n,t,b);
  R ct=std::cos(R(azimuth)),st=std::sin(R(azimuth));
  for(int i=0;i<3;++i){cm[i].mass_kg=masses[i].mass_kg;for(int j=0;j<3;++j){R x=R(seq.momentum_z_kg_m_s[i])*n[j]+R(seq.momentum_x_kg_m_s[i])*(ct*t[j]+st*b[j]);if(!put(x,cm[i].momentum_kg_m_s[j]))return PB11_STATUS_NUMERICAL_FAILURE;}
   if(!put(energy(cm[i].mass_kg,cm[i].momentum_kg_m_s),cm[i].kinetic_energy_J))return PB11_STATUS_NUMERICAL_FAILURE;
  }
 }
 fusion_boost_ledger_v1 boost{};status=fusion_c_boost_particles(p.product_count,p.boost_velocity_m_s,cm,lab,&boost);if(status)return status;
 R K=0,P[3]{},pscale=0;for(int i=0;i<p.product_count;++i){K+=lab[i].kinetic_energy_J;for(int j=0;j<3;++j){P[j]+=lab[i].momentum_kg_m_s[j];pscale+=std::abs(R(lab[i].momentum_kg_m_s[j]));}}
 R pr2=0;for(int j=0;j<3;++j){R residual=P[j]-R(pa[j])-pb[j];pr2+=residual*residual;pscale+=std::abs(R(pa[j]))+std::abs(R(pb[j]));}
 R er=K-p.expected_product_lab_kinetic_J,pr=std::sqrt(pr2);fusion_boost_ledger_v1 result{};
 result.expected_kinetic_energy_J=p.expected_product_lab_kinetic_J;
 if(!put(K,result.output_kinetic_energy_J)||!put(er,result.energy_residual_J)||!put(pr,result.momentum_residual_kg_m_s)||
    std::abs(er)>1e-10L*(K+p.expected_product_lab_kinetic_J)||pr>1e-12L*pscale)return PB11_STATUS_NUMERICAL_FAILURE;
 std::copy(lab,lab+3,out);*parent=p;*ledger=result;return PB11_STATUS_OK;
}
