#include "fusion_nuclear_data.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
namespace {
using R=long double;
constexpr R c=299792458.L,c2=c*c,e=1.602176634e-19L;
constexpr R u=1.66053906892e-27L,su=.00000000052e-27L;
// CODATA2022 relative nuclear masses (u), and their quoted standard errors.
// B slot is the AME2020 NEUTRAL atomic mass, converted explicitly below.
constexpr R ratio[8]={1.0072764665789L,2.013553212544L,3.01550071597L,
 3.014932246932L,4.001506179129L,11.009305166L,1.00866491606L,5.485799090441e-4L};
constexpr R sr[8]={.0000000000083L,.000000000015L,.00000000010L,
 .000000000074L,.000000000062L,.000000013L,.00000000040L,.000000000097e-4L};
constexpr int Z[8]={1,1,1,2,2,5,0,-1},A[8]={1,2,3,3,4,11,1,0};
constexpr R binding_eV=8.298019L+25.15483L+37.93059L+259.374379L+340.2260225L;
constexpr R binding_scale_eV=.000003L+.00005L+.00007L+.000009L+.0000006L;
constexpr int reactants[5][2]={{0,5},{1,1},{1,1},{1,2},{1,3}};
constexpr int products[5][3]={{4,4,4},{2,0,-1},{3,6,-1},{4,6,-1},{4,0,-1}};
R u_part(int id){return ratio[id]-(id==5?5*ratio[7]:0);}
R relative_scale(int id){return sr[id]+(id==5?5*sr[7]:0);}
R nominal_mass(int id){return u_part(id)*u+(id==5?binding_eV*e/c2:0);}
double canonical_mass(int id){return static_cast<double>(nominal_mass(id));}
}
extern "C" int fusion_c_nuclear_mass(int id,fusion_nuclear_mass_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};if(id<0||id>7)return PB11_STATUS_OUT_OF_RANGE;
 out->mass_kg=canonical_mass(id);
 out->rest_energy_J=static_cast<double>(R(out->mass_kg)*c2);
 out->reference_mass_u=static_cast<double>(nominal_mass(id)/u);
 out->known_uncertainty_scale_kg=static_cast<double>(relative_scale(id)*u+
     std::abs(u_part(id))*su+(id==5?binding_scale_eV*e/c2:0));
 out->nuclear_charge=Z[id];out->mass_number=A[id];return PB11_STATUS_OK;
}
extern "C" int fusion_c_nuclear_channel(int ch,fusion_nuclear_channel_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};if(ch<0||ch>=5)return PB11_STATUS_OUT_OF_RANGE;
 fusion_nuclear_channel_v1 r{};
 std::array<int,8>weight{};R mass_difference=0;
 for(int j=0;j<2;++j){int id=reactants[ch][j];r.reactant_ids[j]=id;++weight[id];}
 r.product_count=ch==0?3:2;
 for(int j=0;j<3;++j){int id=products[ch][j];r.product_ids[j]=id;if(id>=0)--weight[id];}
 R relative_difference=0,scale=0;
 for(int id=0;id<8;++id){
  mass_difference+=R(weight[id])*canonical_mass(id);
  relative_difference+=R(weight[id])*u_part(id);
  scale+=std::abs(weight[id])*relative_scale(id);
 }
 r.q_J=static_cast<double>(mass_difference*c2);
 // Shared u uncertainty multiplies the MASS DIFFERENCE, not each full mass.
 // Unknown correlations among quoted ratios are bounded by a linear sum.
 r.known_uncertainty_scale_J=static_cast<double>((u*scale+std::abs(relative_difference)*su)*c2+
     std::abs(weight[5])*binding_scale_eV*e);
 *out=r;return PB11_STATUS_OK;
}
extern "C" int fusion_c_nuclear_two_body_cm(int ch,double energy,const double*dir,
 fusion_particle_four_vector_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 std::fill(out,out+2,fusion_particle_four_vector_v1{});
 if(!std::isfinite(energy))return PB11_STATUS_INVALID_ARGUMENT;
 if(energy<0)return PB11_STATUS_OUT_OF_RANGE;
 fusion_nuclear_channel_v1 reaction{};
 int status=fusion_c_nuclear_channel(ch,&reaction);if(status)return status;
 if(reaction.product_count!=2)return PB11_STATUS_OUT_OF_RANGE;
 R available=R(energy)+reaction.q_J;
 if(!std::isfinite(available)||available>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
 return fusion_c_two_body_cm(canonical_mass(reaction.product_ids[0]),
    canonical_mass(reaction.product_ids[1]),static_cast<double>(available),dir,out);
}
